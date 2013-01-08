/*
 *  FlowStep_LSRK3.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include <vector>

#include <utility>
#include <iostream>
#include <limits>

#include <Timer.h>
#include <Profiler.h>
#include <Convection_CPP.h>

#ifdef _SSE_
#include <emmintrin.h>
#include <Convection_SSE.h>
#endif

#ifdef _AVX_
#include <Convection_AVX.h>
#endif

#include <Update.h>
#include <MaxSpeedOfSound.h>

using namespace std;

#include "FlowStep_LSRK3.h"
#include "Tests.h"

namespace LSRK3data 
{
	Real gamma1 = -1;
	Real gamma2 = -1;
	Real smoothlength = -1;
	Real pc1 = 0;
    Real pc2 = 0;
	
	int verbosity;
	int  NCORES, TLP;
	float PEAKPERF_CORE, PEAKBAND;
	
	string dispatcher;
	
	int step_id = 0;
    int ReportFreq = 1;
}

template<typename Lab, typename Kernel>
void _process(const Real a, const Real dtinvh, vector<BlockInfo>& myInfo, FluidGrid& grid, const Real t=0, bool tensorial=false) 
{
	const int stencil_start[3] = {-3,-3,-3};
	const int stencil_end[3] = {4,4,4};
	BlockInfo * ary = &myInfo.front();
	const int N = myInfo.size();
    
#pragma omp parallel
	{
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif
		Kernel kernel(a, dtinvh);
		
		Lab mylab;
		mylab.prepare(grid, stencil_start, stencil_end, tensorial);
                
#pragma omp for schedule(runtime)
		for(int i=0; i<N; i++) 
		{
			mylab.load(ary[i], t);
            
            const Real * const srcfirst = &mylab(-3,-3,-3).rho;
            const int labSizeRow = mylab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*mylab.template getActualSize<1>();	
			
			Real * const destfirst = &((FluidBlock*)ary[i].ptrBlock)->tmp[0][0][0][0];

			kernel.compute(srcfirst, FluidBlock::gptfloats, labSizeRow, labSizeSlice, 
						   destfirst, FluidBlock::gptfloats, FluidBlock::sizeX, FluidBlock::sizeX*FluidBlock::sizeY);
		}
	}
}

template < typename TSOS>
Real _computeSOS_OMP(FluidGrid& grid,  bool bAwk)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    const size_t N = vInfo.size();
    const BlockInfo * const ary = &vInfo.front();
  //Real * const  local_sos = (Real *)_mm_malloc(N*sizeof(Real), 16);
    Real * tmp = NULL;
    int error = posix_memalign((void**)&tmp, std::max(8, _ALIGNBYTES_), sizeof(Real) * N);
    assert(error == 0);
    
    Real * const local_sos = tmp;
    
#pragma omp parallel
    {
        
#ifdef _USE_NUMA_
		const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif
        
        TSOS kernel;
#pragma omp for schedule(runtime)
        for (size_t i=0; i<N; ++i)
        {
            FluidBlock & block = *(FluidBlock *)ary[i].ptrBlock;
            local_sos[i] =  kernel.compute(&block.data[0][0][0].rho, FluidBlock::gptfloats);
        }
    }
    
    static const int CHUNKSIZE = 64*1024/sizeof(Real);
    Real global_sos = local_sos[0];
    /*
    #pragma omp parallel
    {
		Real mymax = local_sos[0];
		
		#pragma omp for
		for(int i=0; i<N; ++i)
			mymax = max(local_sos[i], mymax);
			
		#pragma omp critical
		{
			global_sos = max(global_sos, mymax);
		}
	}*/
    
#pragma omp parallel for schedule(static)
    for(size_t i=0; i<N; i+= CHUNKSIZE)
    {
        Real mymax = local_sos[i];
        const size_t e = min(N, i+CHUNKSIZE);
        for(size_t c=i; c<e; ++c)
            mymax = max(mymax, local_sos[c]);
        
#pragma omp critical
        {
            global_sos = max(global_sos, mymax);
        }
    }		
    
    free(tmp);
    
     return global_sos;
}

Real FlowStep_LSRK3::_computeSOS(bool bAwk)
{
    Real sos = -1;
	
    const string kernels = parser("-kernels").asString("cpp");
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
	Timer timer;
    
	timer.start();
    sos = _computeSOS_OMP<MaxSpeedOfSound_CPP>(grid,  bAwk);
    const Real time = timer.stop();
    
    if (LSRK3data::verbosity >= 1 && LSRK3data::step_id % LSRK3data::ReportFreq == 0)
    {
        MaxSpeedOfSound_CPP::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, LSRK3data::TLP, vInfo.size(), time, bAwk);
        
        cout << "MAXSOS: " << time << "s (per substep), " << time/vInfo.size()*1e3 << " ms (per block)" << endl;
    }
    
    return sos;
}

template<typename Kflow, typename Kupdate>
struct LSRKstep
{	
    LSRKstep(FluidGrid& grid, Real dtinvh, const Real current_time, bool bAwk)
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        
        vector< vector<double> > timings;
        
        timings.push_back(step(grid, vInfo, 0      , 1./4, dtinvh, current_time));
        timings.push_back(step(grid, vInfo, -17./32, 8./9, dtinvh, current_time));
        timings.push_back(step(grid, vInfo, -32./27, 3./4, dtinvh, current_time));
        
        const double avg1 = ( timings[0][0] + timings[1][0] + timings[2][0] )/3;
        const double avg2 = ( timings[0][1] + timings[1][1] + timings[2][1] )/3;
        
        if (LSRK3data::verbosity >= 1 && LSRK3data::step_id % LSRK3data::ReportFreq == 0)
        {
            cout << "FLOWSTEP: " << avg1 << "s (per substep), " << avg1/vInfo.size()*1e3 << " ms (per block)" << endl;
            cout << "UPDATE: " << avg2 << "s (per substep), " << avg2/vInfo.size()*1e3 << " ms (per block)" << endl;
            
            Kflow::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg1);
            Kupdate::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg2);    
		 }
    }
    
    vector<double> step(FluidGrid& grid, vector<BlockInfo>& vInfo, Real a, Real b, Real dtinvh, const Real current_time)
    {
        Timer timer;
        vector<double> res;
        
        LSRK3data::FlowStep<Kflow, Lab> rhs(a, dtinvh);
       
        timer.start();     
        _process<Lab, Kflow>(a, dtinvh, vInfo, grid, current_time);
        const double t1 = timer.stop();
        
        LSRK3data::Update<Kupdate> update(b, &vInfo.front());
        
        timer.start();
        update.omp(vInfo.size());
        const double t2 = timer.stop();
        
		res.push_back(t1);
		res.push_back(t2);
		
		return res;
    }
};

void FlowStep_LSRK3::set_constants()
{
    LSRK3data::gamma1 = gamma1;
    LSRK3data::gamma2 = gamma2;
    LSRK3data::smoothlength = smoothlength;
    LSRK3data::verbosity = verbosity;
    LSRK3data::PEAKBAND = PEAKBAND;
    LSRK3data::PEAKPERF_CORE = PEAKPERF_CORE;
    LSRK3data::NCORES = parser("-ncores").asInt(1);
    LSRK3data::TLP = parser("-nthreads").asInt(1);
    LSRK3data::pc1 = pc1;
    LSRK3data::pc2 = pc2;
    LSRK3data::dispatcher = blockdispatcher;
    LSRK3data::ReportFreq = parser("-report").asInt(20);
 }

Real FlowStep_LSRK3::operator()(const Real max_dt)
{
    set_constants();
    
    const Real maxSOS = _computeSOS(bAwk);
    double dt = min(max_dt, CFL*h/maxSOS);
    cout << "sos max is " << setprecision(8) << maxSOS << ", " << "dt is "<< dt << "\n";
    
    if (maxSOS > 1e6)
    {
        cout << "Speed of sound is too high. Is it realistic?" << endl;
        abort();
    }
    
    if (dt<std::numeric_limits<double>::epsilon()*1e1)
    {
        cout << "Last time step encountered." << endl;
        return 0;
    }

    if (LSRK3data::verbosity >= 1)
        cout << "Dispatcher is " << LSRK3data::dispatcher << endl;
    
    if (parser("-kernels").asString("cpp")=="cpp")
        LSRKstep<Convection_CPP, Update_CPP>(grid, dt/h, current_time, bAwk);
#ifdef _SSE_
    else if (parser("-kernels").asString("cpp")=="sse")
        LSRKstep<Convection_SSE, Update_CPP>(grid, dt/h, current_time, bAwk);
#endif
#ifdef _AVX_
    else if (parser("-kernels").asString("cpp")=="avx")
        LSRKstep<Convection_AVX, Update_AVX>(grid, dt/h, current_time, bAwk);
#endif
    else
    {
        cout << "combination not supported yet" << endl;
        abort();
    }
    
    LSRK3data::step_id++;
    
    return dt;
}
