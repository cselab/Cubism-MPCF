/*
 *  FlowStep_LSRK3.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include <vector>

#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>

#include <utility>
#include <iostream>

#include <emmintrin.h>

#include <Timer.h>
#include <Profiler.h>
#include <Convection_CPP.h>

#ifdef _SSE_
#include <Convection_SSE.h>
#include <SurfaceTension_SSE.h>
#include <Diffusion_SSE.h>
#endif

#ifdef _AVX_
#include <Convection_AVX.h>
#include <SurfaceTension_AVX.h>
#include <Diffusion_AVX.h>
#endif

#include <Update.h>
#include <MaxSpeedOfSound.h>
#include <SurfaceTension_CPP.h>
#include <Diffusion_CPP.h>

using namespace tbb;
using namespace std;

#include "FlowStep_LSRK3.h"
#include "Tests.h"

namespace LSRK3data {
	Real gamma1 = -1;
	Real gamma2 = -1;
	Real smoothlength = -1;
	int verbosity;
	Profiler * profiler = NULL;
	float PEAKPERF_CORE, PEAKBAND;
	int  NCORES, TLP;
    Real pc1 = 0;
    Real pc2 = 0;
	string dispatcher;
	bool bAffinity = false;
	
	tbb::affinity_partitioner affinitypart;
    int step_id = 0;
    int ReportFreq = 1;
    Real sten_sigma = 0;
    Real nu1, nu2 = 0;
}

template < typename Kernel >
struct MaxSOS
{
	Real sos;
	BlockInfo * ary;
	
	MaxSOS(BlockInfo * ary): ary(ary), sos(0) { }
	
	MaxSOS(const MaxSOS& c, split): ary(c.ary), sos(0) { }
	
	void join(const MaxSOS& c)
	{
		sos = max(sos, c.sos);
	}
	
	void operator()(blocked_range<int> range)
	{
		Kernel kernel;
		
		for(int r=range.begin(); r<range.end(); ++r)
		{
			FluidBlock & block = *(FluidBlock *)ary[r].ptrBlock;
			sos = max(sos, kernel.compute(&block.data[0][0][0].rho, FluidBlock::gptfloats));
		}
	}
};

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
                
#pragma omp for schedule(static)
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

template<typename Lab, typename Kernel>
void _process_surface_tension(const Real sigma, const Real dtinvh, vector<BlockInfo>& myInfo, FluidGrid& grid, const Real t=0, bool tensorial=true) 
{
    const int stencil_start[3] = {-1,-1,-1};
    const int stencil_end[3] = {2,2,2};
    BlockInfo * ary = &myInfo.front();
    const int N = myInfo.size();
    
#pragma omp parallel
    {
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif
        Kernel kernel(1, dtinvh, max((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), min((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), myInfo.front().h_gridpoint, LSRK3data::smoothlength, sigma);
        
        Lab mylab;
        mylab.prepare(grid, stencil_start, stencil_end, tensorial);
        
#pragma omp for schedule(static)
        for(int i=0; i<N; i++) {
            mylab.load(ary[i], t);
            
			const Real * const srcfirst = &mylab(-1,-1,-1).rho;
			const int labSizeRow = mylab.template getActualSize<0>();
			const int labSizeSlice = labSizeRow*mylab.template getActualSize<1>();
			
			Real * const destfirst = &((FluidBlock*)ary[i].ptrBlock)->tmp[0][0][0][0];
            
            kernel.compute(srcfirst, FluidBlock::gptfloats, labSizeRow, labSizeSlice, 
                           destfirst, FluidBlock::gptfloats, FluidBlock::sizeX, FluidBlock::sizeX*FluidBlock::sizeY);
        }
    }
}

template<typename Lab, typename Kernel>
void _process_diffusion(const Real dtinvh, vector<BlockInfo>& myInfo, FluidGrid& grid, const Real t=0, bool tensorial=true)
{
    const int stencil_start[3] = {-1,-1,-1};
    const int stencil_end[3] = {2,2,2};
    BlockInfo * ary = &myInfo.front();
    const int N = myInfo.size();
    
#pragma omp parallel
    {
#ifdef _USE_NUMA_
		const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
		const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif
		Kernel kernel(dtinvh, LSRK3data::nu1, LSRK3data::nu2, max((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), min((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), myInfo[0].h_gridpoint, LSRK3data::smoothlength, dtinvh);
        
        Lab mylab;
        mylab.prepare(grid, stencil_start, stencil_end, tensorial);
        
#pragma omp for schedule(static)
        for(int i=0; i<N; i++) 
		{
            mylab.load(ary[i], t);
            
			const Real * const srcfirst = &mylab(stencil_start[0],stencil_start[1],stencil_start[2]).rho;
			const int labSizeRow = mylab.template getActualSize<0>();
			const int labSizeSlice = labSizeRow*mylab.template getActualSize<1>();
			
			Real * const destfirst = &((FluidBlock*)ary[i].ptrBlock)->tmp[0][0][0][0];
            
            kernel.compute(srcfirst, FluidBlock::gptfloats, labSizeRow, labSizeSlice, 
                           destfirst, FluidBlock::gptfloats, FluidBlock::sizeX, FluidBlock::sizeX*FluidBlock::sizeY);
        }
    }
}

template<typename Lab, typename Operator>
struct TBBWorker {
    const vector<BlockInfo>myInfo;
    Operator myrhs;
    FluidGrid* grid;
    Real t;
    const bool tensorial;
    
    TBBWorker(vector<BlockInfo>vInfo, Operator rhs, FluidGrid& grid, const Real t, const bool tensorial): 
    myrhs(rhs), grid(&grid), myInfo(vInfo), t(t), tensorial(tensorial){}
    TBBWorker(const TBBWorker& c): myrhs(c.myrhs), grid(c.grid), myInfo(c.myInfo), t(c.t), tensorial(c.tensorial){}
    
    void operator()(blocked_range<int> range) const {
        Lab mylab;
        mylab.prepare(*grid, myrhs.stencil_start, myrhs.stencil_end, tensorial);
        
        const BlockInfo * ary = &myInfo.front();
        for(int i=range.begin(); i<range.end(); i++){
            mylab.load(ary[i], t);
            myrhs(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
        }
    }
    
    static void _process(vector<BlockInfo>& vInfo, Operator rhs, FluidGrid& grid, const Real t=0,  const bool tensorial=false) {
        if (!LSRK3data::bAffinity)
            tbb::parallel_for(blocked_range<int>(0, vInfo.size()), TBBWorker(vInfo, rhs, grid, t, tensorial), auto_partitioner() );
        else
            tbb::parallel_for(blocked_range<int>(0, vInfo.size()), TBBWorker(vInfo, rhs, grid, t, tensorial), LSRK3data::affinitypart );
    }
};

template < typename TSOS>
Real _computeSOS_OMP(FluidGrid& grid,  bool bAwk)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    const size_t N = vInfo.size();
    const BlockInfo * const ary = &vInfo.front();
    Real * const  local_sos = (Real *)_mm_malloc(N*sizeof(Real), 16);
    
#pragma omp parallel
    {
        
#ifdef _USE_NUMA_
		const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif
        
        TSOS kernel;
#pragma omp for schedule(static)
        for (size_t i=0; i<N; ++i)
        {
            FluidBlock & block = *(FluidBlock *)ary[i].ptrBlock;
            local_sos[i] =  kernel.compute(&block.data[0][0][0].rho, FluidBlock::gptfloats);
        }
    }
    
    static const int CHUNKSIZE = 64*1024/sizeof(Real);
    Real global_sos = local_sos[0];
    
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
    
    _mm_free(local_sos);
    
    return global_sos;
}

Real FlowStep_LSRK3::_computeSOS(bool bAwk)
{
    if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("MAXSOS");
    
    Real sos = -1;
	
    const string kernels = parser("-kernels").asString("cpp");
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
	Timer timer;
    
	timer.start();
	
    if (LSRK3data::dispatcher == "tbblight")
    {
        if (kernels=="sse" || kernels=="avx")
        {
#if defined(_SSE_) || defined(_AVX_)
            MaxSOS<MaxSpeedOfSound_SSE> compute_sos(&vInfo.front());
            if (!LSRK3data::bAffinity)
                parallel_reduce(blocked_range<int>(0, vInfo.size()), compute_sos, auto_partitioner());
            else
                parallel_reduce(blocked_range<int>(0, vInfo.size()), compute_sos, LSRK3data::affinitypart);
            
            sos = compute_sos.sos;
#endif
        }
        else
        {
            MaxSOS<MaxSpeedOfSound_CPP> compute_sos(&vInfo.front());
            if (!LSRK3data::bAffinity)
                parallel_reduce(blocked_range<int>(0, vInfo.size()), compute_sos, auto_partitioner());
            else
                parallel_reduce(blocked_range<int>(0, vInfo.size()), compute_sos, LSRK3data::affinitypart);
            
            sos = compute_sos.sos;
        }
    }
    else if (kernels=="sse" || kernels=="avx")
    {
#if defined(_SSE_) || defined(_AVX_)
        sos = _computeSOS_OMP<MaxSpeedOfSound_SSE>(grid,  bAwk);
#endif
    }
    else
        sos = _computeSOS_OMP<MaxSpeedOfSound_CPP>(grid,  bAwk);
    
    const Real time = timer.stop();
    
    if (LSRK3data::verbosity >= 1 && LSRK3data::step_id % LSRK3data::ReportFreq == 0)
    {
        if (kernels=="sse" || kernels=="avx")
        {
#if defined(_SSE_) || defined(_AVX_)
            MaxSpeedOfSound_SSE::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, LSRK3data::TLP, vInfo.size(), time, bAwk);
#endif
        }
        else
            MaxSpeedOfSound_CPP::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, LSRK3data::TLP, vInfo.size(), time, bAwk);
        
        cout << "MAXSOS: " << time << "s (per substep), " << time/vInfo.size()*1e3 << " ms (per block)" << endl;
    }
    if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();
    
    return sos;
}

template<typename Kflow, typename Kupdate, typename Ksten, typename Kdiff>
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
        const double avg3 = ( timings[0][2] + timings[1][2] + timings[2][2] )/3;
		const double avg4 = ( timings[0][3] + timings[1][3] + timings[2][3] )/3;
        
        if (LSRK3data::verbosity >= 1 && LSRK3data::step_id % LSRK3data::ReportFreq == 0)
        {
            cout << "FLOWSTEP: " << avg1 << "s (per substep), " << avg1/vInfo.size()*1e3 << " ms (per block)" << endl;
            cout << "UPDATE: " << avg2 << "s (per substep), " << avg2/vInfo.size()*1e3 << " ms (per block)" << endl;
            cout << "DIFFUSION: " << avg3 << "s (per substep), " << avg3/vInfo.size()*1e3 << " ms (per block)" << endl;
			cout << "SUTENSION: " << avg4 << "s (per substep), " << avg4/vInfo.size()*1e3 << " ms (per block)" << endl;
            
            Kflow::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg1);
            Kupdate::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg2);
            
			if(LSRK3data::sten_sigma!=0)
		    {
				Ksten sten(1, dtinvh, max((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), min((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), vInfo.front().h_gridpoint, LSRK3data::smoothlength, LSRK3data::sten_sigma);
				sten.printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg4, bAwk);
			}
            
			if(LSRK3data::nu1!=0)
			{
				Kdiff diffusion(dtinvh, LSRK3data::nu1, LSRK3data::nu2, max((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), min((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), vInfo.front().h_gridpoint, LSRK3data::smoothlength, dtinvh);
                
				diffusion.printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size(), avg3, bAwk);
			}		
        }
    }
    
    vector<double> step(FluidGrid& grid, vector<BlockInfo>& vInfo, Real a, Real b, Real dtinvh, const Real current_time)
    {
        Timer timer;
        vector<double> res;
        
        if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("FLOWSTEP");
        LSRK3data::FlowStep<Kflow, Lab> rhs(a, dtinvh);
        timer.start();
        
        if (LSRK3data::verbosity >= 1)
            cout << "Dispatcher is " << LSRK3data::dispatcher << endl;
        
        if (LSRK3data::dispatcher == "omp")
            _process<Lab, Kflow>(a, dtinvh, vInfo, grid, current_time);
        else if (LSRK3data::dispatcher == "tbblight")
            TBBWorker<Lab, LSRK3data::FlowStep<Kflow, Lab> >::_process(vInfo, rhs, grid, current_time);
        else
            BlockProcessing::process<Lab>(vInfo, rhs, grid, current_time); 
        
        const double t1 = timer.stop();
        if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();
                
		double tsurface_tension = 0;
        if(LSRK3data::sten_sigma!=0)
        {
            LSRK3data::SurfaceTension<Ksten, Lab> sten(dtinvh, LSRK3data::sten_sigma);
            
            if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("SURFACE TENSION");
            
			timer.start();
            
            if (LSRK3data::dispatcher == "omp")
                _process_surface_tension<Lab, Ksten>(LSRK3data::sten_sigma, dtinvh, vInfo, grid, current_time, true);
            else if (LSRK3data::dispatcher == "tbblight")
                TBBWorker<Lab, LSRK3data::SurfaceTension<Ksten, Lab> >::_process(vInfo, sten, grid, current_time, true);
            else               
                BlockProcessing::process<Lab>(vInfo, sten, grid, current_time, true);            
            
			tsurface_tension = timer.stop();
            
            if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();
        }
		
        
		double t3 = 0.;
        if(LSRK3data::nu1!=0)
        {
            if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("DIFFUSION");
            LSRK3data::Diffusion<Kdiff, Lab> diffusion(dtinvh, LSRK3data::nu1, LSRK3data::nu2);
			timer.start();
			if (LSRK3data::dispatcher == "omp")
				_process_diffusion<Lab, Kdiff>(dtinvh, vInfo, grid, current_time, true);
			else if (LSRK3data::dispatcher == "tbblight")
                TBBWorker<Lab, LSRK3data::Diffusion<Kdiff, Lab> >::_process(vInfo, diffusion, grid, current_time, true);
            else 
				BlockProcessing::process<Lab>(vInfo, diffusion, grid, current_time, true); 
			
			t3 = timer.stop();
            if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();    
        }
        
        if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("UPDATE");
        LSRK3data::Update<Kupdate> update(b, &vInfo.front());
        timer.start();
        if (LSRK3data::dispatcher == "omp")
            update.omp(vInfo.size());
        else 
        {
            if (!LSRK3data::bAffinity)
                parallel_for(blocked_range<int>(0, vInfo.size()), update, auto_partitioner());
            else 
                parallel_for(blocked_range<int>(0, vInfo.size()), update, LSRK3data::affinitypart);
        }
        
        const double t2 = timer.stop();
        if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();
        
		res.push_back(t1);
		res.push_back(t2);
		res.push_back(t3);
		res.push_back(tsurface_tension);
        return res;
    }
};

void FlowStep_LSRK3::set_constants()
{
    LSRK3data::gamma1 = gamma1;
    LSRK3data::gamma2 = gamma2;
    LSRK3data::smoothlength = smoothlength;
    LSRK3data::verbosity = verbosity;
    LSRK3data::profiler = profiler;
    LSRK3data::PEAKBAND = PEAKBAND;
    LSRK3data::PEAKPERF_CORE = PEAKPERF_CORE;
    LSRK3data::NCORES = parser("-ncores").asInt(1);
    LSRK3data::TLP = parser("-nthreads").asInt(1);
    LSRK3data::pc1 = pc1;
    LSRK3data::pc2 = pc2;
    LSRK3data::dispatcher = blockdispatcher;
    LSRK3data::bAffinity = parser("-affinity").asBool(false);
    LSRK3data::ReportFreq = parser("-report").asInt(20);
    LSRK3data::sten_sigma = parser("-sten").asDouble(0);
    LSRK3data::nu1    = parser("-nu1").asDouble(0);
    LSRK3data::nu2    = parser("-nu2").asDouble(LSRK3data::nu1);
}

Real FlowStep_LSRK3::operator()(const Real max_dt)
{
    set_constants();
    
    const Real maxSOS = _computeSOS(bAwk);
    Real dt = min(max_dt, CFL*h/maxSOS);
    cout << "sos max is " << setprecision(8) << maxSOS << ", " << "dt is "<< dt << "\n";
    
    if (LSRK3data::sten_sigma>0) 
    {
        const Real sumrho = parser("-sumrho").asDouble(HUGE_VAL);
        dt = min(dt, (Real)sqrt(sumrho*h*h*h/(4*M_PI*LSRK3data::sten_sigma)));
    }
    
    if (LSRK3data::nu1>0)
        dt = min(dt, (Real)(h*h/(12*max(LSRK3data::nu1,LSRK3data::nu2))) );
    
    if (dt<1e-10 || maxSOS<1e-10)
    {
        cout << "Something went wrong. Zero time step detected." << endl;
        abort();
    }
    
    if (parser("-kernels").asString("cpp")=="cpp")
        LSRKstep<Convection_CPP, Update_CPP, SurfaceTension_CPP, Diffusion_CPP>(grid, dt/h, current_time, bAwk);
#ifdef _SSE_
    else if (parser("-kernels").asString("cpp")=="sse")
        LSRKstep<Convection_SSE, Update_SSE, SurfaceTension_SSE, Diffusion_SSE>(grid, dt/h, current_time, bAwk);
#endif
#ifdef _AVX_
    else if (parser("-kernels").asString("cpp")=="avx")
        LSRKstep<Convection_AVX, Update_AVX, SurfaceTension_AVX, Diffusion_AVX>(grid, dt/h, current_time, bAwk);
#endif
    else
    {
        cout << "combination not supported yet" << endl;
        abort();
    }
    
    LSRK3data::step_id++;
    return dt;
}
