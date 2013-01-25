/*
 *  FlowStep_LSRK3MPI.h
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/28/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <limits>

//#include <BlockProcessingMPI.h>
#include <BlockLabMPI.h>
#include <Histogram.h>

#include <FlowStep_LSRK3.h>
#include <Convection_CPP.h>
#include <Diffusion_CPP.h>
#include <Update.h>

#ifdef _SSE_
#include <Convection_SSE.h>
#include <Diffusion_SSE.h>
#endif
#ifdef _AVX_
#include <Convection_AVX.h>
#endif


typedef BlockLabMPI<Lab> LabMPI;

//profile information
namespace LSRK3MPIdata
{
    double t_fs = 0, t_up = 0;
    double t_synch_fs = 0, t_bp_fs = 0;
    int counter = 0, GSYNCH = 0, nsynch = 0;
    
    Histogram histogram;
    
    template<typename Kflow, typename Kupdate, typename Diffusion_CPP>
    void notify(double avg_time_rhs, double avg_time_update, const size_t NBLOCKS, const size_t NTIMES)
    {
		if(LSRK3data::step_id % LSRK3data::ReportFreq == 0 && LSRK3data::step_id > 0)
            histogram.consolidate();
        
		histogram.notify("FLOWSTEP", (float)avg_time_rhs);
		histogram.notify("UPDATE", (float)avg_time_update);
		histogram.notify("STEPID", (float)LSRK3data::step_id);
		
		histogram.notify("NSYNCH", (float)nsynch/NTIMES);
		nsynch = 0;
        
		if(LSRK3data::step_id % LSRK3data::ReportFreq == 0 && LSRK3data::step_id > 0)
		{
			double global_t_synch_fs = 0, global_t_bp_fs = 0;
			double global_avg_time_rhs = 0, global_avg_time_update = 0;
			double global_t_fs, global_t_up;
			
			int global_counter = 0;
			
			MPI::COMM_WORLD.Allreduce(&t_synch_fs, &global_t_synch_fs, 1, MPI::DOUBLE, MPI::SUM);
			MPI::COMM_WORLD.Allreduce(&t_bp_fs, &global_t_bp_fs, 1, MPI::DOUBLE, MPI::SUM);
			MPI::COMM_WORLD.Allreduce(&counter, &global_counter, 1, MPI::INT, MPI::SUM);
			MPI::COMM_WORLD.Allreduce(&avg_time_rhs, &global_avg_time_rhs, 1, MPI::DOUBLE, MPI::SUM);
			MPI::COMM_WORLD.Allreduce(&avg_time_update, &global_avg_time_update, 1, MPI::DOUBLE, MPI::SUM);
			MPI::COMM_WORLD.Allreduce(&t_fs, &global_t_fs, 1, MPI::DOUBLE, MPI::SUM);
			MPI::COMM_WORLD.Allreduce(&t_up, &global_t_up, 1, MPI::DOUBLE, MPI::SUM);
			
			t_synch_fs = t_bp_fs = t_fs = t_up = counter = 0;
			
			global_t_synch_fs /= NTIMES;
			global_t_bp_fs /= NTIMES;
			global_counter /= NTIMES;
			global_t_fs /= NTIMES;
			global_t_up /= NTIMES;
			
			const size_t NRANKS = MPI::COMM_WORLD.Get_size();
			
			if (LSRK3data::verbosity >= 1)
			{
				cout << "FLOWSTEP: " << avg_time_rhs << " s (per substep), " << avg_time_rhs/NBLOCKS*1e3 << " ms (per block) " << global_avg_time_rhs/NBLOCKS*1e3/NRANKS << " ms (per block per node)" << endl;
				
				cout << "TIME LOCALLY AVERAGED FLOWSTEP: " << global_avg_time_rhs/NRANKS <<" s (per substep per node), " << global_avg_time_rhs/NBLOCKS*1e3/NRANKS << " ms (per substep per node per block)" << endl;
				
				cout << "TIME GLOBALLY AVERAGED FLOWSTEP: " << global_t_fs/NRANKS/(double)LSRK3data::ReportFreq << " s (per substep)" << endl;
				
				cout << "===========================STAGE===========================" << endl;
				cout << "Synch done in "<< global_counter/NRANKS/(double)LSRK3data::ReportFreq << " passes" << endl;
				cout << "SYNCHRONIZER FLOWSTEP "<< global_t_synch_fs/NRANKS/(double)LSRK3data::ReportFreq << " s" << endl;
				cout << "BP FLOWSTEP "<< global_t_bp_fs/NRANKS/(double)LSRK3data::ReportFreq << " s" << endl;
				cout << "======================================================" << endl;
				
				Kflow::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, LSRK3data::TLP, NBLOCKS*NRANKS, global_t_fs/(double)LSRK3data::ReportFreq/NRANKS);
				Kupdate::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, LSRK3data::TLP, NBLOCKS*NRANKS, global_t_up/(double)LSRK3data::ReportFreq/NRANKS);
			}
		}
	}
}

template<typename Lab, typename Operator, typename TGrid>
void _process(vector<BlockInfo>& vInfo, Operator rhs, TGrid& grid, const Real t=0, bool tensorial=false)
{
#pragma omp parallel
    {
        vector<BlockInfo>myInfo = vInfo;
        BlockInfo * ary = &myInfo.front();
        
        Operator myrhs = rhs;
        Lab mylab;
        
        const SynchronizerMPI& synch = grid.get_SynchronizerMPI(myrhs);
        
        const int N = vInfo.size();
        mylab.prepare(grid, synch);
        
#pragma omp for schedule(runtime)
        for(int i=0; i<N; i++)
        {
            mylab.load(ary[i], t);
            myrhs(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
        }
    }
}

template<typename TGrid>
class FlowStep_LSRK3MPI : public FlowStep_LSRK3
{
    TGrid & grid;
    
	Real _computeSOS()
	{
		double maxSOS;
        
		const double local_maxSOS = FlowStep_LSRK3::_computeSOS();
		
		MPI::COMM_WORLD.Allreduce(&local_maxSOS, &maxSOS, 1, MPI::DOUBLE, MPI::MAX);
		
		return maxSOS;
	}
	
	template<typename Kflow, typename Kupdate>
	struct LSRKstepMPI
	{
		LSRKstepMPI(TGrid& grid, Real dtinvh, const Real current_time)
		{
			vector<BlockInfo> vInfo = grid.getBlocksInfo();
			
            vector< pair<double, double> > timings;
            
			timings.push_back(step(grid, vInfo, 0      , 1./4, dtinvh, current_time));
			timings.push_back(step(grid, vInfo, -17./32, 8./9, dtinvh, current_time));
			timings.push_back(step(grid, vInfo, -32./27, 3./4, dtinvh, current_time));
            
			double avg1 = ( timings[0].first  + timings[1].first  + timings[2].first  )/3;
			double avg2 = ( timings[0].second + timings[1].second + timings[2].second )/3;
            
            LSRK3MPIdata::notify<Kflow, Kupdate>(avg1, avg2, vInfo.size(), 3);
		}
		
		pair<double, double> step(TGrid& grid, vector<BlockInfo>& vInfo, Real a, Real b, Real dtinvh, const Real current_time)
		{
			Timer timer;
            LSRK3data::FlowStep<Kflow, Lab> rhs(a, dtinvh);
            
            timer.start();
			
            SynchronizerMPI& synch = ((TGrid&)grid).sync(rhs);
            
			while (!synch.done())
			{
				Timer timer2;
                
				timer2.start();
				vector<BlockInfo> avail = synch.avail(LSRK3MPIdata::GSYNCH);
				LSRK3MPIdata::t_synch_fs += timer2.stop();
				
				timer2.start();
				_process< LabMPI >(avail, rhs, (TGrid&)grid, current_time);
				LSRK3MPIdata::t_bp_fs += timer2.stop();
				
				LSRK3MPIdata::counter++;
				LSRK3MPIdata::nsynch++;
			}
            
            const double totalRHS = timer.stop();
			
            if(LSRK3data::nu1>0)
            {
                LSRK3data::Diffusion<Kdiff, Lab> diffusion(dtinvh, LSRK3data::nu1, LSRK3data::nu2);
				
                SynchronizerMPI& synch = ((TGrid&)grid).sync(diffusion);
                
				while (!synch.done())
				{
                    vector<BlockInfo> avail = synch.avail(LSRK3MPIdata::GSYNCH);
                    
                    _process< LabMPI >(avail, diffusion, (TGrid&)grid, current_time);
                }
            }
            
			LSRK3data::Update<Kupdate> update(b, &vInfo.front());
			timer.start();
			update.omp(vInfo.size());
			const double totalUPDATE = timer.stop();
			
			LSRK3MPIdata::t_fs += totalRHS;
			LSRK3MPIdata::t_up += totalUPDATE;
            
			return pair<double, double>(totalRHS,totalUPDATE);
		}
	};
	
public:
    
	FlowStep_LSRK3MPI(TGrid & grid, const Real CFL, const Real gamma1, const Real gamma2, ArgumentParser& parser, const int verbosity, Profiler* profiler=NULL, const Real pc1=0, const Real pc2=0):
    FlowStep_LSRK3(grid, CFL, gamma1, gamma2, parser, verbosity, profiler, pc1, pc2), grid(grid)
    {
        if (verbosity>=1) cout << "GSYNCH " << parser("-gsync").asInt(LSRK3data::TLP) << endl;
    }
	
	Real operator()(const Real max_dt)
	{
		set_constants();
        
        //here we just check stuff and compute the next dt
        if (verbosity>=1 && LSRK3data::step_id==0)
        {
            cout << "Grid spacing and smoothing length are: " << h << ", " << smoothlength << endl;
            }
            
            LSRK3MPIdata::GSYNCH = parser("-gsync").asInt(LSRK3data::TLP);
            
            const Real maxSOS = _computeSOS();
            double dt = min(max_dt, CFL*h/maxSOS);
            
            if (MPI::COMM_WORLD.Get_rank()==0)
            {
                cout << "sos max is " << maxSOS << ", " << "advection dt is "<< dt << "\n";
            }
            
            if (LSRK3data::nu1>0)
            dt = min(dt, (Real)(h*h/(12*max(LSRK3data::nu1,LSRK3data::nu2))) );
            
            if (verbosity>=1)
            {
                cout << "dt is "<< dt << "\n";
                cout << "Dispatcher is " << LSRK3data::dispatcher << endl;
            }
            
            if (maxSOS>1e6)
            {
                cout << "Speed of sound is too high. Is it realistic?" << endl;
                MPI::COMM_WORLD.Abort(1);
            }
            
            if (dt<std::numeric_limits<double>::epsilon() * 1e1)
            {
                cout << "Last time step encountered." << endl;
                
                return 0;
            }
            
            //now we perform an entire RK step
            if (parser("-kernels").asString("cpp")=="cpp")
            {
                LSRKstepMPI<Convection_CPP, Update_CPP, Diffusion_CPP>(grid, dt/h, current_time);
            }
#ifdef _SSE_
            else if (parser("-kernels").asString("cpp")=="sse")
            {
                LSRKstepMPI<Convection_SSE, Update_SSE, Diffusion_SSE>(grid, dt/h, current_time);
            }
#endif
#ifdef _AVX_
            else if (parser("-kernels").asString("cpp")=="avx")
            {
                LSRKstepMPI<Convection_AVX, Update_AVX>(grid, dt/h, current_time);
            }
#endif
            else
            {
                cout << "combination not supported yet" << endl;
                MPI::COMM_WORLD.Abort(1);
            }
            
            LSRK3data::step_id++; current_time+=dt;
            
            return dt;
            }
            };
