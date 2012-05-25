/*
 *  FlowStep_LSRK3MPI.h
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/28/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <tbb/task.h>

#include <BlockProcessingMPI.h>
#include <BlockLabMPI.h>

#include <FlowStep_LSRK3.h>
#include <Convection_CPP.h>
#include <Convection_SSE.h>
#ifdef _AVX_
#include <Convection_AVX.h>
#include <SurfaceTension_AVX.h>
#include <Diffusion_AVX.h>
#endif
#include <Update.h>
#include <SurfaceTension_CPP.h>
#include <SurfaceTension_SSE.h>
#include <Diffusion_CPP.h>
#include <Diffusion_SSE.h>

#include "Histogram.h"

typedef BlockLabMPI<Lab> LabMPI;

namespace LSRK3MPIdata {
    double t_fs=0, t_up=0;
    double t_synch_fs=0, t_bp_fs=0;
    int counter=0, GSYNCH=0, nsynch=0;
    Histogram histogram;
    double t_sten=0, t_diff=0;
}

template<typename Lab, typename Operator, typename TGrid>
void _process(vector<BlockInfo>& vInfo, Operator rhs, TGrid& grid, const Real t=0, bool tensorial=false) {
    
#pragma omp parallel
    {       
        vector<BlockInfo>myInfo = vInfo;
        BlockInfo * ary = &myInfo.front();
        
        Operator myrhs = rhs;
        Lab mylab;
        
        const SynchronizerMPI& Synch = grid.get_SynchronizerMPI(myrhs);
        
        const int N = vInfo.size();
        mylab.prepare(grid, Synch);
        
#pragma omp for schedule(dynamic, 1)
        for(int i=0; i<N; i++) {
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
	
	template<typename Kflow, typename Kupdate, typename Ksten, typename Kdiff>
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
            
			if(LSRK3data::step_id%LSRK3data::ReportFreq==0 && LSRK3data::step_id>0)
				LSRK3MPIdata::histogram.consolidate();
			
            LSRK3MPIdata::histogram.notify("FLOWSTEP", (float)avg1);
            LSRK3MPIdata::histogram.notify("UPDATE", (float)avg2);
            LSRK3MPIdata::histogram.notify("DIFFUSION", (float)LSRK3MPIdata::t_diff/3);
            LSRK3MPIdata::histogram.notify("SURFACETENSION", (float)LSRK3MPIdata::t_sten/3);
			
	    LSRK3MPIdata::histogram.notify("NSYNCH", (float)LSRK3MPIdata::nsynch/3);
	    LSRK3MPIdata::histogram.notify("STEPID", (float)LSRK3data::step_id);
	    LSRK3MPIdata::nsynch = 0;

            if(LSRK3data::step_id%LSRK3data::ReportFreq==0 && LSRK3data::step_id>0)
            {
                double global_t_synch_fs=0, global_t_bp_fs=0;
                double global_t_fs, global_t_up, global_t_diff, global_t_sten;
                double global_avg1=0, global_avg2=0;
                int global_counter=0;
				
                MPI::COMM_WORLD.Allreduce(&LSRK3MPIdata::t_synch_fs, &global_t_synch_fs, 1, MPI::DOUBLE, MPI::SUM);
                MPI::COMM_WORLD.Allreduce(&LSRK3MPIdata::t_bp_fs, &global_t_bp_fs, 1, MPI::DOUBLE, MPI::SUM);
                MPI::COMM_WORLD.Allreduce(&LSRK3MPIdata::counter, &global_counter, 1, MPI::INT, MPI::SUM);
                MPI::COMM_WORLD.Allreduce(&avg1, &global_avg1, 1, MPI::DOUBLE, MPI::SUM);
                MPI::COMM_WORLD.Allreduce(&avg2, &global_avg2, 1, MPI::DOUBLE, MPI::SUM);
                MPI::COMM_WORLD.Allreduce(&LSRK3MPIdata::t_fs, &global_t_fs, 1, MPI::DOUBLE, MPI::SUM);
                MPI::COMM_WORLD.Allreduce(&LSRK3MPIdata::t_up, &global_t_up, 1, MPI::DOUBLE, MPI::SUM);
				MPI::COMM_WORLD.Allreduce(&LSRK3MPIdata::t_diff, &global_t_diff, 1, MPI::DOUBLE, MPI::SUM);
				MPI::COMM_WORLD.Allreduce(&LSRK3MPIdata::t_sten, &global_t_sten, 1, MPI::DOUBLE, MPI::SUM);
				
                global_t_synch_fs/=3.;
                global_t_bp_fs/=3.;
                global_counter/=3.;
                global_t_fs/=3.;
				global_t_up/=3;
				global_t_diff/=3;
				global_t_sten/=3;
                
                if (LSRK3data::verbosity >= 1)
                {
                    cout << "FLOWSTEP: " << avg1 << " s (per substep), " << avg1/vInfo.size()*1e3 << " ms (per block) " << global_avg1/vInfo.size()*1e3/MPI::COMM_WORLD.Get_size() << " ms (per block per node)" << endl;
                    
                    cout << "TIME LOCALLY AVERAGED FLOWSTEP: " << global_avg1/MPI::COMM_WORLD.Get_size() <<" s (per substep per node), " << global_avg1/vInfo.size()*1e3/MPI::COMM_WORLD.Get_size() << " ms (per substep per node per block)" << endl;
                    
                    cout << "TIME GLOBALLY AVERAGED FLOWSTEP: " << global_t_fs/MPI::COMM_WORLD.Get_size()/(double)LSRK3data::ReportFreq << " s (per substep)" << endl;
                    
                    cout << "===========================STAGE===========================" << endl;
                    cout << "Synch done in "<< global_counter/MPI::COMM_WORLD.Get_size()/(double)LSRK3data::ReportFreq << " passes" << endl;
                    cout << "SYNCHRONIZER FLOWSTEP "<< global_t_synch_fs/MPI::COMM_WORLD.Get_size()/(double)LSRK3data::ReportFreq << " s" << endl;
                    cout << "BP FLOWSTEP "<< global_t_bp_fs/MPI::COMM_WORLD.Get_size()/(double)LSRK3data::ReportFreq << " s" << endl;      
                    cout << "======================================================" << endl;
					
                    Kflow::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, LSRK3data::TLP, vInfo.size()*MPI::COMM_WORLD.Get_size(), global_t_fs/(double)LSRK3data::ReportFreq/MPI::COMM_WORLD.Get_size());
					Kupdate::printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, LSRK3data::TLP, vInfo.size()*MPI::COMM_WORLD.Get_size(), global_t_up/(double)LSRK3data::ReportFreq/MPI::COMM_WORLD.Get_size());
					
					if(LSRK3data::sten_sigma!=0)
					{
						Ksten sten(1, dtinvh, max((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), min((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), vInfo.front().h_gridpoint, LSRK3data::smoothlength, LSRK3data::sten_sigma);
						sten.printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size()*MPI::COMM_WORLD.Get_size(), global_t_sten/(double)LSRK3data::ReportFreq/MPI::COMM_WORLD.Get_size());
					}
					
					if(LSRK3data::nu1!=0)
					{
						Kdiff diffusion(dtinvh, LSRK3data::nu1, LSRK3data::nu2, max((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), min((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), vInfo.front().h_gridpoint, LSRK3data::smoothlength, dtinvh);
						diffusion.printflops(LSRK3data::PEAKPERF_CORE*1e9, LSRK3data::PEAKBAND*1e9, LSRK3data::NCORES, 1, vInfo.size()*MPI::COMM_WORLD.Get_size(), global_t_diff/(double)LSRK3data::ReportFreq/MPI::COMM_WORLD.Get_size());
					}
                }
				
				LSRK3MPIdata::t_synch_fs=0;
				LSRK3MPIdata::t_bp_fs=0;
				LSRK3MPIdata::t_fs=0;
				LSRK3MPIdata::t_up=0;
				LSRK3MPIdata::counter = 0;
				LSRK3MPIdata::t_diff = 0;
				LSRK3MPIdata::t_sten = 0;
            }
		}
		
		template< typename KernelType > 
		class IncomingBlocks : public tbb::task
		{
			TGrid& grid;
			Real current_time;
			SynchronizerMPI& synch;
			KernelType& rhs;
			
			bool master;
			vector<BlockInfo> avail;
			
		public:
			
			tbb::task* execute()
			{
				Timer timer;
				
				if (master)
				{
					tbb::empty_task& c = *new (allocate_continuation()) tbb::empty_task;
					c.set_ref_count(1);
					
					while(true)
					{
						timer.start();
						avail = synch.avail(LSRK3MPIdata::GSYNCH);			
						LSRK3MPIdata::t_synch_fs += timer.stop();
						LSRK3MPIdata::counter++;
						LSRK3MPIdata::nsynch++;
						
						const bool isdone = synch.done();
						
						if (!isdone)
						{
							IncomingBlocks& slave = * new (allocate_additional_child_of(c)) IncomingBlocks(synch, rhs, grid, current_time, false);
							slave.avail = avail;
							spawn(slave);
						}
						else 
						{
							recycle_as_child_of(c);
							master = false;
							
							return this;
						}
						
						if (isdone) break;
					}
					
					return NULL;
				}
				else
				{
					if (!LSRK3data::bAffinity)
                        BlockProcessingMPI::process< LabMPI >(avail, rhs, (TGrid&)grid, current_time);
				    else
                        BlockProcessingMPI::process< LabMPI >(avail, rhs, (TGrid&)grid, current_time, LSRK3data::affinitypart);
					
					return NULL;
				}
			}
			
			IncomingBlocks(SynchronizerMPI& synch, KernelType& rhs, TGrid& grid, Real current_time, bool master): 
			synch(synch), rhs(rhs), grid(grid), current_time(current_time), master(master)
			{
			}
		};
		
		template<typename KType > 
		void compute_asynch(SynchronizerMPI& synch, KType& rhs, TGrid& grid, Real current_time)
		{
			//printf("hello compute_asynch with %s\n", typeid(rhs).name());
			IncomingBlocks<KType> & dispatcher = *new (tbb::task::allocate_root()) IncomingBlocks<KType>(synch, rhs, grid, current_time, true);
			tbb::task::spawn_root_and_wait(dispatcher);
		}
		
		pair<double, double> step(TGrid& grid, vector<BlockInfo>& vInfo, Real a, Real b, Real dtinvh, const Real current_time)
		{
			Timer timer1, timer2;	
            
            if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("FLOWSTEP");                        
			LSRK3data::FlowStep<Kflow, Lab> rhs(a, dtinvh);           
            timer1.start();            
			
            SynchronizerMPI& synch = ((TGrid&)grid).sync(rhs);
            
			if (LSRK3data::dispatcher != "overlap")
				while (!synch.done())
				{			  
					timer2.start();
					vector<BlockInfo> avail = synch.avail(LSRK3MPIdata::GSYNCH);
					LSRK3MPIdata::t_synch_fs += timer2.stop();                           
					
					timer2.start();
					
					if (LSRK3data::dispatcher == "omp")
						_process< LabMPI >(avail, rhs, (TGrid&)grid, current_time);
					else
					{
						if (!LSRK3data::bAffinity)
							BlockProcessingMPI::process< LabMPI >(avail, rhs, (TGrid&)grid, current_time);
						else
							BlockProcessingMPI::process< LabMPI >(avail, rhs, (TGrid&)grid, current_time, LSRK3data::affinitypart);
					}
					
					LSRK3MPIdata::t_bp_fs += timer2.stop();
					
					LSRK3MPIdata::counter++;
					LSRK3MPIdata::nsynch++;
				}
            else 
			{
				timer2.start();
				compute_asynch(synch, rhs, grid, current_time);
				LSRK3MPIdata::t_bp_fs += timer2.stop();
			}
			
            double t1 = timer1.stop();
			if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();
            
            timer1.start();
            if(LSRK3data::sten_sigma>0)
            {
                if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("SURFACE TENSION");                        
                LSRK3data::SurfaceTension<Ksten, Lab> sten(dtinvh, LSRK3data::sten_sigma);                     
                
                SynchronizerMPI& synch = ((TGrid&)grid).sync(sten);
                
				if (LSRK3data::dispatcher != "overlap")
					while (!synch.done())
					{                
						vector<BlockInfo> avail = synch.avail(LSRK3MPIdata::GSYNCH);
						
						if (LSRK3data::dispatcher == "omp")
							_process< LabMPI >(avail, sten, (TGrid&)grid, current_time);
						else
							BlockProcessingMPI::process< LabMPI >(avail, sten, (TGrid&)grid, current_time);
					}
				else
					compute_asynch(synch, sten, grid, current_time);
                
                if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();
            }
            LSRK3MPIdata::t_sten += timer1.stop();
            
            timer1.start();
            if(LSRK3data::nu1>0)
            {
                if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("Diffusion");                        
                LSRK3data::Diffusion<Kdiff, Lab> diffusion(dtinvh, LSRK3data::nu1, LSRK3data::nu2);
				
                SynchronizerMPI& synch = ((TGrid&)grid).sync(diffusion);
                
				if (LSRK3data::dispatcher != "overlap")
					while (!synch.done())
					{                
						vector<BlockInfo> avail = synch.avail(LSRK3MPIdata::GSYNCH);
						
						if (LSRK3data::dispatcher == "omp")
							_process< LabMPI >(avail, diffusion, (TGrid&)grid, current_time);
						else
							BlockProcessingMPI::process< LabMPI >(avail, diffusion, (TGrid&)grid, current_time);
					}
				else
					compute_asynch(synch, diffusion, grid, current_time);
                
                if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();
            }
            LSRK3MPIdata::t_diff += timer1.stop();
            
			if (LSRK3data::profiler != NULL) LSRK3data::profiler->push_start("UPDATE");
			LSRK3data::Update<Kupdate> update(b, &vInfo.front());
			timer1.start();
			if (LSRK3data::dispatcher == "omp")
				update.omp(vInfo.size());
			else
			{
				if (!LSRK3data::bAffinity)
					parallel_for(blocked_range<int>(0, vInfo.size()), update, auto_partitioner());
				else
					parallel_for(blocked_range<int>(0, vInfo.size()), update, LSRK3data::affinitypart);
			} 
			
			double t2 = timer1.stop();
			if (LSRK3data::profiler != NULL) LSRK3data::profiler->pop_stop();
            
			LSRK3MPIdata::t_fs+=t1;
			LSRK3MPIdata::t_up+=t2;
            
			return pair<double, double>(t1,t2);
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
        
        if (verbosity>=1 && LSRK3data::step_id==0)
            cout << "Grid spacing and smoothing length are: " << h << ", " << smoothlength << endl; 
        
		LSRK3MPIdata::GSYNCH = parser("-gsync").asInt(LSRK3data::TLP);
        
		const Real maxSOS = _computeSOS();
		Real dt = min(max_dt, CFL*h/maxSOS);
		
		if (verbosity>=1)
			cout << "sos max is " << maxSOS << ", " << "advection dt is "<< dt << "\n";
		
        if (LSRK3data::nu1>0)
            dt = min(dt, (Real)(h*h/(12*max(LSRK3data::nu1,LSRK3data::nu2))) );
        
        if (LSRK3data::sten_sigma>0) 
        {
            const Real sumrho = parser("-sumrho").asDouble(HUGE_VAL);
            dt = min(dt, (Real)sqrt(sumrho*h*h*h/(4*M_PI*LSRK3data::sten_sigma)));
        }
        
        if (verbosity>=1)
			cout << "dt is "<< dt << "\n";
        
        if (dt<1e-10 || maxSOS<1e-10)
        {
            cout << "Something went wrong. Zero time step detected." << endl;
            abort();
        }
		
		if (parser("-kernels").asString("cpp")=="cpp")
			LSRKstepMPI<Convection_CPP, Update_CPP, SurfaceTension_CPP, Diffusion_CPP>(grid, dt/h, current_time);
		else if (parser("-kernels").asString("cpp")=="sse")
			LSRKstepMPI<Convection_SSE, Update_SSE, SurfaceTension_SSE, Diffusion_SSE>(grid, dt/h, current_time);
#ifdef _AVX_
    	else if (parser("-kernels").asString("cpp")=="avx")
			LSRKstepMPI<Convection_AVX, Update_AVX, SurfaceTension_AVX, Diffusion_AVX>(grid, dt/h, current_time);
#endif
		else
		{
			cout << "combination not supported yet" << endl;
			abort();
		}
		
		LSRK3data::step_id++; current_time+=dt;
		return dt;
	}
};

