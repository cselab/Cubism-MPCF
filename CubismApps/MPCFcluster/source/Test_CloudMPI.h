/*
 *  Test_CloudMPI.h
 *  MPCFcluster
 *
 *  Created by Babak Hejazialhosseini on 2/25/13.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <limits>
#include <Test_Cloud.h>
#include "Test_SICMPI.h"

#define Tshape shape

class Test_CloudMPI: public Test_Cloud
{
   	Test_SteadyStateMPI * t_ssmpi;
    Test_ShockBubbleMPI * t_sbmpi;

protected:
	int XPESIZE, YPESIZE, ZPESIZE;
    
	G * grid;
	FlowStep_LSRK3MPI<G> * stepper;
    
public:
	bool isroot;
    
	Test_CloudMPI(const bool isroot, const int argc, const char ** argv):
    Test_Cloud(argc, argv), isroot(isroot)
	{
        t_ssmpi = new Test_SteadyStateMPI(isroot, argc, argv);
        t_sbmpi = new Test_ShockBubbleMPI(isroot, argc, argv);
	}
    
	void setup()
	{
		if (isroot && VERBOSITY)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("///////////               TEST Cloud MPI         ///////////\n");
			printf("////////////////////////////////////////////////////////////\n");
		}
        
		_setup_constants();
        t_ssmpi->setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);
        
		if (!isroot)
			VERBOSITY = 0;
        
		grid = new G(XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ);
        
		assert(grid != NULL);
        
		stepper = new FlowStep_LSRK3MPI<G>(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, parser, VERBOSITY, &profiler,  Simulation_Environment::PC1, Simulation_Environment::PC2);
        
		if(bRESTART)
		{
			t_ssmpi->restart(*grid);
			t = t_ssmpi->get_time();
			step_id = t_ssmpi->get_stepid();
		}
		else
        {
            bRestartedSeed = parser("-seed").asBool(0);
            
            const int n_shapes = CloudData::n_shapes;
            
            Real bcast_buffer[4*n_shapes];
            vector< shape * > v_shapes;
            
            if (isroot)
            {
                Seed my_seed(CloudData::seed_s, CloudData::seed_e, CloudData::n_shapes);
                my_seed.make_shapes(bRestartedSeed);
                
                v_shapes = my_seed.get_vshapes();
                
                for(int i=0; i< n_shapes; i++)
                {
                    Real c[3*n_shapes];
                    bcast_buffer[0+i*4] = v_shapes[i]->get_rad();
                    v_shapes[i]->get_center(c);
                    bcast_buffer[1+i*4] = c[0];
                    bcast_buffer[2+i*4] = c[1];
                    bcast_buffer[3+i*4] = c[2];
                }
            }
            
            MPI::Cartcomm cartcomm = grid->getCartComm();
            cartcomm.Bcast(bcast_buffer, 4*n_shapes, MPI_REAL, 0);
            
            if (!isroot)
            {
                for(int i=0; i< n_shapes; i++)
                {
                    Real c[3], rad;
                    rad  = bcast_buffer[0+i*4];
                    c[0] = bcast_buffer[1+i*4];
                    c[1] = bcast_buffer[2+i*4];
                    c[2] = bcast_buffer[3+i*4];
                    
                    shape * cur_shape = new shape(c, rad);
                    v_shapes.push_back(cur_shape);
                }
            }
            
           //_my_ic(*grid, v_shapes);
           _my_ic_quad(*grid, v_shapes);
            
            v_shapes.clear();
        }
	}
    
	void run()
	{
		if (isroot) printf("HELLO RUN\n");
		bool bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
        
		while(bLoop)
		{
			if (isroot) printf("Step id %d,Time %f\n", step_id, t);
            
			if (step_id%DUMPPERIOD==0)
			{
				std::stringstream streamer;
				streamer<<"data-"<<step_id;
                t_ssmpi->dump(*grid, step_id, streamer.str());
				t_ssmpi->vp(*grid, step_id, bVP);
			}
            
			if (step_id%SAVEPERIOD==0)
				t_ssmpi->save(*grid, step_id, t);
            
			const Real dt = (*stepper)(TEND-t);
            
			if(step_id%10 == 0 && isroot && step_id > 0)
				profiler.printSummary();
            
            profiler.push_start("DUMP STATISTICS");
            if (step_id%ANALYSISPERIOD==0)
                t_sbmpi->dumpStatistics(*grid, step_id, t, dt);
            profiler.pop_stop();
            
            profiler.push_start("DUMP ANALYSIS");
            if (step_id%ANALYSISPERIOD==0)
                t_sbmpi->dumpAnalysis(*grid, step_id, t, dt);
            profiler.pop_stop();
            
			t+=dt;
			step_id++;
            
			bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
            
            if (dt==0)
                break;
		}
        
		std::stringstream streamer;
		streamer<<"data-"<<step_id;;
		t_ssmpi->dump(*grid, step_id, streamer.str());
        
		delete stepper;
		if (isroot) printf("Finishing RUN\n");

		MPI_Finalize();
	}
};

