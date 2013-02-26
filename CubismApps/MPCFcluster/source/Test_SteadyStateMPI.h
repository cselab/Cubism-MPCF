/*
 *  Test_SteadyStateMPI.h
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <GridMPI.h>
#include <HDF5Dumper_MPI.h>

#include "FlowStep_LSRK3MPI.h"
#include <Test_SteadyState.h>

typedef GridMPI< FluidGrid > G;

class Test_SteadyStateMPI: public Test_SteadyState
{
protected:
	int XPESIZE, YPESIZE, ZPESIZE;
	
    G * grid;
	FlowStep_LSRK3MPI<G> * mystepper;
    
public:
	
	const bool isroot;
	
	Test_SteadyStateMPI(const bool isroot, const int argc, const char ** argv):
		Test_SteadyState(argc, argv), isroot(isroot) { }
        
    void setup_mpi_constants(int& xpesize, int& ypesize, int& zpesize)
	{	
		xpesize = parser("-xpesize").asInt(2);
		ypesize = parser("-ypesize").asInt(2);
		zpesize = parser("-zpesize").asInt(2);
	}
    
    void dump(G& grid, const int step_id, const string filename)
    {	
        if (isroot) cout << "Dumping " << "..." ;

        const string path = parser("-fpath").asString(".");
		DumpHDF5_MPI<G, StreamerGamma_HDF5>(grid, step_id, filename, path);
        
        if (isroot) cout << "done." << endl;
    }

    void restart(G& grid)
    {
		const string path = parser("-fpath").asString(".");
        
        {
            const string restart_status = path+"/restart.status";
            ifstream status(restart_status.c_str());
            assert(status.good());
           
            status >> t;
            status.ignore(1, ' ');
            status >> step_id;
            
            assert(t>=0);
            assert(step_id >= 0);
       }
        
        if (isroot) printf("DESERIALIZATION: time is %f and step id is %d\n", t, step_id);
        
        ReadHDF5_MPI<G, StreamerDummy_HDF5>(grid, "data_restart", path.c_str());
        DumpHDF5_MPI<G, StreamerDummy_HDF5>(grid, 0, "data_restart_restarted", path.c_str());
    }
    
    void save(G& grid, const int step_id, const Real t)
    {
      if (isroot) cout << "Saving...";
        
        const string path = parser("-fpath").asString(".");
       
        if (isroot)
        {
            const string restart_status = path+"/restart.status";
            ofstream status(restart_status.c_str());
            
            status << t << " " << step_id;
           
            printf( "time: %20.20e\n", t);
            printf( "stepid: %d\n", step_id);
        }
        
        DumpHDF5_MPI<G, StreamerDummy_HDF5>(grid, step_id, "data_restart", path.c_str());
        
         if (isroot) cout << "done" <<endl;
     
      }
    
    void vp(G& grid, const int step_id, const bool bVP)
    {
		if (bVP)
		{
			if (isroot) cout << "dumping MPI VP ..." ;
		   
			const string path = parser("-fpath").asString(".");
	 
			std::stringstream streamer;
			streamer<<path;
			streamer<<"/";
			streamer<<"data";
			streamer.setf(ios::dec | ios::right);
			streamer.width(5);
			streamer.fill('0');
			streamer<<MPI::COMM_WORLD.Get_rank();
			streamer<<"_";
			streamer.setf(ios::dec | ios::right);
			streamer.width(5);
			streamer.fill('0');
			streamer<<step_id;
			
			SerializerIO_VP<G, StreamerDensity> vp(BPDX, BPDY, BPDZ);
			vp.Write(grid, streamer.str());

			if (isroot) cout << "done" << endl;
		}
    }
    
	void setup()
	{
		if (VERBOSITY)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("////////////      TEST STEADY STATE   MPI    ///////////////\n");
			printf("////////////////////////////////////////////////////////////\n");
		}
		
		_setup_constants();
		setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);

		if (!isroot)
			VERBOSITY = 0;
		
		grid = new G(XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ);
		
		assert(grid != NULL);
		
		mystepper = new FlowStep_LSRK3MPI< G >(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, parser, VERBOSITY);
		
		//NUMA touch - was screwing up things on BGQ
		//(*mystepper)(0);
		
		if(bRESTART)
		  restart(*grid);
		else
		{
			_ic(*grid);
			dump(*grid, step_id, "mpi_initialcondition");
		}
	}	
	
	void run()
	{
		printf("HELLO RUN\n");
		
		while(t < TEND)
		{
            if (isroot) printf("Time is %f\n", t);
			const Real dt = (*mystepper)(TEND-t);
			t += dt;
		}
		
        printf("Finishing RUN\n");
	}
};
