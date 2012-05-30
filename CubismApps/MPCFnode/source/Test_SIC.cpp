/*
 *  Test_ShockBubble.cpp
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include <limits>
#include <sstream>

#ifdef _USE_NUMA_
#include <numa.h>
#include <omp.h>
#endif

#include <Profiler.h>

#include "Test_SIC.h"
#include "Tests.h"

void Test_SIC::_ic(FluidGrid& grid)
{
	cout << "SIC Initial condition..." ;
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    const double G1 = Simulation_Environment::GAMMA1-1;
    const double G2 = Simulation_Environment::GAMMA2-1;
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;
    
#pragma omp parallel
	{	
#ifdef _USE_NUMA_
		const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
		const int mynode = omp_get_thread_num() / cores_per_node;
		numa_run_on_node(mynode);
#endif
		
		for(int i=0; i<(int)vInfo.size(); i++)
		{
            BlockInfo info = vInfo[i];
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;
            
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                for(int iy=0; iy<FluidBlock::sizeY; iy++)
                    for(int ix=0; ix<FluidBlock::sizeX; ix++)
                    {
                        Real p[3], post_shock[3];
                        info.pos(p, ix, iy, iz);
                        const double r1 = sqrt(pow(p[0]-Simulation_Environment::shock_pos-1.2*radius,2)+pow(p[1]-bubble_pos[1],2));
                        const double r2 = sqrt(pow(p[0]-Simulation_Environment::shock_pos-3.5*radius,2)+pow(p[1]-bubble_pos[1],2));
                        
                        const double bubble = Simulation_Environment::heaviside_smooth(min(r1-radius, r2-radius));                                                                        
                                                
                        const Real pre_shock[3] = {10,0,10};
                        Simulation_Environment::getPostShockRatio(pre_shock, Simulation_Environment::mach, Simulation_Environment::GAMMA1, Simulation_Environment::PC1, post_shock);	      
                        const double shock = Simulation_Environment::heaviside_smooth(p[0]-Simulation_Environment::shock_pos);                                           
                        
                        b(ix, iy, iz).rho      = shock*post_shock[0] + (1-shock)*(0.01*bubble+pre_shock[0]*(1-bubble));
                        b(ix, iy, iz).u        = (shock*post_shock[1] + (1-shock)*pre_shock[1])*b(ix, iy, iz).rho;
                        b(ix, iy, iz).v        = 0;
                        b(ix, iy, iz).w        = 0;
                        
                        const double pressure  = shock*post_shock[2] + (1-shock)*pre_shock[2];
                        
                        SETUP_MARKERS_IC
                    }
        }		
	}	
	cout << "done." << endl;
}
/*
void Test_SIC::setup()
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("////////////         TEST SHOCK BUBBLE       ///////////////\n");
	printf("////////////////////////////////////////////////////////////\n");
	
	_setup_constants();
	parser.mute();
	
	if (parser("-morton").asBool(0))
		grid = new GridMorton<FluidGrid>(BPDX, BPDY, BPDZ);
	else
		grid = new FluidGrid(BPDX, BPDY, BPDZ);
	
	assert(grid != NULL);
	
	stepper = new FlowStep_LSRK3(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, parser, VERBOSITY, &profiler, Simulation_Environment::PC1, Simulation_Environment::PC2, bAWK);
	
    
	if(bRESTART)
	{
		_restart();
		_dump("restartedcondition.vti");
	}
	else
	{
		_ic(*grid);
	}
}*/
