/*
 *  Test_ShockTube.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <limits>
#include <sstream>
#include "Test_ShockTube.h"

void Test_ShockTube::_ic(FluidGrid& grid)
{
	cout << "ShockTube Initial condition..." ;
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	
    const double G1 = Simulation_Environment::GAMMA1-1;
    const double G2 = Simulation_Environment::GAMMA2-1;
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;
    
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					Real p[3];
					info.pos(p, ix, iy, iz);
                    
                    const bool mask = p[0]<0.5 && p[1]<0.5 || p[0]>0.5 && p[1]>0.5;
                    
					b(ix, iy, iz).rho      = mask? 1:0.125;
					b(ix, iy, iz).u        = 0;
					b(ix, iy, iz).v        = 0;
					b(ix, iy, iz).w        = 0;
					b(ix, iy, iz).energy   = mask? 2.5:0.25;
					b(ix, iy, iz).G = Simulation_Environment::GAMMA1;
                    
#ifdef _LIQUID_
                    b(ix, iy, iz).P = Simulation_Environment::PC1;
#endif
				}
	}
	
	cout << "done." << endl;
}

void Test_ShockTube::run()
{
	int step_id=0;	
	bool bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
    
    while(bLoop)
	{
		cout << "time is " << t << endl;
		cout << "step_id is " << step_id << endl;
        
		if(step_id%DUMPPERIOD == 0)
		{
			std::stringstream streamer;
			streamer<<"data-"<<step_id<<".vti";
			_dump(streamer.str());
		}
		
		const Real dt = (*stepper)(TEND-t);
		t+=dt;
		step_id++;
        bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
	}
}

void Test_ShockTube::setup()
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("////////////         TEST SHOCK TUBE         ///////////////\n");
	printf("////////////////////////////////////////////////////////////\n");
	
	_setup_constants();
	
	grid = new FluidGrid(BPDX, BPDY, BPDZ);
	assert(grid != NULL);
	
	stepper = new FlowStep_LSRK3(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA1, parser);
	
	_ic(*grid);
	_dump("initialcondition.vti");
}
