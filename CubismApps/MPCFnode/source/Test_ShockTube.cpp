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
					b(ix, iy, iz).rho      = p[0]<0.5? 1:0.125;
					b(ix, iy, iz).u        = 0;
					b(ix, iy, iz).v        = 0;
					b(ix, iy, iz).w        = 0;
					b(ix, iy, iz).energy   = p[0]<0.5? 2.5:0.25;
					b(ix, iy, iz).levelset = 1;
				}
	}
	
	cout << "done." << endl;
}

void Test_ShockTube::run()
{
	int counter=0;
	
	while (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1)
	{
		cout << "time is " << t << endl;
		
		if(counter%SAVEPERIOD == 0)
		{
			std::stringstream streamer;
			streamer<<"data-"<<counter<<".vti";
			_dump(streamer.str());
		}
		
		const Real dt = (*stepper)(TEND-t);
		t+=dt;
		counter++;
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
