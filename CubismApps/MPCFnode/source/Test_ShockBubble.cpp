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

#include "Test_ShockBubble.h"
#include "Tests.h"

void Test_ShockBubble::_ic(FluidGrid& grid)
{
	cout << "ShockBubble Initial condition..." ;
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
                        const double r = sqrt(pow(p[0]-Simulation_Environment::shock_pos-1.2*radius,2)+pow(p[1]-bubble_pos[1],2));
                        const double bubble = Simulation_Environment::heaviside_smooth(r-radius);                                                                        
                        
                        const Real pre_shock[3] = {1,0,1};
                        Simulation_Environment::getPostShockRatio(pre_shock, Simulation_Environment::mach, Simulation_Environment::GAMMA1, Simulation_Environment::PC1, post_shock);	      
                        const double shock = Simulation_Environment::heaviside_smooth(p[0]-Simulation_Environment::shock_pos);                                           
                        
                        b(ix, iy, iz).rho      = shock*post_shock[0] + (1-shock)*(0.138*bubble+pre_shock[0]*(1-bubble));
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

Real _findzero(const Real f0, const Real f1, const Real f2, const Real h)
{
    const Real a = (f0-2*f1+f2)/(2*h*h);
    const Real b = (f2-f1-a*h*h)/h;
    const Real c = f1;
    const Real delta = b*b- 4*a*c;
    return (f1<=f2? (-b+sqrt(delta))/(2*a) : (-b-sqrt(delta))/(2*a));
}

void Test_ShockBubble::_dumpStatistics(FluidGrid& grid, const int step_id, const Real t, const Real dt)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
	Real rInt=0., uInt=0., vInt=0., wInt=0., eInt=0., vol=0., ke=0.;
    Real x[3];
    const Real h = vInfo[0].h_gridpoint;
    const Real h3 = h*h*h;
    
    Lab lab;
    const int ss[3] = {-1,-1,-1};
    const int se[3] = {2,2,2};
    lab.prepare(grid, ss, se, false);
    
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
                    info.pos(x, ix, iy, iz);
                    rInt += b(ix, iy, iz).rho;
                    uInt += b(ix, iy, iz).u;
                    vInt += b(ix, iy, iz).v;
                    wInt += b(ix, iy, iz).w;
                    eInt += b(ix, iy, iz).energy;
                    vol  += b(ix,iy,iz).G;
                    ke   += 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w);
                }
    }
    
    FILE * f = fopen("integrals.dat", "a");
    fprintf(f, "%d  %f  %f  %f  %f  %f  %f  %f %f   %f\n", step_id, t, dt, rInt*h3, uInt*h3, 
            vInt*h3, wInt*h3, eInt*h3, vol*h3, ke);
    fclose(f);
}

void Test_ShockBubble::run()
{	
	bool bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
	while (bLoop)
	{
		cout << "time is " << t << endl;
		cout << "step_id is " << step_id << endl;
        
		if(step_id%DUMPPERIOD == 0 && DUMPPERIOD < 1e5)
		{
			profiler.push_start("DUMP");
            
			std::stringstream streamer;
			streamer<<"data-"<<step_id<<".vti";
			_dump(streamer.str());
			//_vp(*grid);			
			profiler.pop_stop();
		}
        
		if (step_id%SAVEPERIOD == 0 && SAVEPERIOD < 1e5) _save();
		
		profiler.push_start("EVOLVE");
		const Real dt = (*stepper)(TEND-t);
		profiler.pop_stop();
		
		if(step_id%10 == 0)
			profiler.printSummary();			
		
        profiler.push_start("DUMP STATISTICS");
		if (DUMPPERIOD < 1e5) _dumpStatistics(*grid, step_id, t, dt);
		profiler.pop_stop();
        
		t+=dt;
		step_id++;
		bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
	}
    
	std::stringstream streamer;
	streamer<<"data-"<<step_id<<".vti";
	if (DUMPPERIOD < 1e5) _dump(streamer.str());
}

void Test_ShockBubble::_setup_constants()
{
	Test_SteadyState::_setup_constants();
	
	parser.set_strict_mode();
	Simulation_Environment::mach          = parser("-mach").asDouble();
	Simulation_Environment::shock_pos     = parser("-shockpos").asDouble();
	bubble_pos[0] = parser("-bubx").asDouble();
	bubble_pos[1] = parser("-buby").asDouble();
	bubble_pos[2] = parser("-bubz").asDouble();
	radius        = parser("-rad").asDouble();
	parser.unset_strict_mode();
}

void Test_ShockBubble::setup()
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
}
