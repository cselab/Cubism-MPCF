/*
 *  Test_TG.cpp
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

#include "Test_TG.h"
#include "Tests.h"

void Test_TG::_ic(FluidGrid& grid)
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
											
						const double gamma  = Simulation_Environment::GAMMA1;
						
						//b(ix, iy, iz).G = 1/G1;
                        
						b(ix, iy, iz).rho      = 1;
						b(ix, iy, iz).u        = b(ix, iy, iz).rho*sin(2*M_PI*p[0])*cos(2*M_PI*p[1])*cos(2*M_PI*p[2])*(2*M_PI);
						b(ix, iy, iz).v        = -b(ix, iy, iz).rho*cos(2*M_PI*p[0])*sin(2*M_PI*p[1])*cos(2*M_PI*p[2])*(2*M_PI);
						b(ix, iy, iz).w        = 0;
						
						const double pressure  = 4*M_PI*M_PI/pow(Simulation_Environment::mach,2)/gamma + 1./16*(cos(4*M_PI*p[2])+2*cos(4*M_PI*p[0])+cos(4*M_PI*p[1])-2); // since rho=1
						b(ix, iy, iz).energy   = pressure/(gamma-1.) + 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;
					}
		}
		
	}	
	cout << "done." << endl;
}

void Test_TG::_dumpStatistics(FluidGrid& grid, const int step_id, const Real t, const Real dt)
{
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
	Real rInt=0., uInt=0., vInt=0., wInt=0., eInt=0., vol=0., ke=0.;
    Real ens=0, maxvor = 0;
    Real x[3];
    const Real h = vInfo[0].h_gridpoint;
    const Real h2 = h*h;
    const Real h3 = h2*h;
    
    Lab lab;
    const int ss[3] = {-1,-1,-1};
    const int se[3] = {2,2,2};
    lab.prepare(grid, ss, se, false);
    
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
        lab.load(info);
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
                    
                    ke   += 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w);
                    
                    const Real omega[3] = {lab(ix, iy+1, iz).w/lab(ix, iy+1, iz).rho-lab(ix, iy-1, iz).w/lab(ix, iy-1, iz).rho+lab(ix, iy, iz-1).v/lab(ix, iy, iz-1).rho-lab(ix, iy, iz+1).v/lab(ix, iy, iz+1).rho,
                        lab(ix, iy, iz+1).u/lab(ix, iy, iz+1).rho-lab(ix, iy, iz-1).u/lab(ix, iy, iz-1).rho+lab(ix-1, iy, iz).w/lab(ix-1, iy, iz).rho-lab(ix+1, iy, iz).w/lab(ix+1, iy, iz).rho,
                        lab(ix+1, iy, iz).v/lab(ix+1, iy, iz).rho-lab(ix-1, iy, iz).v/lab(ix-1, iy, iz).rho+lab(ix, iy-1, iz).u/lab(ix, iy-1, iz).rho-lab(ix, iy+1, iz).u/lab(ix, iy+1, iz).rho};
                    
                    ens += omega[0]*omega[0]+omega[1]*omega[1]+omega[2]*omega[2];
                    
                    maxvor = max(maxvor, (Real)(sqrt(omega[0]*omega[0]+omega[1]*omega[1]+omega[2]*omega[2])));
                }
    }
    
    ens *= 0.25*h;
    maxvor *= 0.5/h;
    
    FILE * f = fopen("integrals.dat", "a");
    fprintf(f, "%d  %f  %f  %f  %f  %f %f  %f  %f %f %f\n", step_id, t, dt, rInt*h3, uInt*h3, 
            vInt*h3, wInt*h3, eInt*h3, ke*h3, ens, maxvor);
    fclose(f);
}

void Test_TG::run()
{	
    Real dt=0;
	bool bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
	while (bLoop)
	{
		cout << "time is " << t << endl;
		cout << "step_id is " << step_id << endl;
        
		if(step_id%DUMPPERIOD == 0)
		{
			profiler.push_start("DUMP");
            
			std::stringstream streamer;
			streamer<<"data-"<<step_id<<".vti";
			_dump(streamer.str());
			//_vp(*grid);			
			profiler.pop_stop();
            
            profiler.push_start("DUMP STATISTICS");
            _dumpStatistics(*grid, step_id, t, dt);
            profiler.pop_stop();
		}
        
		if (step_id%SAVEPERIOD == 0) _save();
		
		profiler.push_start("EVOLVE");
		dt = (*stepper)(TEND-t);
		profiler.pop_stop();
		
		if(step_id%10 == 0)
			profiler.printSummary();			
		        
		t+=dt;
		step_id++;
		bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
	}
    
	std::stringstream streamer;
	streamer<<"data-"<<step_id<<".vti";
	if (DUMPPERIOD < 1e5) _dump(streamer.str());
}

void Test_TG::_setup_constants()
{
	Test_SteadyState::_setup_constants();
	
	parser.set_strict_mode();
	Simulation_Environment::mach          = parser("-mach").asDouble();
	parser.unset_strict_mode();
}

void Test_TG::setup()
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
