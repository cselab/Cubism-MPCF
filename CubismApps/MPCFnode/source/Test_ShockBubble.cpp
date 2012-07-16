/*
 *  Test_ShockBubble.cpp
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include <sstream>
#include <algorithm>

#ifdef _USE_NUMA_
#include <numa.h>
#include <omp.h>
#endif

#include <Profiler.h>

#include "Test_ShockBubble.h"
#include "Tests.h"

using namespace std;

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
    double rInt=0., uInt=0., vInt=0., wInt=0., eInt=0., vol=0., ke=0., r2Int=0., mach_max=-HUGE_VAL, p_max=-HUGE_VAL;
    const double h = vInfo[0].h_gridpoint;
    const double h3 = h*h*h;
    
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
                    rInt += b(ix, iy, iz).rho;
                    uInt += b(ix, iy, iz).u;
                    vInt += b(ix, iy, iz).v;
                    wInt += b(ix, iy, iz).w;
                    eInt += b(ix, iy, iz).energy;
                    vol  += b(ix,iy,iz).G>0.5*(1/(Simulation_Environment::GAMMA1-1)+1/(Simulation_Environment::GAMMA2-1))? 1:0;
                    r2Int += b(ix, iy, iz).rho*(1-min(max((b(ix,iy,iz).G-1/(Simulation_Environment::GAMMA1-2))/(1/(Simulation_Environment::GAMMA1-1)-1/(Simulation_Environment::GAMMA2-1)),(Real)0),(Real)1));
                    ke   += 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w);
                    
#ifndef _LIQUID_
                    const double pressure = (b(ix, iy, iz).energy - 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w))/b(ix,iy,iz).G;
                    const double c = sqrt((1/b(ix,iy,iz).G+1)*pressure/b(ix, iy, iz).rho);
#else
                    const double pressure = (b(ix, iy, iz).energy - 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w) - b(ix, iy, iz).P)/b(ix,iy,iz).G;
                    const double c = sqrt((1/b(ix,iy,iz).G+1)*(pressure+b(ix,iy,iz).P/b(ix,iy,iz).G/(1/b(ix,iy,iz).G+1))/b(ix, iy, iz).rho);
#endif
                    
                    const double velmag = sqrt(b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w)/b(ix, iy, iz).rho;
                    mach_max = max(mach_max, velmag/c);
                    p_max = max(p_max, pressure);
                }
    }
    
    FILE * f = fopen("integrals.dat", "a");
    fprintf(f, "%d  %e  %e  %e  %e  %e  %e  %e %e   %e %e   %e  %e\n", step_id, t, dt, rInt*h3, uInt*h3, 
            vInt*h3, wInt*h3, eInt*h3, vol*h3, ke, r2Int*h3, mach_max, p_max);
    fclose(f);
    
    
    FILE * f2 = fopen("centerline_velocities.dat", "a");
    FILE * f3 = fopen("p_center.dat", "a");
    FILE * f4 = fopen("interface.dat","a");
    
    Lab lab;
    const int ss[3]={0,0,0};
    const int se[3]={2,2,2};
    lab.prepare(grid, ss, se, false);
    
    vector<pair<Real,Real> > velocities;
    vector<Real> iso_gamma;
    Real p_wall=0;
    
    Real x[3];
    
    for(int i=0; i<(int)vInfo.size(); i++)
    {        
        BlockInfo info = vInfo[i];
        
        if (info.index[1]!=grid.getBlocksPerDimension(1)/2) continue;
        
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        lab.load(info);
        
        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    if (iy!=0 || iz!=0) continue;
                    
                    info.pos(x,ix,iy,iz);
                    
                    if((ix+1)==FluidBlock::sizeX && (info.index[0]+1)==grid.getBlocksPerDimension(0)){
                        
                        const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;
#ifndef _LIQUID_
                        p_wall = (b(ix, iy, iz).energy - ke)/b(ix, iy, iz).G;
#else                    
                        p_wall = (b(ix, iy, iz).energy - ke -  b(ix, iy, iz).P)/b(ix, iy, iz).G;
#endif
                    }
                    
                    if (abs(b(ix,iy,iz).G-2.35)<0.1) {
                        iso_gamma.push_back(x[0]);
                    }
                    
                    if( Simulation_Environment::heaviside(lab(ix, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))*                        Simulation_Environment::heaviside(lab(ix+1, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1))) == 0 &&
                       !(Simulation_Environment::heaviside(lab(ix, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))==0 && 
                         Simulation_Environment::heaviside(lab(ix+1, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))==0 ) ) 
                    {                   
#ifndef _LIQUID_
                        const double pressure = (b(ix, iy, iz).energy - 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w))/b(ix,iy,iz).G;                        
#else
                        const double pressure = (b(ix, iy, iz).energy - 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w) - b(ix, iy, iz).P)/b(ix,iy,iz).G;
#endif
                        velocities.push_back(pair<Real,Real>(sqrt(6.59*b(ix,iy,iz).rho*pressure),b(ix, iy, iz).u/b(ix, iy, iz).rho));
                    }
                }
    }
    
    std::sort(velocities.begin(), velocities.end(), sort_pred());
    std::sort(iso_gamma.begin(),iso_gamma.end());
    
    if (velocities.size()>0) fprintf(f2, "%d %f %e %e %f\n", step_id, t, velocities.front().second, velocities.back().second, velocities.front().first);
    if (iso_gamma.size()>0) fprintf(f4, "%d %e %f %f\n", step_id, t, iso_gamma.front(), iso_gamma.back());
    fprintf(f3, "%d %e %e\n", step_id, t, p_wall);
    
    fclose(f4);
    fclose(f3);
    fclose(f2);
}

void Test_ShockBubble::_analysis(FluidGrid& grid, const int step_id)
{
    std::stringstream streamer;
    streamer<<"centerline_pressure-"<<step_id<<".dat";
    FILE * f = fopen(streamer.str().c_str(), "w");
    streamer.str("");
    
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    Real x[3];
    
    for(int i=0; i<(int)vInfo.size(); i++)
    {        
        BlockInfo info = vInfo[i];
        
        if (info.index[1]!=grid.getBlocksPerDimension(1)/2) continue;
        
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        
        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    if (iy!=0 || iz!=0) continue;
                    
                    info.pos(x,ix,iy,iz);
                    
                    const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;
                    
#ifndef _LIQUID_
                    const double pressure = (b(ix, iy, iz).energy - ke)/b(ix, iy, iz).G;
#else                    
                    const double pressure = (b(ix, iy, iz).energy - ke -  b(ix, iy, iz).P)/b(ix, iy, iz).G;
#endif                    
                    fprintf(f, "%e %e\n", x[0], pressure);
                }
    }    
    
    fclose(f); 
    
    streamer<<"wall_pressure-"<<step_id<<".dat";
    f = fopen(streamer.str().c_str(), "w");
    streamer.str("");
    
    for(int i=0; i<(int)vInfo.size(); i++)
    {        
        BlockInfo info = vInfo[i];
        
        if (info.index[0]!=grid.getBlocksPerDimension(0)-1 || info.index[2]!=grid.getBlocksPerDimension(2)/2) continue;
        
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        
        for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
                    if (ix!=FluidBlock::sizeX-1 || iz!=0) continue;
                    
                    info.pos(x,ix,iy,iz);
                    
                    const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;
                    const double pressure = (b(ix, iy, iz).energy - ke -  b(ix, iy, iz).P)/b(ix, iy, iz).G;
                    
                    fprintf(f, "%e %e\n", x[1], pressure);
                }
    }
    
    fclose(f);
}

void Test_ShockBubble::run()
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
            
            profiler.push_start("DUMP ANALYSIS");
            _analysis(*grid, step_id);
            profiler.pop_stop();
		}
        
		if (step_id%SAVEPERIOD == 0) _save();
		
		profiler.push_start("EVOLVE");

        stepper->set_current_time(t);
		dt = (*stepper)(TEND-t);

		profiler.pop_stop();
		
		if(step_id%10 == 0)
			profiler.printSummary();			
		
        profiler.push_start("DUMP STATISTICS");
        _dumpStatistics(*grid, step_id, t, dt);
        profiler.pop_stop();
                
		t+=dt;
		step_id++;
		bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
        if (dt==0) break;
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
