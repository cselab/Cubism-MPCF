/*
 *  Test_ShockBubbleMPI.h
 *  MPCFcluster
 *
 *  Created by Babak Hejazialhosseini on 11/29/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <limits>

#include "Test_SteadyStateMPI.h"
#include <Test_ShockBubble.h>

class Test_ShockBubbleMPI: public Test_ShockBubble
{
	Test_SteadyStateMPI * t_ssmpi;
    
protected:
	int XPESIZE, YPESIZE, ZPESIZE;
    
	G * grid;
	FlowStep_LSRK3MPI<G> * stepper;
    
public:
	bool isroot;
    
	Test_ShockBubbleMPI(const bool isroot, const int argc, const char ** argv):
    Test_ShockBubble(argc, argv), isroot(isroot)
	{
		t_ssmpi = new Test_SteadyStateMPI(isroot, argc, argv);
	}
    
	void setup()
	{
		if (isroot && VERBOSITY)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("///////////   TEST SHOCK BUBBLE INTERACTION MPI  ///////////\n");
			printf("////////////////////////////////////////////////////////////\n");
		}
        
		_setup_constants();
		t_ssmpi->setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);
        
		if (!isroot)
			VERBOSITY = 0;
        
		grid = new G(XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ);
        
		assert(grid != NULL);
        
		stepper = new FlowStep_LSRK3MPI<G>(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, parser, VERBOSITY, &profiler);
        
		if(bRESTART)
		{
			t_ssmpi->restart(*grid);
			t = t_ssmpi->get_time();
			step_id = t_ssmpi->get_stepid();
		}
		else
			_ic(*grid);
	}	
    
    void dumpStatistics(G& grid, const int step_id, const Real t, const Real dt)
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        double rInt=0., uInt=0., vInt=0., wInt=0., eInt=0., vol=0., ke=0., r2Int=0., mach_max=-HUGE_VAL, p_max=-HUGE_VAL;
        double g_rInt=0., g_uInt=0., g_vInt=0., g_wInt=0., g_eInt=0., g_vol=0., g_ke=0., g_r2Int=0., g_mach_max=-HUGE_VAL, g_p_max=-HUGE_VAL;
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
        
        MPI::COMM_WORLD.Reduce(&rInt, &g_rInt, 1, MPI::FLOAT, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&uInt, &g_uInt, 1, MPI::FLOAT, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&vInt, &g_vInt, 1, MPI::FLOAT, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&wInt, &g_wInt, 1, MPI::FLOAT, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&eInt, &g_eInt, 1, MPI::FLOAT, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&vol, &g_vol, 1, MPI::FLOAT, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&ke, &g_ke, 1, MPI::FLOAT, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&r2Int, &g_r2Int, 1, MPI::FLOAT, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&mach_max, &g_mach_max, 1, MPI::FLOAT, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce(&p_max, &g_p_max, 1, MPI::FLOAT, MPI::MAX, 0);
        
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            FILE * f = fopen("integrals.dat", "a");
            fprintf(f, "%d  %e  %e  %e  %e  %e  %e  %e %e   %e %e   %e  %e\n", step_id, t, dt, g_rInt*h3, g_uInt*h3, 
                    g_vInt*h3, g_wInt*h3, g_eInt*h3, g_vol*h3, g_ke, g_r2Int*h3, g_mach_max, g_p_max);
            fclose(f);
        }
              
        Lab lab;
        const int ss[3]={0,0,0};
        const int se[3]={2,2,2};
        lab.prepare(grid, ss, se, false);
        
        vector<pair<Real,Real> > velocities;
        
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
                        
                        if( Simulation_Environment::heaviside(lab(ix, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))*                        Simulation_Environment::heaviside(lab(ix+1, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1))) == 0 &&
                           !(Simulation_Environment::heaviside(lab(ix, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))==0 && 
                             Simulation_Environment::heaviside(lab(ix+1, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))==0 ) ) 
                        {                        
                            velocities.push_back(pair<Real,Real>(x[0],b(ix, iy, iz).u/b(ix, iy, iz).rho));
                        }
                    }
        }
        
        std::sort(velocities.begin(), velocities.end(), sort_pred());
        
        FILE * f2;
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            std::stringstream streamer;
            streamer<<"centerline_velocities-"<<step_id<<".dat";
            f2 = fopen(streamer.str().c_str(), "a");
            streamer.str("");
        }
        
         vector<pair<Real,Real> > g_velocities;
        if (MPI::COMM_WORLD.Get_rank()==0) g_velocities.resize(velocities.size()*MPI::COMM_WORLD.Get_size());
        
         MPI::COMM_WORLD.Gather(&velocities.front(), velocities.size(), MPI_FLOAT, &g_velocities.front(), velocities.size(), MPI_FLOAT, 0);                
        
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            if (g_velocities.size()>0) fprintf(f2, "%d %e %e\n", step_id, g_velocities.front().second, g_velocities.back().second);
        
            fclose(f2);
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
				streamer<<"data-"<<step_id;;
				t_ssmpi->dump(*grid, step_id, streamer.str());
				t_ssmpi->vp(*grid, step_id, bVP);
			}
            
			if (step_id%SAVEPERIOD==0)
				t_ssmpi->save(*grid, step_id, t);
            
			const Real dt = (*stepper)(TEND-t);
            
			if(step_id%10 == 0 && isroot && step_id > 0)
				profiler.printSummary();
            
			dumpStatistics(*grid, step_id, t, dt);
            
			t+=dt;
			step_id++;
            
			bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
            
            if (dt==0)
                break;
		}
        
		std::stringstream streamer;
		streamer<<"data-"<<step_id;;
		t_ssmpi->dump(*grid, step_id, streamer.str());
        
		if (isroot) printf("Finishing RUN\n");
		MPI_Finalize();
		exit(0);
	}
};

