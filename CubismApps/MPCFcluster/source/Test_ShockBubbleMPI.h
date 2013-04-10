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
#include <SynchronizerMPI.h>

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
        _setup_constants();
        t_ssmpi->setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);
        
        if (!isroot)
			VERBOSITY = 0;
        
        if (VERBOSITY)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("///////////   TEST SHOCK BUBBLE INTERACTION MPI  ///////////\n");
			printf("////////////////////////////////////////////////////////////\n");
		}
        
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
        vector<BlockInfo> vInfo = grid.getResidentBlocksInfo();
        double rInt=0., uInt=0., vInt=0., wInt=0., eInt=0., vol=0., ke=0., r2Int=0., mach_max=-HUGE_VAL, p_max=-HUGE_VAL;
        double g_rInt=0., g_uInt=0., g_vInt=0., g_wInt=0., g_eInt=0., g_vol=0., g_ke=0., g_r2Int=0., g_mach_max=-HUGE_VAL, g_p_max=-HUGE_VAL;
        double wall_p_max=-HUGE_VAL, g_wall_p_max=-HUGE_VAL;
        
        vector<BlockInfo> g_vInfo = grid.getBlocksInfo();
        const double h = g_vInfo.front().h_gridpoint;
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
                        vol  += b(ix,iy,iz).G>0.9*(1/(Simulation_Environment::GAMMA1-1)+1/(Simulation_Environment::GAMMA2-1))? 1:0;
                        r2Int += b(ix, iy, iz).rho*(1-min(max((b(ix,iy,iz).G-1/(Simulation_Environment::GAMMA1-2))/(1/(Simulation_Environment::GAMMA1-1)-1/(Simulation_Environment::GAMMA2-1)),(Real)0),(Real)1));
                        ke   += 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w);
                        
                        const double pressure = (b(ix, iy, iz).energy - 0.5/b(ix, iy, iz).rho * (b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w) - b(ix, iy, iz).P)/b(ix,iy,iz).G;
                        const double c = sqrt((1/b(ix,iy,iz).G+1)*(pressure+b(ix,iy,iz).P/b(ix,iy,iz).G/(1/b(ix,iy,iz).G+1))/b(ix, iy, iz).rho);
                        
                        const double velmag = sqrt(b(ix, iy, iz).u*b(ix, iy, iz).u+b(ix, iy, iz).v*b(ix, iy, iz).v+b(ix, iy, iz).w*b(ix, iy, iz).w)/b(ix, iy, iz).rho;
                        mach_max = max(mach_max, velmag/c);
                        p_max = max(p_max, pressure);
                        
                        if (info.index[2]==0 && iz==0)
                            wall_p_max = max(wall_p_max, pressure);
                    }
        }
        
        MPI::COMM_WORLD.Reduce(&rInt, &g_rInt, 1, MPI::DOUBLE, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&uInt, &g_uInt, 1, MPI::DOUBLE, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&vInt, &g_vInt, 1, MPI::DOUBLE, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&wInt, &g_wInt, 1, MPI::DOUBLE, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&eInt, &g_eInt, 1, MPI::DOUBLE, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&vol, &g_vol, 1, MPI::DOUBLE, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&ke, &g_ke, 1, MPI::DOUBLE, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&r2Int, &g_r2Int, 1, MPI::DOUBLE, MPI::SUM, 0);
        MPI::COMM_WORLD.Reduce(&mach_max, &g_mach_max, 1, MPI::DOUBLE, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce(&p_max, &g_p_max, 1, MPI::DOUBLE, MPI::MAX, 0);
        MPI::COMM_WORLD.Reduce(&wall_p_max, &g_wall_p_max, 1, MPI::DOUBLE, MPI::MAX, 0);
        
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            FILE * f = fopen("integrals.dat", "a");
            fprintf(f, "%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", step_id, t, dt, g_rInt*h3, g_uInt*h3,
                    g_vInt*h3, g_wInt*h3, g_eInt*h3, g_vol*h3, g_ke*h3, g_r2Int*h3, g_mach_max, g_p_max, pow(0.75*g_vol*h3/M_PI,1./3.), g_wall_p_max);
            fclose(f);
        }
        
    }
    
    struct Dummy
    {
        StencilInfo stencil;
        
        Dummy(): stencil(-1,-1,-1,2,2,2,false,1, 5) { }
        
        Dummy(const Dummy& c): stencil(c.stencil) { }
    };
    
    void dumpAnalysis(G& grid, const int step_id, const Real t, const Real dt)
    {
        std::stringstream streamer;

        vector<pair<Real,Real> > velocities;        
        vector<pair<Real,Real> > pressures;        
        
        Real x[3];
        
        Dummy dummy;
        
        SynchronizerMPI& synch = grid.sync(dummy);
        
        while (!synch.done())
        {
            vector<BlockInfo> avail = synch.avail(LSRK3MPIdata::GSYNCH);
            
            LabMPI lab;
            const SynchronizerMPI& Synch = grid.get_SynchronizerMPI(dummy);
            lab.prepare(grid, Synch);
            
            for(int i=0; i<(int)avail.size(); i++)
            {        
                BlockInfo info = avail[i];            
                
                if (info.index[1]!=grid.getBlocksPerDimension(1)/2) continue;
                
                FluidBlock& b = *(FluidBlock*)info.ptrBlock;
                lab.load(info, t);
                
                for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                    for(int iy=0; iy<FluidBlock::sizeY; iy++)
                        for(int ix=0; ix<FluidBlock::sizeX; ix++)
                        {
                            if (iy!=0 || iz!=0) continue;
                            
                            info.pos(x,ix,iy,iz);
                            
                            if( Simulation_Environment::heaviside(lab(ix, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))*
                               Simulation_Environment::heaviside(lab(ix+1, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1))) == 0 &&
                               !(Simulation_Environment::heaviside(lab(ix, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))==0 && 
                                 Simulation_Environment::heaviside(lab(ix+1, iy, iz).G-0.5*(1/(Simulation_Environment::GAMMA2-1)+1/(Simulation_Environment::GAMMA1-1)))==0 ) ) 
                            {                        
                                velocities.push_back(pair<Real,Real>(x[0],b(ix, iy, iz).u/b(ix, iy, iz).rho));
                            }
                            
                            const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;
                            
                            const double pressure = (b(ix, iy, iz).energy - ke -  b(ix, iy, iz).P)/b(ix, iy, iz).G;
                            pressures.push_back(pair<Real,Real>(x[0], pressure));
                        }
            }
        }
        
        velocities.resize(2, pair<Real,Real>(0,0));
        
        FILE * f2;
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            streamer<<"centerline_velocities.dat";
            f2 = fopen(streamer.str().c_str(), "a");
            streamer.str("");
        }
        
        vector<pair<Real,Real> > g_velocities;
        if (MPI::COMM_WORLD.Get_rank()==0) g_velocities.resize(velocities.size()*MPI::COMM_WORLD.Get_size());
        
        MPI::COMM_WORLD.Gather(&velocities.front(), 2*velocities.size(), MPI_FLOAT, &g_velocities.front(), 2*velocities.size(), MPI_FLOAT, 0);                
        
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            std::sort(g_velocities.begin(), g_velocities.end(), sort_pred());
            
            if (g_velocities.size()>0) fprintf(f2, "%d %f %e %e\n", step_id, t, g_velocities[velocities.size()*MPI::COMM_WORLD.Get_size()-2].second, g_velocities[velocities.size()*MPI::COMM_WORLD.Get_size()-1].second);
            
            fclose(f2);
        }        
        
        FILE * f3;
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            streamer<<"centerline_pressure-"<<step_id<<".dat";
            f3 = fopen(streamer.str().c_str(), "w");
            streamer.str("");
        }
        
        int pressures_lsize;
        int pressures_size = pressures.size();
        MPI::COMM_WORLD.Allreduce(&pressures_size, &pressures_lsize, 1, MPI::INT, MPI::MAX);
        pressures.resize(pressures_lsize, pair<Real,Real>(-1,0));

        vector<pair<Real,Real> > g_pressures;
        if (MPI::COMM_WORLD.Get_rank()==0) g_pressures.resize(pressures.size()*MPI::COMM_WORLD.Get_size());
        
        MPI::COMM_WORLD.Gather(&pressures.front(), 2*pressures.size(), MPI_FLOAT, &g_pressures.front(), 2*pressures.size(), MPI_FLOAT, 0);                
        
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            std::sort(g_pressures.begin(), g_pressures.end(), sort_pred());
            
            for(int i=0; i< g_pressures.size(); ++i)
            {
                if(g_pressures[i].first==-1) continue;
                fprintf(f3, "%e %e\n", g_pressures[i].first, g_pressures[i].second);
            }
            
            fclose(f3);
        }
        
        pressures.clear();
        g_pressures.clear();
        
        FILE * f4;
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            streamer<<"wall_pressure-"<<step_id<<".dat";
            f4 = fopen(streamer.str().c_str(), "w");
            streamer.str("");
        }
        
        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        
        for(int i=0; i<(int)vInfo.size(); i++)
        {        
            BlockInfo info = vInfo[i];
            
            if (info.index[2]!=0) continue;
            
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;
            
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                for(int iy=0; iy<FluidBlock::sizeY; iy++)
                    for(int ix=0; ix<FluidBlock::sizeX; ix++)
                    {
                        if (iz!=0) continue;
                        
                        info.pos(x,ix,iy,iz);
                        
                        const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho;
                        const double pressure = (b(ix, iy, iz).energy - ke -  b(ix, iy, iz).P)/b(ix, iy, iz).G;
                        
                        pressures.push_back(pair<Real,Real>(x[1], pressure));
                    }
        }
        
        pressures_size = pressures.size();
        MPI::COMM_WORLD.Allreduce(&pressures_size, &pressures_lsize, 1, MPI::INT, MPI::MAX);
        pressures.resize(pressures_lsize, pair<Real,Real>(-1,0));
        
        if (MPI::COMM_WORLD.Get_rank()==0) g_pressures.resize(pressures.size()*MPI::COMM_WORLD.Get_size());
        
        MPI::COMM_WORLD.Gather(&pressures.front(), 2*pressures.size(), MPI_FLOAT, &g_pressures.front(), 2*pressures.size(), MPI_FLOAT, 0);
        
        if (MPI::COMM_WORLD.Get_rank()==0)
        {
            std::sort(g_pressures.begin(), g_pressures.end(), sort_pred());
            
            for(int i=0; i< g_pressures.size(); ++i)
            {
                if(g_pressures[i].first==-1) continue;
                fprintf(f4, "%e %e\n", g_pressures[i].first, g_pressures[i].second);
            }
            
            fclose(f4);
        }
    }
    
	void run()
	{
		if (isroot) printf("HELLO RUN\n");
		bool bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
        
		while(bLoop)
		{
			if (isroot) printf("Step id %d,Time %f\n", step_id, t);

#ifdef _USE_HPM_
                        HPM_Start("Dumping");
#endif            
			if (step_id%DUMPPERIOD==0)
			{
				std::stringstream streamer;
				streamer<<"data-"<<step_id;;
				t_ssmpi->dump(*grid, step_id, streamer.str());
				t_ssmpi->vp(*grid, step_id, bVP);
			}
            
			if (step_id%SAVEPERIOD==0)
				t_ssmpi->save(*grid, step_id, t);
#ifdef _USE_HPM_
                        HPM_Stop("Dumping");
#endif            
			const Real dt = (*stepper)(TEND-t);
            
			if(step_id%10 == 0 && isroot && step_id > 0)
				profiler.printSummary();
#ifdef _USE_HPM_
			HPM_Start("Diagnostics");
#endif
            if (step_id%ANALYSISPERIOD==0)
            {
	      ;//dumpStatistics(*grid, step_id, t, dt);
	      ;// dumpAnalysis(*grid, step_id, t, dt);
            }
#ifdef _USE_HPM_
	    HPM_Stop("Diagnostics");
#endif            
			t+=dt;
			step_id++;
            
			bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
            
            if (dt==0)
                break;
		}
        
		std::stringstream streamer;
		streamer<<"data-"<<step_id;;
		t_ssmpi->dump(*grid, step_id, streamer.str());

		if (stepper != NULL) {
			delete stepper;        
			stepper = NULL;
		}
		if (isroot) printf("Finishing RUN\n");
	}
};

