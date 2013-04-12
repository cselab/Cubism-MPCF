/*
 *  Test_CloudMPI.h
 *  MPCFcluster
 *
 *  Created by Babak Hejazialhosseini on 2/25/13.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 *  *
 *	PATER NOSTER, qui es in caelis, sanctificetur nomen tuum. 
 *	Adveniat regnum tuum. 
 *	Fiat voluntas tua, sicut in caelo et in terra. 
 *	Panem nostrum quotidianum da nobis hodie, 
 *	et dimitte nobis debita nostra sicut et nos dimittimus debitoribus nostris. 
 *	Et ne nos inducas in tentationem, 
 *	sed libera nos a malo. 
 *	Amen.
 *
 */
 
#pragma once

#include <limits>
#include <Test_Cloud.h>
#include "Test_SICMPI.h"

#define Tshape shape

//typedef BlockLabMPI< BlockLabCloudLaplace< FluidBlock, std::allocator> > LabLaplace;


typedef BlockLabMPI< BlockLab< FluidBlock, std::allocator> > LabLaplace;

template<typename TLab>
struct GaussSeidel // , not. just 2ndorder bspline convolution
{
    StencilInfo stencil;
    
    int stencil_start[3];
    int stencil_end[3];
    
    GaussSeidel(): stencil(-1, -1, -1, +2, +2, +2, true, 2, 4, 7)
    {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }
    
    GaussSeidel(const GaussSeidel& c): stencil(-1, -1, -1, +2, +2, +2, true, 2, 4, 7)
    {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }
    
    inline void operator()(TLab& lab, const BlockInfo& info, FluidBlock& o) const
    {     
		const Real w[3] = { 0.25f, 0.5f, 0.25f };
		
		  for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
                for(int ix=0; ix<FluidBlock::sizeX; ix++)
                {
					double s = 0;
					
					for(int dz=-1; dz < 2; ++dz)
					for(int dy=-1; dy < 2; ++dy)
					for(int dx=-1; dx < 2; ++dx)
						s += w[dx + 1]* w[dy + 1]* w[dz + 1] * lab(ix + dx, iy + dy, iz + dz).energy;
					
					o(ix,iy,iz).dummy = s;
				}
    }
};

template<typename Lab, typename Operator, typename TGrid>
void _process_laplace(vector<BlockInfo>& vInfo, Operator rhs, TGrid& grid, const Real t=0, bool tensorial=false)
{
#pragma omp parallel
    {
        vector<BlockInfo>myInfo = vInfo;
        BlockInfo * ary = &myInfo.front();
        
        Operator myrhs = rhs;
        Lab mylab;
        
        const SynchronizerMPI& synch = grid.get_SynchronizerMPI(myrhs);
        
        const int N = vInfo.size();
        mylab.prepare(grid, synch);
        
#pragma omp for schedule(runtime)
        for(int i=0; i<N; i++)
        {
            mylab.load(ary[i], t);
            myrhs(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
        }
    }
}

template <typename TGrid>
void _process_update(const vector<BlockInfo>& vInfo, TGrid& grid)
{
#pragma omp parallel
    {
        vector<BlockInfo>myInfo = vInfo;
        const int N = vInfo.size();
        
#pragma omp for schedule(runtime)
        for(int i=0; i<N; i++)
        {
			FluidBlock& b = *(FluidBlock*)myInfo[i].ptrBlock;
			
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
						b(ix,iy,iz).energy = b(ix,iy,iz).dummy;
        }
    }
}

class BS4
{
public:
    static inline Real eval(Real x)
    {
        const Real t = fabs(x);
        
        if (t>2) return 0;
        
        if (t>1) return pow(2-t,3)/6;
        
        return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
    }
};

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
    
    void _relax_pressure(G & grid)
    {
        GaussSeidel< LabLaplace > gs;
               
        SynchronizerMPI& synch = ((G&)grid).sync(gs);
        
        while (!synch.done())
        {
            //if (isroot) printf("One avail loop of gs ");
            vector<BlockInfo> avail = synch.avail(1);
            
            _process_laplace< LabLaplace >(avail, gs, (G&)grid);
            //if (isroot) printf("... done\n");            
        }
        
        _process_update(grid.getBlocksInfo(), grid);
    }
    
	void setup()
	{
		_setup_constants();
        t_ssmpi->setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);

        if (!isroot)
			VERBOSITY = 0;
    
        if (isroot)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("///////////                                      ///////////\n");
			printf("///////////               TEST Cloud MPI         ///////////\n");
			printf("///////////                                      ///////////\n");
			printf("PATER NOSTER, qui es in caelis, sanctificetur nomen tuum.\n"); 
			printf("Adveniat regnum tuum.\n");
			printf("Fiat voluntas tua, sicut in caelo et in terra.\n");
			printf("Panem nostrum quotidianum da nobis hodie,\n");
			printf("et dimitte nobis debita nostra sicut et\nnos dimittimus debitoribus nostris.\n");
			printf("Et ne nos inducas in tentationem,\n");
			printf("sed libera nos a malo.\n");
			printf("Amen.\n");
			printf("////////////////////////////////////////////////////////////\n");
		}

	const double extent = parser("-extent").asDouble(1.0);
	grid = new G(XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ, extent);

	printf("rank %d local bpd %d %d %d\n", isroot, grid->getResidentBlocksPerDimension(0), grid->getResidentBlocksPerDimension(1), grid->getResidentBlocksPerDimension(2));
        
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
			const MPI::Cartcomm& mycartcomm = grid->getCartComm();
			//load the Cloud namespace and broadcast it
            {
				if(isroot)
					_initialize_cloud();
					
				mycartcomm.Bcast(&CloudData::n_shapes, 1, MPI::INT, 0);
				mycartcomm.Bcast(&CloudData::n_small, 1, MPI::INT, 0);
				mycartcomm.Bcast(&CloudData::small_count, 1, MPI::INT, 0);
				mycartcomm.Bcast(&CloudData::min_rad, 1, MPI_REAL, 0);
				mycartcomm.Bcast(&CloudData::max_rad, 1, MPI_REAL, 0);
				mycartcomm.Bcast(CloudData::seed_s, 3, MPI::DOUBLE, 0);
				mycartcomm.Bcast(CloudData::seed_e, 3, MPI::DOUBLE, 0);
			}
                                 
            Seed myseed(CloudData::seed_s, CloudData::seed_e);
            		
            //load the shapes from file and broadcast them		
            {
				vector<shape> v(CloudData::n_shapes);
				
				if(isroot)
				{
					myseed.make_shapes(CloudData::n_shapes);
					v = myseed.get_shapes();
				}
				
				mycartcomm.Bcast(&v.front(), v.size() * sizeof(shape), MPI::CHAR, 0);
				
				if (!isroot)
					myseed.set_shapes(v);
			}
            
			assert(myseed.get_shapes_size()>0 && myseed.get_shapes_size() == CloudData::n_shapes);
			
			int peidx[3];
			grid->peindex(peidx);
			const double spacing = grid->getH()*_BLOCKSIZE_;
			const double mystart[3] = {peidx[0]*BPDX*spacing, peidx[1]*BPDY*spacing, peidx[2]*BPDZ*spacing};
			const double myextent[3] = {BPDX*spacing, BPDY*spacing, BPDZ*spacing};
						
			myseed = myseed.retain_shapes(mystart,myextent);
			MPI::COMM_WORLD.Barrier();
            if (isroot) 
				cout << "Setting ic now...\n";
				
            _my_ic_quad(*grid, myseed);
            
            if (isroot) 
				cout << "done!"<< endl;
            
            if (isroot) 
				cout << "relaxing pressure a little bit..."<< endl;
            
            for(int i = 0; i < 2; ++i)
				_relax_pressure(*grid);
            
            if (isroot) 
				cout << "done!"<< endl;
            
            if (isroot) 
				printf("Setting energy now...");
				
            _set_energy(*grid);
            
            if (isroot) 
				printf("done\n");
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
                profiler.push_start("IO HDF");
                t_ssmpi->dump(*grid, step_id, streamer.str());
                profiler.pop_stop();
                profiler.push_start("IO WAVELET");
				t_ssmpi->vp(*grid, step_id, bVP);
                profiler.pop_stop();
			}
            
            profiler.push_start("SAVE");
			if (step_id%SAVEPERIOD==0)
				t_ssmpi->save(*grid, step_id, t);
            profiler.pop_stop();
            
            profiler.push_start("STEP");
            //stepper->set_CFL(CFL/(1+15*BS4::eval((t-0.3)/0.025)));
            if (isroot) printf("CFL is %f, original CFL is %f\n", stepper->CFL, CFL);
			const Real dt = (*stepper)(TEND-t);
            profiler.pop_stop();
            
            profiler.push_start("DUMP STATISTICS");
            if (step_id%ANALYSISPERIOD==0)
                t_sbmpi->dumpStatistics(*grid, step_id, t, dt);
            profiler.pop_stop();
            
            profiler.push_start("DUMP ANALYSIS");
            if (step_id%ANALYSISPERIOD==0)
                t_sbmpi->dumpAnalysis(*grid, step_id, t, dt);
            profiler.pop_stop();
            
            if(step_id%10 == 0 && isroot && step_id > 0)
				profiler.printSummary();
            
			t+=dt;
			step_id++;
            
			bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
            
            if (dt==0)
                break;
		}
        
		std::stringstream streamer;
		streamer<<"data-"<<step_id;;
		t_ssmpi->dump(*grid, step_id, streamer.str());
        
		if (isroot) 
			printf("Finishing RUN\n");
	}

	void dispose()
	{
			t_ssmpi->dispose();
			delete t_ssmpi;
			t_ssmpi = NULL;

			t_sbmpi->dispose();
			delete t_sbmpi;
			t_sbmpi = NULL;

			if (stepper != NULL)
			{
					delete stepper;
					stepper = NULL;
			}

			if (grid != NULL)
			{
					delete grid;
					grid = NULL;
			}
	}
};

