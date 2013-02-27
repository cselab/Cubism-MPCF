/*
 *  Test_Cloud.cpp
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 2/24/13.
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

#include "Test_Cloud.h"
#include "Tests.h"

void Test_Cloud::_ic(FluidGrid& grid)
{
    cout << "Not supposed to call _ic in Cloud\n" ;
    cout << "Should call _my_ic instead\n" ;
    abort();
}

void Test_Cloud::_my_ic(FluidGrid& grid, const vector< shape * > v_shapes)
{
	if (VERBOSITY > 0)
		cout << "Cloud Initial condition..." ;
    
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
    const double G1 = Simulation_Environment::GAMMA1-1;
    const double G2 = Simulation_Environment::GAMMA2-1;
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;
    
    const double h = vInfo.front().h_gridpoint;

#pragma omp parallel
	{
#ifdef _USE_NUMA_
		const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
		const int mynode = omp_get_thread_num() / cores_per_node;
		numa_run_on_node(mynode);
#endif
        
#pragma omp for
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
                        
                        const double bubble = eval(v_shapes, p);
                        
                        const Real pre_shock[3] = {1000,0,100};
                        
                        b(ix, iy, iz).rho      = 1.0*bubble+pre_shock[0]*(1-bubble);
                        b(ix, iy, iz).u        = 0;
                        b(ix, iy, iz).v        = 0;
                        b(ix, iy, iz).w        = 0;
                        
                        const double pressure  = 0.0234*bubble+pre_shock[2]*(1-bubble);
                        
                        SETUP_MARKERS_IC                       
			  }
        }
	}
	
	if (VERBOSITY > 0)
		cout << "done." << endl;
}

void Test_Cloud::setup()
{
  printf("////////////////////////////////////////////////////////////\n");
  printf("////////////             TEST CLOUD          ///////////////\n");
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
      bRestartedSeed = parser("-seed").asBool(0);
      
      Seed my_seed(CloudData::seed_s, CloudData::seed_e, CloudData::n_shapes);
      my_seed.make_shapes(bRestartedSeed);
      
      _my_ic(*grid, my_seed.get_vshapes());
  }
}
