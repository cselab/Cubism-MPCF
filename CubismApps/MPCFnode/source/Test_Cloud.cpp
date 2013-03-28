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

namespace CloudData
{
    int n_shapes = 0;
    int n_small = 0;
    int small_count = 0;
    Real min_rad = 0;
    Real max_rad = 0;
    Real seed_s[3], seed_e[3];
}

Test_Cloud::Test_Cloud(const int argc, const char ** argv): Test_ShockBubble(argc, argv)
{
    _setup_constants();
}

void Test_Cloud::_initialize_cloud()
{
    ifstream f_read("cloud_config.dat");
    if(f_read)
    {
        if (VERBOSITY) cout << "cloud config file is there" << endl;
        f_read >> CloudData::n_shapes >> CloudData::n_small;
        f_read >> CloudData::min_rad >> CloudData::max_rad;
        f_read >> CloudData::seed_s[0] >> CloudData::seed_s[1] >> CloudData::seed_s[2];
        f_read >> CloudData::seed_e[0] >> CloudData::seed_e[1] >> CloudData::seed_e[2];
        f_read.close();
    }
    else
    {
        if (VERBOSITY) cout << "cloud config file not there...aborting" << endl;
        abort();
    }
    
    if (VERBOSITY)
        printf("cloud data: N %d Nsmall %d Rmin %f Rmax %f s=%f,%f,%f e=%f,%f,%f\n", CloudData::n_shapes, CloudData::n_small, CloudData::min_rad, CloudData::max_rad,
               CloudData::seed_s[0], CloudData::seed_s[1], CloudData::seed_s[2],
               CloudData::seed_e[0], CloudData::seed_e[1], CloudData::seed_e[2]);
}

void Test_Cloud::_ic(FluidGrid& grid)
{
    if (VERBOSITY) cout << "Not supposed to call _ic in Cloud\n" ;
    if (VERBOSITY) cout << "Should call _my_ic instead\n" ;
    abort();
}

FluidElement operator * (Real a, FluidElement gp)
{
    FluidElement out;
 
    out.rho = gp.rho * a;
    out.u = gp.u * a;
    out.v = gp.v * a;
    out.w = gp.w * a;
    out.energy = gp.energy * a;
    out.G = gp.G * a;
    out.P = gp.P * a;
    
    return out;
}

FluidElement operator + (FluidElement gpa, FluidElement gpb)
{
    FluidElement out;
    
    out.rho = gpa.rho+gpb.rho;
    out.u = gpa.u+gpb.u;
    out.v = gpa.v+gpb.v;
    out.w = gpa.w+gpb.w;
    out.energy = gpa.energy+gpb.energy;
    out.G = gpa.G+gpb.G;
    out.P = gpa.P+gpb.P;
    
    return out;
}

//ALERT: this method puts pressure into energy since we
//will solve a laplace p = 0 without going back and forth
//between energy and pressure
template<typename T>
T get_ic(const Real p[3], const vector< shape * > v_shapes)
{
    T out;
    
    const double G1 = Simulation_Environment::GAMMA1-1;
    const double G2 = Simulation_Environment::GAMMA2-1;
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;
    
    const double bubble = eval(v_shapes, p);
    
    const Real pre_shock[3] = {1000,0,100};
    
    out.rho      = 1.0*bubble+pre_shock[0]*(1-bubble);
    out.u        = 0;
    out.v        = 0;
    out.w        = 0;
    
    const double pressure  = 0.0234*bubble+pre_shock[2]*(1-bubble);
    
    const double mix_gamma = 1 + (G2*G1)/(G1*bubble+G2*(1-bubble));
    const double mix_pinf  = (mix_gamma-1)/mix_gamma * (F1/G1*(1-bubble) + F2/G2*bubble);
    out.G  = 1./(mix_gamma-1);
    out.P = mix_gamma*mix_pinf/(mix_gamma-1);
    out.energy   = pressure;
    
    return out;    
}

template<typename T>
T integral(const Real p[3], float h, const vector< shape * > v_shapes) // h should be cubism h/2
{
	T samples[3][3][3];
	T zintegrals[3][3];
	T yzintegrals[3];
	
	const Real x0[3] = {p[0] -h, p[1] - h, p[2] - h};
	
	for(int iz=0; iz<3; iz++)
		for(int iy=0; iy<3; iy++)
			for(int ix=0; ix<3; ix++)
            {
                const Real mypos[3] = {x0[0]+ix*h, x0[1]+iy*h, x0[2]+iz*h};
				samples[iz][iy][ix] = get_ic<T>(mypos, v_shapes);
            }
	
	for(int iy=0; iy<3; iy++)
		for(int ix=0; ix<3; ix++)
			zintegrals[iy][ix] = (1/6.) * samples[0][iy][ix]+(2./3) * samples[1][iy][ix]+(1/6.)* samples[2][iy][ix];
	
	for(int ix=0; ix<3; ix++)
		yzintegrals[ix] = (1./6)*zintegrals[0][ix] + (2./3)*zintegrals[1][ix]+(1./6)*zintegrals[2][ix];
	
	return (1./6) * yzintegrals[0]+(2./3) * yzintegrals[1]+(1./6)* yzintegrals[2];
}

void Test_Cloud::_set_energy(FluidGrid& grid)
{    
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
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
                        
                        b(ix,iy,iz).energy = b(ix,iy,iz).energy*b(ix,iy,iz).G+b(ix,iy,iz).P;
                    }
        }
	}
}

void Test_Cloud::_my_ic_quad(FluidGrid& grid, const vector< shape * > v_shapes)
{
	if (VERBOSITY)
		cout << "Cloud Initial condition..." ;
    
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
    
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
                        
                        b(ix,iy,iz) = integral<FluidElement>(p,0.5*h,v_shapes);
                    }
        }
	}
	
	if (VERBOSITY)
		cout << "done." << endl;
}

void Test_Cloud::_my_ic(FluidGrid& grid, const vector< shape * > v_shapes)
{
	if (VERBOSITY)
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
	
	if (VERBOSITY)
		cout << "done." << endl;
}

void Test_Cloud::setup()
{
  printf("////////////////////////////////////////////////////////////\n");
  printf("////////////             TEST CLOUD          ///////////////\n");
  printf("////////////////////////////////////////////////////////////\n");
    
  _setup_constants();
    printf("cloud verb is %d\n", VERBOSITY);
  
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
      
      //_my_ic(*grid, my_seed.get_vshapes());
      _my_ic_quad(*grid, my_seed.get_vshapes());
  }
}
