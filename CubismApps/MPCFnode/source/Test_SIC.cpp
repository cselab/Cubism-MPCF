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

#include "Test_SIC.h"
#include "Tests.h"

static const Real w_G[3]={ 5./9., 8./9., 5./9. };

template<int PointNumber>
static Real _getGaussPoints(Real x, Real h)
{
    Real xG[3]= {
        -0.5*sqrt(3./5.)*h+x,
        x,
        0.5*sqrt(3./5.)*h+x
    };
    return xG[PointNumber];
}

void Test_SIC::_ic_gauss(FluidGrid& grid)
{
	if (VERBOSITY > 0)
		cout << "SIC Initial condition..." ;
    
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
            
            b.clear_tmp();
            
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                for(int iy=0; iy<FluidBlock::sizeY; iy++)
                    for(int ix=0; ix<FluidBlock::sizeX; ix++)
                    {
                        Real p[3], post_shock[3];
                        info.pos(p, ix, iy, iz);
                        
                        const Real xG1D[3] = {
                            _getGaussPoints<0>(p[0], h),
                            _getGaussPoints<1>(p[0], h),
                            _getGaussPoints<2>(p[0], h)
                        };
                        
                        const Real yG1D[3] = {
                            _getGaussPoints<0>(p[1], h),
                            _getGaussPoints<1>(p[1], h),
                            _getGaussPoints<2>(p[1], h)
                        };
                        
                        
                        const Real zG1D[3] = {
                            _getGaussPoints<0>(p[2], h),
                            _getGaussPoints<1>(p[2], h),
                            _getGaussPoints<2>(p[2], h)
                        };
                        
                        for (int iG=0; iG<3; iG++)
                            for(int jG=0; jG<3; jG++)
                                for(int kG=0; kG<3; kG++)
                                {
                                    const Real xG3D[3]=
                                    {
                                        xG1D[iG],
                                        yG1D[jG],
                                        zG1D[kG]
                                    };
                                    
                                    const double r = sqrt(pow(xG3D[0]-bubble_pos[0],2)+pow(xG3D[1]-bubble_pos[1],2)+pow(xG3D[2]-bubble_pos[2],2));
                                    const double bubble =   Simulation_Environment::heaviside_smooth(r-radius);
                                    const double shock_pressure = 353;
                                    
                                    const Real pre_shock[3] = {10,0,1};
                                    Simulation_Environment::getPostShockRatio(pre_shock, shock_pressure, post_shock);
                                    
                                    const double shock = 0;//Simulation_Environment::heaviside_smooth(xG3D[0]-Simulation_Environment::shock_pos);
                                    
                                    const Real rho = shock*post_shock[0] + (1-shock)*(0.01*bubble+pre_shock[0]*(1-bubble));
                                    const Real u   = post_shock[0]*post_shock[1]*shock;
                                    const Real v   = 0;
                                    const Real w   = 0;
                                    const double pressure  = post_shock[2]*shock + pre_shock[2]*(1-shock);
                                    
                                    const double mix_gamma = 1 + (G2*G1)/(G1*bubble+G2*(1-bubble));
                                    const double mix_pinf  = (mix_gamma-1)/mix_gamma * (F1/G1*(1-bubble) + F2/G2*bubble);
                                    const Real G  = 1./(mix_gamma-1);
                                    const Real P = mix_gamma*mix_pinf/(mix_gamma-1);
                                    const double ke = 0.5*(pow(u,2)+pow(v,2)+pow(w,2))/rho;
                                    const Real energy   = pressure*G + P + ke;
                                    
                                    b(ix, iy, iz).rho			+=	(1./pow((Real)2,3))*w_G[iG]*w_G[jG]*w_G[kG]*rho;
                                    b(ix, iy, iz).u   			+=	(1./pow((Real)2,3))*w_G[iG]*w_G[jG]*w_G[kG]*u;
                                    b(ix, iy, iz).v             +=	(1./pow((Real)2,3))*w_G[iG]*w_G[jG]*w_G[kG]*v;
                                    b(ix, iy, iz).w             +=	(1./pow((Real)2,3))*w_G[iG]*w_G[jG]*w_G[kG]*w;
                                    b(ix, iy, iz).energy        +=	(1./pow((Real)2,3))*w_G[iG]*w_G[jG]*w_G[kG]*energy;
                                    b(ix, iy, iz).G             +=	(1./pow((Real)2,3))*w_G[iG]*w_G[jG]*w_G[kG]*G;
                                    b(ix, iy, iz).P             +=	(1./pow((Real)2,3))*w_G[iG]*w_G[jG]*w_G[kG]*P;
                                }
                        
                        //**************************
                        //Let's do 1D in Colonius
                        //**************************
                        /*Real p[3];
                         info.pos(p, ix, iy, iz);
                         
                         //test 5.3 equation 32
                         //const Real pre_shock[3] = {10,0,10};
                         //const Real post_shock[3] = {0.125,0,0.1};
                         //const double bubble = Simulation_Environment::heaviside(0.5-p[0]);
                         
                         //test 5.3 equation 31
                         const Real pre_shock[3] = {1.241,0,2.753};
                         const Real post_shock[3] = {0.991,0,3.059e-4};
                         const double bubble = Simulation_Environment::heaviside(0.5-p[0]);
                         
                         const double shock = 1-bubble;
                         b(ix, iy, iz).rho      =  shock*pre_shock[0] + (1-shock)*post_shock[0];
                         b(ix, iy, iz).u        = pre_shock[1]*b(ix, iy, iz).rho*shock;
                         b(ix, iy, iz).v        = 0;
                         b(ix, iy, iz).w        = 0;
                         
                         const double pressure  = pre_shock[2]*shock+post_shock[2]*(1-shock);
                         
                         SETUP_MARKERS_IC*/
                    }
        }
	}
	
	if (VERBOSITY > 0)
		cout << "done." << endl;
}

void Test_SIC::_ic(FluidGrid& grid)
{
	if (VERBOSITY > 0)
		cout << "SIC Initial condition..." ;
    
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
                        
                        const double r = sqrt(pow(p[0]-bubble_pos[0],2)+pow(p[1]-bubble_pos[1],2)+pow(p[2]-bubble_pos[2],2));
                        
                        const double bubble =   Simulation_Environment::heaviside_smooth(r-radius);
                        
                        const double shock_pressure = 3530;
                        
                        const Real pre_shock[3] = {10,0,10};

                        Simulation_Environment::getPostShockRatio(pre_shock, shock_pressure, post_shock);
                        
                        const double shock = Simulation_Environment::heaviside_smooth(p[0]-Simulation_Environment::shock_pos);           
                                                
                        b(ix, iy, iz).rho      = shock*post_shock[0] + (1-shock)*(0.01*bubble+pre_shock[0]*(1-bubble));
                        b(ix, iy, iz).u        = post_shock[0]*post_shock[1]*shock;
                        b(ix, iy, iz).v        = 0;
                        b(ix, iy, iz).w        = 0;
                        
                        const double pressure  = post_shock[2]*shock + pre_shock[2]*(1-shock);

                        SETUP_MARKERS_IC
                        
                        //**************************
                        //Let's do 1D in Colonius
                        //**************************
                        /*Real p[3];
                         info.pos(p, ix, iy, iz);
                         
                         //test 5.3 equation 32
                         //const Real pre_shock[3] = {10,0,10};
                         //const Real post_shock[3] = {0.125,0,0.1};
                         //const double bubble = Simulation_Environment::heaviside(0.5-p[0]);     
                         
                         //test 5.3 equation 31
                         const Real pre_shock[3] = {1.241,0,2.753};
                         const Real post_shock[3] = {0.991,0,3.059e-4};                       
                         const double bubble = Simulation_Environment::heaviside(0.5-p[0]);                                                                        
                         
                         const double shock = 1-bubble;
                         b(ix, iy, iz).rho      =  shock*pre_shock[0] + (1-shock)*post_shock[0];
                         b(ix, iy, iz).u        = pre_shock[1]*b(ix, iy, iz).rho*shock;
                         b(ix, iy, iz).v        = 0;
                         b(ix, iy, iz).w        = 0;
                         
                         const double pressure  = pre_shock[2]*shock+post_shock[2]*(1-shock);
                         
                         SETUP_MARKERS_IC*/
                    }
        }		
	}	
	
	if (VERBOSITY > 0)
		cout << "done." << endl;
}

void Test_SIC::setup()
{
    printf("////////////////////////////////////////////////////////////\n");
    printf("////////////             TEST SIC            ///////////////\n");
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
        _ic(*grid);
}
