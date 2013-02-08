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

#define E_SCALE 1e5
#define R_SCALE 10

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
                                    const double shock_pressure = 3530;
                                    
                                    const Real pre_shock[3] = {10,0,10};
                                    Simulation_Environment::getPostShockRatio(pre_shock, shock_pressure, post_shock);
                                    
                                    const double shock = Simulation_Environment::heaviside_smooth(xG3D[0]-Simulation_Environment::shock_pos);
                                    
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
                        
                        const double r = sqrt(pow(p[0]-bubble_pos[0],2));//sqrt(pow(p[0]-bubble_pos[0],2)+pow(p[1]-bubble_pos[1],2)+pow(p[2]-bubble_pos[2],2));
                        
                        const double bubble =   Simulation_Environment::heaviside_smooth(r-radius);
                        
                        const double shock_pressure = 3530;
                        
                        const Real pre_shock[3] = {10,0,10};//{10,0,10};
                        
                        Simulation_Environment::getPostShockRatio(pre_shock, shock_pressure, post_shock);
                        
                        const double shock = 0;//Simulation_Environment::heaviside_smooth(p[0]-Simulation_Environment::shock_pos);//Simulation_Environment::heaviside_smooth(p[0]-Simulation_Environment::shock_pos);
                        
                        b(ix, iy, iz).rho      = shock*post_shock[0] + (1-shock)*(0.3*bubble+pre_shock[0]*(1-bubble));
                        b(ix, iy, iz).u        = post_shock[0]*post_shock[1]*shock;
                        b(ix, iy, iz).v        = 0;
                        b(ix, iy, iz).w        = 0;
                        
                        const double pressure  = post_shock[2]*shock + (1-shock)*(0.294*bubble+pre_shock[2]*(1-bubble));//pre_shock[2]*(1-shock);
                        
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

void Test_SIC::_ic_LUT(FluidGrid& grid)
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
                        
                        const double r = sqrt(pow(p[0]-bubble_pos[0],2));
                        
                        const double bubble =   Simulation_Environment::heaviside_smooth(r-radius);
                        
                        Temperature T_liquid = fromcelsius(27.0);
                        Pressure p_liquid = 1.0 * bar;//Pa
                        SteamCalculator S_liquid;
                        S_liquid.set_pT(p_liquid,T_liquid);
                        Density rho_liquid = S_liquid.dens();
                        SpecificEnergy e_liquid = S_liquid.specienergy();
                        
                        Pressure p_vapor = Boundaries::getSatPres_T(T_liquid);
                        SteamCalculator S_vapor;
                        S_vapor.set_pT(p_vapor,T_liquid);
                        Density rho_vapor = S_vapor.dens();
                        SpecificEnergy e_vapor = S_vapor.specienergy();
                        
                        const double rho_l = *reinterpret_cast<const double*>(&rho_liquid);
                        const double rho_v = *reinterpret_cast<const double*>(&rho_vapor);
                        const double e_l = *reinterpret_cast<const double*>(&e_liquid);
                        const double e_v = *reinterpret_cast<const double*>(&e_vapor);
                        
                        
                        b(ix, iy, iz).rho      = (bubble*rho_v + (1-bubble)*rho_l)/R_SCALE;//freesteam_rho(S_mixture)/R_SCALE;
                        b(ix, iy, iz).u        = 0;
                        b(ix, iy, iz).v        = 0;
                        b(ix, iy, iz).w        = 0;
                        
                        b(ix, iy, iz).P        = 0;
                        b(ix, iy, iz).energy   = (bubble*e_v*rho_v + (1-bubble)*e_l*rho_l)/E_SCALE;
                        
                        SpecificVolume v_mix = 1.0/(b(ix, iy, iz).rho * R_SCALE) * m3_kg;
                        SpecificEnergy e_mix = b(ix, iy, iz).energy * E_SCALE * J_kg;
                        
                        Solver2<SpecificEnergy,Density> SUR;
                        SteamCalculator S_mixture;
                        S_mixture = SUR.solve(e_mix,v_mix);
                        b(ix,iy,iz).G = v_mix/;//vapor mass fraction
                        
//                        SETUP_MARKERS_IC
                        
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
        _ic_LUT(*grid);
}

/*SteamCalculator S1, S2, S3;
 
 // turn on display of units
 cerr.flags(ios_base::showbase);
 
 // initialise T1, p1, p2
 Temperature T1 = fromcelsius(200.0);
 Pressure p1 = 5.0 * bar;
 Pressure p2 = 10.0 * bar;
 
 // Part 1
 cerr << endl << "Part (1) - density at 10 bar, 200?C" << endl;
 cerr << "p1 = " << p1/bar << " bar" << endl;
 cerr << "T1 = " << T1;
 cerr << " (" << tocelsius(T1) << "?C)" << endl;
 S1.set_pT(p1, T1);
 Density rho1 = S1.dens();
 cerr << "rho1 = " << rho1 << endl;
 
 // Part 2
 cerr << endl << "Part (2) - isentropic compression to 10 bar" << endl;
 SpecificEntropy s1 = S1.specentropy();
 cerr << "s1 = " << s1 << endl;
 cerr << "p2 = " << p2/bar << " bar" << endl;
 Solver<Pressure,SpecificEntropy,Temperature> PS(p2, s1);
 S2 = PS.solve(0.0001 * kJ_kgK, 0.0001 * Kelvin);
 Temperature T2 = S2.temp();
 cerr << "T2 = " << T2;
 cerr << " (" << tocelsius(T2) << "?C)" << endl;
 
 // part (3) - Finding p,T for v as above and u increased by 200 kJ_kg
 cerr << endl << "Part (3) - Finding p,T for v as above and u increased by 200 kJ_kg" << endl;
 SpecificVolume v2 = S2.specvol();
 SpecificEnergy u2 = S2.specienergy();
 cerr << "v2 = " << v2 << endl;
 cerr << "u2 = " << u2 << endl;
 
 SpecificEnergy u3 = u2 + 200.0 * kJ_kg;
 cerr << "u3 = " << u3 << endl;
 
 Solver2<SpecificEnergy,SpecificVolume> SUV;
 S3 = SUV.solve(u3,v2);
 Temperature T3 = S3.temp();
 Pressure p3 = S3.pres();
 cerr << "p3 = " << p3/bar << " bar" << endl;
 cerr << "T3 = " << T3;
 cerr << " (" << tocelsius(T3) << "?C)" << endl;*/

