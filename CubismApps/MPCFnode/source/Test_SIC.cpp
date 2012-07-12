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

void Test_SIC::_ic(FluidGrid& grid)
{
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

			const double shockvalue = 1e4;

                        const double r1 = sqrt(pow(p[0]-bubble_pos[0],2)+pow(p[1]-bubble_pos[1],2));
                        
                        const double bubble = Simulation_Environment::heaviside_smooth(r1-radius);                                                                        
                        
                        const Real pre_shock[3] = {10,0,10};
                        Simulation_Environment::getPostShockRatio(pre_shock, Simulation_Environment::mach, Simulation_Environment::GAMMA1, Simulation_Environment::PC1, post_shock);	      
                        const double shock = Simulation_Environment::heaviside_smooth(p[0]-Simulation_Environment::shock_pos);//*Simulation_Environment::heaviside_smooth(0.05-p[0]);                                           
                        const double shock2 = Simulation_Environment::heaviside(p[0]-Simulation_Environment::shock_pos);
                        
                        b(ix, iy, iz).rho      =  shock*pre_shock[0] + (1-shock)*(0.01*bubble+pre_shock[0]*(1-bubble));
                        
                        b(ix, iy, iz).v        = 0;
                        b(ix, iy, iz).w        = 0;
                        
                        const double pulse_period = 0.5;
                        const double c_liquid = sqrt(Simulation_Environment::GAMMA1*(pre_shock[2]+Simulation_Environment::PC1)/pre_shock[0]);
                        const double pulse_omega = 7.34;//(1.5*M_PI-M_PI/3)/pulse_period;
                        const double pulse_decay = 8.85;//pulse_omega * 11;
                        //const double pulse_amp = 1e4;
                        const double ramp = 1;//1.03*(1-exp(-742.87*max((Real)0.,(Real)(Simulation_Environment::shock_pos-p[0]) ) ) );
                        const double p_front = shockvalue;//pre_shock[2]+2*pulse_amp*exp(-pulse_decay*(Simulation_Environment::shock_pos-p[0]))*cos(pulse_omega*(Simulation_Environment::shock_pos-p[0])+M_PI/3);
                        const double pressure  = p_front*ramp*shock+pre_shock[2]*(1-shock);
                        const double c2_liquid = sqrt(Simulation_Environment::GAMMA1*(pressure+Simulation_Environment::PC1)/pre_shock[0]);
                        
                        b(ix, iy, iz).u        = shockvalue/pre_shock[0]/c2_liquid* b(ix, iy, iz).rho*shock;//5.0*shock* b(ix, iy, iz).rho;// shockval*2/pre_shock[0]/c_liquid/3
                        
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
    {
        _ic(*grid);
    }
}

