/*
 *  Test_CVT.cpp
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include <limits>
#include <sstream>
#include <Profiler.h>

#include "Test_CVT.h"
#include "Tests.h"

void Test_CVT::_ic(FluidGrid& grid)
{
    cout << "CVT Initial condition..." ;
    vector<BlockInfo> vInfo = grid.getBlocksInfo();
    BlockInfo info = vInfo[0];
	
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
						Real p[3];
						info.pos(p, ix, iy, iz);
						
						const double X = p[0];
						const double Y = p[1];
						const double Z = p[2];
						
						const double factK = 0.5*exp(2.0)*log(2.0);
						const double ampl = 0.2/(2.0*M_PI);
						const double k = 1.0*(2.0*M_PI);
						const double incl = 60.0 * M_PI/180.0;
						const double vsrad = 0.4/(2.0*M_PI);
						const double vtdist = 1.732/(2.0*M_PI);
						const double rthresh = 0.666/(2.0*M_PI);
						const double center[3]={0.5,0.5,0.5};
						
						const double omega0 = re*nu1/(0.3832515/pow(2*M_PI,2));
						const double wavefact = k;
						
						const double vxt = ampl*cos(incl)*(1-cos(wavefact*Z));
						const double vyt = ampl*sin(incl)*(1-cos(wavefact*Z));
						const double vpxt = ampl*wavefact*cos(incl)*sin(wavefact*Z);
						const double vpyt = ampl*wavefact*sin(incl)*sin(wavefact*Z);
						const double rad = sqrt( pow(X - center[0]-vxt+vtdist*0.5,2.0)+pow(Y-center[1]-vyt,2.0) );
						const double rad_rel = rad/rthresh;
						const double wz0 = (rad_rel<=1.0) ? omega0*(1.0 - exp(-factK/rad_rel* exp(1./(rad_rel-1.)) ) ): 0.0;
						
						const double vxt2 = ampl*cos(M_PI-incl)*(1-cos(wavefact*Z));
						const double vyt2 = ampl*sin(M_PI-incl)*(1-cos(wavefact*Z));
						const double vpxt2 = ampl*wavefact*cos(M_PI-incl)*sin(wavefact*Z);
						const double vpyt2 = ampl*wavefact*sin(M_PI-incl)*sin(wavefact*Z);
						const double rad2 = sqrt( pow(X- center[0]-vxt2-vtdist*0.5,2.0)+pow(Y- center[1]-vyt2,2.0));
						const double rad_rel2 = rad2/rthresh;
						const double wz02 = (rad_rel2<=1.0) ? omega0*(1.0 - exp(-factK/rad_rel2* exp(1./(rad_rel2-1.)) ) ): 0.0;
						
						b(ix, iy, iz).levelset = 1;
						
						b(ix, iy, iz).rho = 1.0;
						b(ix, iy, iz).u =  (wz0*vpxt - wz02*vpxt2)*b(ix, iy, iz).rho;
						b(ix, iy, iz).v =  (wz0*vpyt - wz02*vpyt2)*b(ix, iy, iz).rho;
						b(ix, iy, iz).w =  (wz0      - wz02)*b(ix, iy, iz).rho;
						
						b(ix, iy, iz).energy = 0;
					}
		}
		
		cout << "done." << endl;
	}
}
void Test_CVT::_dumpStatistics(FluidGrid& grid, const int step_id, const Real time, const Real dt)
{
	const string path = parser("-fpath").asString(".");
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	
	double rInt=0, uInt=0, vInt=0, wInt=0, sInt=0, omegaMin=HUGE_VAL, omegaMax=-HUGE_VAL, densityMin=HUGE_VAL, densityMax=-HUGE_VAL, omegaMaxPiS=-HUGE_VAL, omegaMaxPiD=-HUGE_VAL;
	double gammaPiS=0, gammaPiD=0, kenergy=0, kenstrophy=0;
	double x[3];        
	
	EvaluateVorticity< FluidGrid, EvaluateVorticity_CPP, Lab > eval_vort(grid);        
	eval_vort.execute();
	
	const Real h = vInfo[0].h_gridpoint;
    const Real h2 = h*h;
	const Real h3 = h*h2;
    
    for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
                    info.pos(x, ix, iy, iz);
                    
					rInt += b(ix,iy,iz).rho*h3;
					uInt += b(ix,iy,iz).u*h3;
					vInt += b(ix,iy,iz).v*h3;
					wInt += b(ix,iy,iz).w*h3;
					sInt += b(ix,iy,iz).energy*h3;
					
					const Real  vorticity_magnitude =  sqrt(pow(b.tmp[iz][iy][ix][0],2)+pow(b.tmp[iz][iy][ix][1],2)+pow(b.tmp[iz][iy][ix][2],2));
					omegaMin = min((Real)omegaMin,(Real)vorticity_magnitude);
					omegaMax = max((Real)omegaMax, (Real)vorticity_magnitude);
					
					const bool bLastZ = info.index[2]==BPDZ-1 && iz==FluidBlock::sizeZ-1;
					const bool bLastX = info.index[0]==BPDX-1 && ix==FluidBlock::sizeX-1;
					
					const bool bExcludeX = false;//info_global.index[2]==0 && iz<2;
					
					if (bLastX && !bExcludeX)
					{
						omegaMaxPiS = max((Real)omegaMaxPiS, (Real)vorticity_magnitude);
						gammaPiS += b.tmp[iz][iy][ix][0];
					}
					
					if (bLastZ)
					{
						omegaMaxPiD = max((Real)omegaMaxPiD, (Real)vorticity_magnitude);
						gammaPiD += b.tmp[iz][iy][ix][2];
					}
					
					densityMin = min((Real)densityMin, (Real)b(ix,iy,iz).rho);
					densityMax = max((Real)densityMax, (Real)b(ix,iy,iz).rho);
					
					kenergy    += (pow(b(ix,iy,iz).u,2)+pow(b(ix,iy,iz).v,2)+pow(b(ix,iy,iz).w,2))/pow(b(ix,iy,iz).rho,2) * h3;
					kenstrophy += pow(vorticity_magnitude,2) * h3;
					
					b.tmp[iz][iy][ix][0]=0;
					b.tmp[iz][iy][ix][1]=0;
					b.tmp[iz][iy][ix][2]=0;
                }
    }
	
	if (step_id==0)
	{
		gamma_max = gammaPiD;
		
		const string restart_status = path+"/restart_omega.status";
		ofstream status(restart_status.c_str());
		assert(status.good());
		status << gamma_max;
	}
	else if (bRESTART)
	{
		const string restart_status = path+"/restart_omega.status";
		ifstream status(restart_status.c_str());
		assert(status.good());
		status >> gamma_max;
	}
	else
		cout << "**** Maximum inital vorticity is not set. Either it is not t=0 or no data was available to the restarted run" << endl;
	
    
    FILE * f = fopen("integrals.dat", "a");
	FILE * f2 = fopen("ranges.dat", "a");
	
	fprintf(f,"%d %e %e %e %e %e %e %e %e %e %e\n", step_id, time, rInt, uInt, vInt, wInt, sInt, gammaPiS/gamma_max, gammaPiD/gamma_max, kenergy, kenstrophy);
	fprintf(f2,"%d %e %e %e %e %e %e %e\n", step_id, time, densityMin, densityMax, omegaMin, omegaMax, omegaMaxPiS, omegaMaxPiD);
	
	fclose(f);
	fclose(f2);
}

void Test_CVT::_restart()
{	
	const string path = parser("-fpath").asString(".");
	
	//read status
	{
		ifstream status("restart.status");
		assert(status.good());
		cout << "restart.status" << endl;
		
		status >> t;
        status.ignore(1, ' ');
		status >> step_id;
		
		assert(t>=0);
		assert(step_id >= 0);
	}
	
	
#ifdef _USE_HDF_
	std::stringstream streamer;
	streamer << "data_restart";
	std::stringstream streamer_restarted;
	streamer_restarted << "data_restart_restared";

	if (VERBOSITY) printf("DESERIALIZATION from HDF5: time is %f and step id is %d\n", t, step_id);
	
	ReadHDF5<FluidGrid, StreamerDummy_HDF5>(*grid, step_id, streamer.str(), path.c_str());
	//DumpHDF5<FluidGrid, StreamerDummy_HDF5>(*grid, step_id, streamer_restarted.str(), path.c_str());
#else
	cout << "Reading IC from " << path+"/restart" << endl;
	
	if (bASCIIFILES)
		SerializerIO<FluidGrid, StreamerGridPointASCII>().Read(*grid, path+"/restart");
	else 
		SerializerIO<FluidGrid, StreamerGridPoint>().Read(*grid, path+"/restart");
#endif
}

void Test_CVT::run()
{	
	const string path = parser("-fpath").asString(".");
	bool bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
	
	while (bLoop)
	{
		cout << "time is " << t << endl;
		cout << "step_id is " << step_id << endl;
		/*
		if(step_id%DUMPPERIOD == 0)
		{
			profiler.push_start("DUMP");
			std::stringstream streamer_vort;
			streamer_vort<<"vorticity-"<<step_id;
			DumpHDF5_vorticity< FluidGrid, Lab > (*grid, step_id, streamer_vort.str(), path);
            
			std::stringstream streamer;
			streamer<<path<<"/data-"<<step_id<<".vti";
			_dump(streamer.str());
			profiler.pop_stop();
		}
        
		if (step_id%SAVEPERIOD == 0) _save();
		*/
		profiler.push_start("EVOLVE");
		const Real dt = (*stepper)(TEND-t);
		profiler.pop_stop();
		
		if(step_id%10 == 0)
			profiler.printSummary();			
		
        profiler.push_start("DUMP STATISTICS");
		_dumpStatistics(*grid, step_id, t, dt);
		profiler.pop_stop();
        
		t+=dt;
		step_id++;
		
		bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
	}
    /*
	std::stringstream streamer;
	streamer<<path<<"/data-"<<step_id<<".vti";
	_dump(streamer.str());
	 */
}

void Test_CVT::_setup_constants()
{
	Test_SteadyState::_setup_constants();
	
	parser.set_strict_mode();
	Simulation_Environment::mach          = parser("-mach").asDouble();
	Simulation_Environment::shock_pos     = parser("-shockpos").asDouble();
	bubble_pos[0] = parser("-bubx").asDouble();
	bubble_pos[1] = parser("-buby").asDouble();
	bubble_pos[2] = parser("-bubz").asDouble();
	radius        = parser("-rad").asDouble();
    nu1    = parser("-nu1").asDouble(0);
    re     = parser("-re").asDouble();
	parser.unset_strict_mode();
}

void Test_CVT::setup()
{
	printf("////////////////////////////////////////////////////////////\n");
	printf("//////     TEST Compressible Vortex Reconnection   /////////\n");
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
/*		
#ifdef _USE_HDF_
		const string path = parser("-fpath").asString(".");
		std::stringstream streamer_vort;
		streamer_vort<<"vorticity-restarted-"<<step_id;
		DumpHDF5_vorticity< FluidGrid, Lab > (*grid, step_id, streamer_vort.str(), path);
#endif
 */
	}
	else
	{
		cout << "Only MPI IC currently supported! Aborting!" << endl;
		abort();
		_ic(*grid);
	}
}
