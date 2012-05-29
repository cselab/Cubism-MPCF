/*
 *  Test_SteadyState.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cmath>
#include <iomanip>
#include <sstream>
using namespace std;

#include "Test_SteadyState.h"

Test_SteadyState::Test_SteadyState(const int argc, const char ** argv):
parser(argc, argv), t(0), step_id(0), grid(NULL), stepper(NULL) { }

void Test_SteadyState::_restart()
{	
	//read status
	{
		ifstream status("restart.status");
		assert(status.good());
		
		status >> t;
        status.ignore(1, ' ');
		status >> step_id;
		
		assert(t>=0);
		assert(step_id >= 0);
	}
	
	if (VERBOSITY) printf("DESERIALIZATION: time is %f and step id is %d\n", t, step_id);
	
	//read grid
	if (bASCIIFILES)
		SerializerIO<FluidGrid, StreamerGridPointASCII>().Read(*grid, "restart");
	else
		SerializerIO<FluidGrid, StreamerGridPoint>().Read(*grid, "restart");
}

void Test_SteadyState::_save()
{	
	cout << "Saving...";
	
	//write status
	{
		ofstream status("restart.status");
		
		status << t << " " << step_id;
		
		printf( "time: %20.20e\n", t);
		printf( "stepid: %d\n", step_id);
	}
	
	
	//write grid
	if (bASCIIFILES)
		SerializerIO<FluidGrid, StreamerGridPointASCII>().Write(*grid, "restart");
	else 
		SerializerIO<FluidGrid, StreamerGridPoint>().Write(*grid, "restart");
	
	cout << "done." << endl;
	
	{
		ofstream history("history.txt", step_id == 0? ios::out : ios::app);
		
		history << setprecision(10) << t << " " << step_id;
	}
}

void Test_SteadyState::_dump(string filename)
{	
	cout << "Dump to " << filename << "..." ;
#ifdef _USE_VTK_ 
	SerializerIO_ImageVTK<FluidGrid, StreamerGridPoint> vtkdumper;
	vtkdumper.Write(*grid, filename);
#else
	#warning VTK WAS DISABLED AT COMPILE TIME
#endif
    
	cout << "done." << endl;
}

void Test_SteadyState::_ic(FluidGrid& grid)
{
	cout << "Initial condition..." ;
	vector<BlockInfo> vInfo = grid.getBlocksInfo();
	
#pragma omp parallel for	
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
					b(ix, iy, iz).rho      = sqrt(pow(p[0]-0.5,2)+pow(p[1]-0.5,2)+pow(p[2]-0.5,2))+0.1;
					b(ix, iy, iz).u        = sqrt(pow(p[0]-0.3,2)+pow(p[1]-0.7,2)+pow(p[2]-0.35,2))-0.1;
					b(ix, iy, iz).v        = 1.32;
					b(ix, iy, iz).w        = -11.2;
					b(ix, iy, iz).energy   = 2.5+0.5*3;
					b(ix, iy, iz).G = Simulation_Environment::GAMMA1;
                    
#ifdef _LIQUID_
                    b(ix, iy, iz).P = Simulation_Environment::PC1;
#endif                    
				}
	}	
	
	cout << "done." << endl;
}

void Test_SteadyState::_vp(FluidGrid& grid)
{
    if (bVP)
    {
        std::stringstream streamer;
        streamer<<"data";
        streamer.setf(ios::dec | ios::right);
        streamer.width(5);
        streamer.fill('0');
        streamer<<0;
        streamer<<"_";
        streamer.setf(ios::dec | ios::right);
        streamer.width(5);
        streamer.fill('0');
        streamer<<step_id;
        
        SerializerIO_VP<FluidGrid, StreamerGridPoint> vp(BPDX, BPDY, BPDZ);
        vp.Write(grid, streamer.str());
    }
}

void Test_SteadyState::run()
{
	for(int i=0; i<NSTEPS; ++i)
	{
		(*stepper)(TEND-t);
	}
	
	_dump("ciao.vti");
	_save();
	_vp(*grid);
}

void Test_SteadyState::paint() { }

void Test_SteadyState::_setup_constants()
{    
    bRESTART = parser("-restart").asBool();
	
	parser.set_strict_mode();
	
	BPDX = parser("-bpdx").asInt();
	TEND = parser("-tend").asDouble();
	DUMPPERIOD = parser("-dumpperiod").asInt();
	SAVEPERIOD = parser("-saveperiod").asInt();
	CFL = parser("-cfl").asDouble();
	MOLLFACTOR = parser("-mollfactor").asInt();
    
	parser.unset_strict_mode();
	
	bASCIIFILES = parser("-ascii").asBool(false);
	BPDY = parser("-bpdy").asInt(BPDX);
	BPDZ = parser("-bpdz").asInt(BPDX);
	Simulation_Environment::GAMMA1 = parser("-g1").asDouble(1.4);
	Simulation_Environment::GAMMA2 = parser("-g2").asDouble(1.4);
	bVP = parser("-vp").asBool(0);
	VERBOSITY = parser("-verb").asInt(0);
	NSTEPS = parser("-nsteps").asInt(0);
	bAWK = parser("-awk").asBool(false);
	
	parser.save_options();
	parser.mute();
	
	assert(TEND >= 0.0);
	assert(BPDX >= 1);
	assert(BPDY >= 1);
	assert(BPDZ >= 1);
	assert(DUMPPERIOD > 0);
	assert(CFL > 0 && CFL<1);
	assert(MOLLFACTOR > 0);
	
    Simulation_Environment::PC1 = parser("-pc1").asDouble(0);
	Simulation_Environment::PC2 = parser("-pc2").asDouble(0);
    REPORT_FREQ = parser("-report").asInt(10);
}

void Test_SteadyState::setup()
{	
	if (VERBOSITY)
	{
		printf("////////////////////////////////////////////////////////////\n");
		printf("////////////         TEST STEADY STATE       ///////////////\n");
		printf("////////////////////////////////////////////////////////////\n");
	}
	
	_setup_constants();
	
	grid = new FluidGrid(BPDX, BPDY, BPDZ);
	
	assert(grid != NULL);
	
	stepper = new FlowStep_LSRK3(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, parser);
	
	if(bRESTART)
	{
		_restart();
		_dump("restartedcondition.vti");
	}
	else
	{
		_ic(*grid);
		_dump("initialcondition.vti");
	}
}
