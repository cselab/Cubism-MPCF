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

#ifdef _USE_HDF_
#include <HDF5Dumper.h>
#endif

using namespace std;

#ifdef _USE_HDF_
#include <HDF5Dumper.h>
#endif

#include "Test_SteadyState.h"
#include "SerializerIO_WaveletCompression.h"
#include "WaveletCompressor.h"

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
    const string path = parser("-fpath").asString(".");
	
#ifdef _USE_HDF_
    cout << "Dump to " << path << filename << "..." ;
    DumpHDF5<FluidGrid, StreamerGamma_HDF5>(*grid, step_id, filename+"-g", path);
    DumpHDF5<FluidGrid, StreamerPressure_HDF5>(*grid, step_id, filename+"-p", path);
    cout << "done." << endl;
#else
#warning HDF WAS DISABLED AT COMPILE TIME
#endif
}

void Test_SteadyState::_ic(FluidGrid& grid)
{
	//cout << "Initial condition..." ;
	
	const double G1 = Simulation_Environment::GAMMA1-1;
	const double G2 = Simulation_Environment::GAMMA2-1;
	const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
	const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;
	
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
					b(ix, iy, iz).rho      = 1.132;//sqrt(pow(p[0]-0.5,2)+pow(p[1]-0.5,2)+pow(p[2]-0.5,2))+0.1;
					b(ix, iy, iz).u        = -2.*b(ix, iy, iz).rho;//sqrt(pow(p[0]-0.3,2)+pow(p[1]-0.7,2)+pow(p[2]-0.35,2))-0.1;
					b(ix, iy, iz).v        = 3*b(ix, iy, iz).rho;//1.32;
					b(ix, iy, iz).w        = -100*b(ix, iy, iz).rho;//-11.2;
					const double pressure = 10;
					const double bubble = 0;

					SETUP_MARKERS_IC
				}
	}	
	
	//cout << "done." << endl;
}

void Test_SteadyState::_vp(FluidGrid& grid)
{
	char bufname[1024];
	sprintf(bufname, "compressed-step%05d.binary", step_id);
	
	_vp_dump(grid, bufname);
}

void Test_SteadyState::_vp_dump(FluidGrid& grid, string filename)
{
    if (bVP)
    {
		static const bool quantization = false;
		static const int NC = StreamerGridPoint::channels;
		static Real gmin[NC], gmax[NC];
		static bool reduced = false;
				
		SerializerIO_WaveletCompression<FluidGrid, StreamerGridPoint> wavelet_serializer;
		
		wavelet_serializer.set_threshold(1e-4);
		
		if (quantization)
			wavelet_serializer.float16();
		//wavelet_serializer.singlethreaded();
		//wavelet_serializer.nozlib();
		wavelet_serializer.Write(grid, filename);
		
		//for consistency checking, uncomment the following line
		wavelet_serializer.Read(grid, filename);
		//exit(0);
    }
}

void Test_SteadyState::run()
{
	for(int i=0; i<NSTEPS; ++i)
	{
	  	(*stepper)(TEND-t);
	}      
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
    ANALYSISPERIOD = parser("-analysisperiod").asInt(std::numeric_limits<int>::max());
    
	if(VERBOSITY) 
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
	_setup_constants();
	    
    if (VERBOSITY)
	{
		printf("////////////////////////////////////////////////////////////\n");
		printf("////////////         TEST STEADY STATE       ///////////////\n");
		printf("////////////////////////////////////////////////////////////\n");
	}
    
	grid = new FluidGrid(BPDX, BPDY, BPDZ);
	
	assert(grid != NULL);
	
	stepper = new FlowStep_LSRK3(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, parser, VERBOSITY);
	
	if(bRESTART)
	{
		_restart();
		_dump("restartedcondition");
	}
	else
	{
		_ic(*grid);
		_dump("initialcondition");
	}
}
