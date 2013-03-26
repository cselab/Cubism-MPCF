/*
 *  main.cpp
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <mpi.h>
#ifdef _SSE_
#include <xmmintrin.h>
#endif

#include <ArgumentParser.h>
#include <Timer.h>
//#include <BlockProcessing.h>

#include "Tests.h"
#include "Test_SteadyStateMPI.h"
#include "Test_ShockBubbleMPI.h"
#include "Test_SICMPI.h"
#include "Test_CloudMPI.h"

using namespace std;

Simulation * sim = NULL;

int main (int argc, const char ** argv) 
{
	MPI::Init();
	
	const bool isroot = MPI::COMM_WORLD.Get_rank() == 0;

	if (isroot)
		cout << "=================  MPCF cluster =================" << endl;
	
	ArgumentParser parser(argc, argv);	
	const bool bFlush2Zero = parser("-f2z").asBool(true);
	
#ifdef _SSE_
	if (bFlush2Zero)
#pragma omp parallel
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	}
	
#endif

	parser.set_strict_mode();
	
	if (!isroot)
		parser.mute();
		
	if (isroot)
	  cout << "Dispatcher: " << parser("-dispatcher").asString() << endl;

	//Environment::setup(max(1, parser("-nthreads").asInt()));

	if( parser("-sim").asString() == "steady" )
		sim = new Test_SteadyStateMPI(isroot, argc, argv);
	else if( parser("-sim").asString() == "sb" )
		sim = new Test_ShockBubbleMPI(isroot, argc, argv);
	else if( parser("-sim").asString() == "sic" )
		sim = new Test_SICMPI(isroot, argc, argv);
	else if( parser("-sim").asString() == "cloud" )
		sim = new Test_CloudMPI(isroot, argc, argv);
	else
		if (isroot)
		{
			printf("-sim value not recognized. Aborting.\n");
			abort();
		}
		else abort();

	sim->setup();
	
	double wallclock;
	
	{
		Timer timer;
		
		timer.start();		
		
		sim->run();
		
		wallclock = timer.stop();
	}
	
	sim->dispose();
	
	delete sim;
	
	sim = NULL;
	
	if (isroot)
		printf("we spent: %2.2f \n", wallclock);
	
	MPI::COMM_WORLD.Barrier();
	MPI::Finalize();
	
	return 0;
}
