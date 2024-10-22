/*
 *  main.cpp
 *  MPCFcluster
 *
 *  Created by Diego Rossinelli on 11/25/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <map>
#include <iomanip>
#include <iostream>
#include <mpi.h>

#include "Test_IO_Compressed.h"

using namespace std;

Simulation * sim = NULL;

int main(const int argc, const char ** argv)
{
	MPI::Init();
	
	const bool isroot = MPI::COMM_WORLD.Get_rank() == 0;
	
	ArgumentParser parser(argc, argv);	

//	parser.set_strict_mode();
	
	if (!isroot)
		parser.mute();
    
	if( parser("-sim").asString() == "io" )
		sim = new Test_IO_Compressed(isroot, argc, argv);
	else
		if (isroot)
		{
			printf("-sim value not recognized. Aborting.\n");
                        printf("usage: %s [-sim io][-bpdx <nbx>][-bpdy <nby>][-bpdz <nbz>][-wtype <type>][-threshold <t>]\n", argv[0]);
			abort();
		}
        else abort();

	sim->setup();
	
	sim->dispose();
	
	delete sim;
	
	sim = NULL;
	
	if (isroot) printf("Finishing...");
	
	MPI::COMM_WORLD.Barrier();
	MPI::Finalize();
	
	return 0;
}
