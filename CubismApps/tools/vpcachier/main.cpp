/*
 *  main.cpp
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <mpi.h>

#include "Reader_WaveletCompression.h"

int main()
{
	std::string path = "../../MPCFcluster/makefiles/datawavelet00000.StreamerGridPointIterative.channel0";
	
	MPI::Init();
	
	MPI::Intracomm& mycomm = MPI::COMM_WORLD;
	
	const int myrank = mycomm.Get_rank();
	Reader_WaveletCompressionMPI myreader(mycomm);
	
	myreader.load_file(path);
	
	MPI::Finalize();
	
	if (myrank == 0)
		std::cout << "Also tchuess zame.\n";
	
	return 0;
}
