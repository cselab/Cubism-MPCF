/*
 *  main.cpp
 *
 *
 *  Created by Panagiotis Chatzidoukas on 3/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <sstream>
#include <mpi.h>

#include <float.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//#include "ParIO.h"
#include <omp.h>
//#define _TRANSPOSE_DATA_

#include "ArgumentParser.h"
#include "Reader_WaveletCompression.h"
#include "Reader_WaveletCompression0.h"

int main(int argc, char **argv)
{
	const double init_t0 = omp_get_wtime();

	/* Initialize MPI */
	MPI_Init(&argc, &argv);

//	blosc_init();

#if defined(_USE_SZ_)||defined(_USE_SZ3_)
	printf("sz.config...\n");
	SZ_Init((char *)"sz.config");
	omp_set_num_threads(1);
#endif

	/* MPI variables */
	MPI_Comm comm  = MPI_COMM_WORLD;
	MPI_Info info  = MPI_INFO_NULL;
	MPI::Intracomm& mycomm = MPI::COMM_WORLD;

	const int mpi_rank = mycomm.Get_rank();
	const int mpi_size = mycomm.Get_size();

	const bool isroot = !mpi_rank;

	ArgumentParser argparser(argc, (const char **)argv);

	if (isroot)
		argparser.loud();
	else
		argparser.mute();

	const string inputfile_name1 = argparser("-simdata1").asString("none");
	const string inputfile_name2 = argparser("-simdata2").asString("none");

//	float threshold = 0.1;	// 0.1%
	const double threshold = (double) argparser("-threshold").asDouble(1);
	if (isroot) printf("threshold = %.2lf%%\n", threshold);
//	printf("give threshold =");
//	scanf("%f", &threshold);

	if ((inputfile_name1 == "none")||(inputfile_name2 == "none"))
	{
		printf("usage: %s -simdata1 <filename>  -simdata2 <filename> [-swap] [-wtype <wtype>]\n", argv[0]);
		exit(1);
	}

	const bool swapbytes = argparser.check("-swap");
	const int wtype = argparser("-wtype").asInt(1);

	Reader_WaveletCompressionMPI  myreader1(mycomm, inputfile_name1, swapbytes, wtype);
//	Reader_WaveletCompressionMPI0 myreader2(mycomm, inputfile_name2, swapbytes, wtype);

	myreader1.load_file();
//	myreader2.load_file();
	const double init_t1 = omp_get_wtime();

	int dim[3], period[3], reorder;
	int coord[3], id;

	int NBX1 = myreader1.xblocks();
	int NBY1 = myreader1.yblocks();
	int NBZ1 = myreader1.zblocks();
	if (isroot) fprintf(stdout, "[A] I found in total %dx%dx%d blocks.\n", NBX1, NBY1, NBZ1);

	static Real targetdata1[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];

	const double t0 = omp_get_wtime();

	const int nblocks = NBX1*NBY1*NBZ1;
	const int b_end = ((nblocks + (mpi_size - 1))/ mpi_size) * mpi_size;

	for (int b = mpi_rank; b < b_end; b += mpi_size)
	{
		int z = b / (NBY1 * NBX1);
		int y = (b / NBX1) % NBY1;
		int x = b % NBX1;

		if (b < nblocks)
		{
#if defined(VERBOSE)
			fprintf(stdout, "loading block( %d, %d, %d )...\n", x, y, z);
#endif
			double zratio1 = myreader1.load_block3(x, y, z, targetdata1);
//			double zratio1 = myreader1.load_block2(x, y, z, targetdata1);

		}
		else {
		}
	}
	const double t1 = omp_get_wtime();

	MPI_Barrier(MPI_COMM_WORLD);

	if (!mpi_rank)
	{
		fprintf(stdout, "Init time = %.3lf seconds\n", init_t1-init_t0);
		fprintf(stdout, "Elapsed time = %.3lf seconds\n", t1-t0);
		myreader1.print_times();
		fflush(0);
	}

	/* Close/release resources */
	MPI_Finalize();

	return 0;
}
