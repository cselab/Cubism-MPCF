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
#include <hdf5.h>
#include <omp.h>
#define _TRANSPOSE_DATA_
#define _COLLECTIVE_IO_

#include "ArgumentParser.h"
#include "Reader_WaveletCompression.h"

int main(int argc, const char **argv)
{
	const double init_t0 = omp_get_wtime();

	/* Initialize MPI */
	MPI::Init_thread(MPI_THREAD_SERIALIZED);

	/* MPI variables */
	MPI_Comm comm  = MPI_COMM_WORLD;
	MPI_Info info  = MPI_INFO_NULL;
	MPI::Intracomm& mycomm = MPI::COMM_WORLD;

	const int mpi_rank = mycomm.Get_rank();
	const int mpi_size = mycomm.Get_size();

	const bool isroot = !mpi_rank;

	ArgumentParser argparser(argc, argv);

	if (isroot)
		argparser.loud();
	else
		argparser.mute();

	const string inputfile = argparser("-simdata").asString("none");
	if (inputfile == "none")
	{
		printf("usage: %s -simdata <filename>\n", argv[0]);
		exit(1);
	}

	/* HDF5 APIs definitions */
	hid_t file_id, dset_id; /* file and dataset identifiers */
	hid_t filespace, memspace;      /* file and memory dataspace identifiers */
	hsize_t dims[4]; /* dataset dimensions */
	hsize_t	count[4];	  /* hyperslab selection parameters */
	hsize_t	offset[4];
	id_t	plist_id; /* property list identifier */
	herr_t	status;

#if 1
	Reader_WaveletCompressionMPI myreader(mycomm, inputfile);
#else
	Reader_WaveletCompression myreader(inputfile);
#endif
	myreader.load_file();
	const double init_t1 = omp_get_wtime();


	const double t0 = omp_get_wtime(); 

	stringstream ss;
	ss << inputfile << "_np"  << mpi_size ;
	string h5file_name = ss.str(); 
	string h5file_fullname = h5file_name + ".h5";;

	int dim[3], period[3], reorder;
	int coord[3], id;

	/*  Set up file access property list with parallel I/O access */
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, comm, info);

	/* Create a new file collectively and release property list identifier. */
	file_id = H5Fcreate(h5file_fullname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	int NBX = myreader.xblocks();
	int NBY = myreader.yblocks();
	int NBZ = myreader.zblocks();
	printf("I found in total %dx%dx%d blocks.\n", NBX, NBY, NBZ);
	
	int NX = NBX*_BLOCKSIZE_;
	int NY = NBY*_BLOCKSIZE_;	
	int NZ = NBZ*_BLOCKSIZE_;

	/* Create the dataspace for the dataset.*/
#if defined(_TRANSPOSE_DATA_)
	dims[0] = NX;
	dims[1] = NY;
	dims[2] = NZ;
	dims[3] = 1;
#else
	dims[0] = NZ;
	dims[1] = NY;
	dims[2] = NX;
	dims[3] = 1;
#endif
	filespace = H5Screate_simple(4, dims, NULL);

	/* Create the dataset with default properties and close filespace.*/
	dset_id = H5Dcreate(file_id, "data", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(filespace);

	count[0] = _BLOCKSIZE_;
	count[1] = _BLOCKSIZE_;
	count[2] = _BLOCKSIZE_;
	count[3] = 1;
	memspace = H5Screate_simple(4, count, NULL);

	/* Create property list for collective dataset write. */
	plist_id = H5Pcreate(H5P_DATASET_XFER);
#if defined(_COLLECTIVE_IO_)
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif	
	static Real targetdata[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
	static Real storedata[_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_];

	const int nblocks = NBX*NBY*NBZ;
	const int b_end = ((nblocks + (mpi_size - 1))/ mpi_size) * mpi_size; 

	for (int b = mpi_rank; b < b_end; b += mpi_size)	
	{
#if defined(_TRANSPOSE_DATA_)
		int x = b / (NBY * NBZ);
		int y = (b / NBZ) % NBY;
		int z = b % NBZ;
#else
		int z = b / (NBY * NBX);
		int y = (b / NBX) % NBY;
		int x = b % NBX;
#endif

		if (b < nblocks)
		{
//			printf("loading block(%d,%d,%d)\n", x, y, z); 
			myreader.load_block(x, y, z, targetdata);
	
#if defined(_TRANSPOSE_DATA_)
			for (int xb = 0; xb < _BLOCKSIZE_; xb++)
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int zb = 0; zb < _BLOCKSIZE_; zb++)
						storedata[xb*_BLOCKSIZE_*_BLOCKSIZE_ + yb*_BLOCKSIZE_ + zb] = targetdata[zb][yb][xb];
			offset[0] = x * count[0];
			offset[1] = y * count[1];
			offset[2] = z * count[2];
			offset[3] = 0;
#else
			offset[0] = z * count[0];
			offset[1] = y * count[1];
			offset[2] = x * count[2];
			offset[3] = 0;
#endif
		

			/* Select hyperslab in the file */
			filespace = H5Dget_space(dset_id);
			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

#if defined(_TRANSPOSE_DATA_)
			status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, storedata);
#else
			status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, targetdata);
#endif
			H5Sclose(filespace);
		}
		else {
			H5Sselect_none(memspace);
			
			filespace = H5Dget_space(dset_id);
			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
			H5Sselect_none(filespace);
			
			status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, NULL);
		}
	}
	const double t1 = omp_get_wtime(); 

	if (!mpi_rank)
	{
		fprintf(stderr, "Init time = %.3lf seconds\n", init_t1-init_t0);
		fprintf(stderr, "Elapsed time = %.3lf seconds\n", t1-t0);
		fflush(0);
	}
	

	/* Close/release resources */
	H5Sclose(memspace);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
	H5Fclose(file_id);

	if (!mpi_rank)
	{
		char wrapper[256];
		sprintf(wrapper, "%s.xmf", h5file_name.c_str());
		FILE *xmf = 0;
		xmf = fopen(wrapper, "w");
		fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
		fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
		fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
		fprintf(xmf, " <Domain>\n");
		fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
		fprintf(xmf, "     <Time Value=\"%05d\"/>\n", 0);
		fprintf(xmf, "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n", (int)dims[0], (int)dims[1], (int)dims[2]);
		fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
		fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", 0.,0.,0.);
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
#if 1
		fprintf(xmf, "        %e %e %e\n", 1./(Real)max(dims[0],max(dims[1],dims[2])),1./(Real)max(dims[0],max(dims[1],dims[2])),1./(Real)max(dims[0],max(dims[1],dims[2])));
#else
		fprintf(xmf, "        %e %e %e\n", grid.getH(), grid.getH(), grid.getH());
#endif
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Geometry>\n");

		fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Node\">\n", "anything");
		fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (int)dims[0], (int)dims[1], (int)dims[2], 1);
		fprintf(xmf, "        %s:/data\n", h5file_fullname.c_str());
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Attribute>\n");

		fprintf(xmf, "   </Grid>\n");
		fprintf(xmf, " </Domain>\n");
		fprintf(xmf, "</Xdmf>\n");
		fclose(xmf);
	}
	
	MPI::Finalize();

	return 0;
}
