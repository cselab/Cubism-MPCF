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

#include "ArgumentParser.h"
#include "Reader_WaveletCompression.h"

int main(int argc, const char **argv)
{
	ArgumentParser parser(argc, argv);
	parser.loud();

	const string inputfile = parser("-input").asString("none");
	if (inputfile == "none")
	{
		printf("usage: %s -input <filename>\n", argv[0]);
		exit(1);
	}

	/* HDF5 APIs definitions */ 	
	hid_t file_id, dset_id; /* file and dataset identifiers */
	hid_t filespace, memspace;      /* file and memory dataspace identifiers */
	hsize_t dims[3]; /* dataset dimensions */
	hsize_t	count[3];	  /* hyperslab selection parameters */
	hsize_t	offset[3];
	id_t	plist_id; /* property list identifier */
	herr_t	status;

	/* MPI variables */
	MPI_Comm comm  = MPI_COMM_WORLD;
	MPI_Info info  = MPI_INFO_NULL;

	/* Initialize MPI */
	MPI::Init();

	MPI::Intracomm& mycomm = MPI::COMM_WORLD;

	const int mpi_rank = mycomm.Get_rank();
	const int mpi_size = mycomm.Get_size();

#if 1
	Reader_WaveletCompressionMPI myreader(mycomm, inputfile);
#else
	Reader_WaveletCompression myreader;
#endif
//	myreader.load_file(inputfile);
	myreader.load_file();

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
	dims[0] = NX;
	dims[1] = NY;
	dims[2] = NZ;
	filespace = H5Screate_simple(3, dims, NULL);

	/* Create the dataset with default properties and close filespace.*/
	dset_id = H5Dcreate(file_id, "data", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(filespace);

	count[0] = _BLOCKSIZE_;
	count[1] = _BLOCKSIZE_;
	count[2] = _BLOCKSIZE_;
	memspace = H5Screate_simple(3, count, NULL);

	/* Create property list for collective dataset write. */
#if 1
	plist_id = H5Pcreate(H5P_DATASET_XFER);
#else
	hsize_t cdims[3];
	cdims[0] = _BLOCKSIZE_/2;
	cdims[1] = _BLOCKSIZE_/2;
	cdims[2] = _BLOCKSIZE_/2;
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	status = H5Pset_chunk (plist_id, 3, cdims);
	status = H5Pset_deflate (plist_id, 6); 
#endif

	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	
	static Real targetdata[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
	static Real storedata[_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_];

	const int nblocks = NBX*NBY*NBZ;
	const int b_end = ((nblocks + (mpi_size - 1))/ mpi_size) * mpi_size; 

	const double t0 = omp_get_wtime(); 
	for (int b = mpi_rank; b < b_end; b += mpi_size)	
	{
		int x = b / (NBY * NBZ);
		int y = (b / NBZ) % NBY;
		int z = b % NBZ;; 
	
		if (b < nblocks)
		{
//			printf("loading block(%d,%d,%d)\n", x, y, z); 
			myreader.load_block(x, y, z, targetdata);
	
#if 1
			for (int xb = 0; xb < _BLOCKSIZE_; xb++)
				for (int yb = 0; yb < _BLOCKSIZE_; yb++)
					for (int zb = 0; zb < _BLOCKSIZE_; zb++)
						storedata[xb*_BLOCKSIZE_*_BLOCKSIZE_ + yb*_BLOCKSIZE_ + zb] = targetdata[zb][yb][xb];
#endif
		
			offset[0] = x * count[0];
			offset[1] = y * count[1];
			offset[2] = z * count[2];

			/* Select hyperslab in the file */
			filespace = H5Dget_space(dset_id);
			H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

#if 1
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
		fprintf(stderr, "Elapsed time = %.3lf seconds\n", t1-t0); fflush(0);
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
		fprintf(xmf, "%e %e %e\n", 0.,0.,0.);
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "%e %e %e\n", 1., 1., 1.);
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Geometry>\n");

		fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Node\">\n", "anything");
		fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (int)dims[0], (int)dims[1], (int)dims[2], 1);
		
		fprintf(xmf, "%s:/data\n", h5file_fullname.c_str());
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
