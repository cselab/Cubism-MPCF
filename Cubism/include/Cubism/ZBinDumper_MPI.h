/*
 *  ZBinDumper_MPI.h
 *  Cubism
 *
 *  Created by Panos Hadjidoukas on 3/20/14.
 *  Copyright 2014 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <iostream>
#include <cstdio>
#include <vector>
#include <string>
#include <cassert>
#include <sstream>

#include "BlockInfo.h"
#include "LosslessCompression.h"

CUBISM_NAMESPACE_BEGIN

#define MAX_MPI_PROCS	(16*1024)	// header: 0.25MB* 8

typedef struct _header
{
    long offset[8];	// 1 for single compression, NCHANNELS otherwise
    long size[8];
} header;


// The following requirements for the data TStreamer are required:
// TStreamer::NCHANNELS        : Number of data elements (1=Scalar, 3=Vector, 9=Tensor)
// TStreamer::operate          : Data access methods for read and write
// TStreamer::getAttributeName : Attribute name of the date ("Scalar", "Vector", "Tensor")
template<typename TStreamer, typename TGrid>
void DumpZBin_MPI(
        const TGrid &grid,
        const int iCounter,
        const typename TGrid::Real t,
        const std::string &f_name,
        const std::string &dump_path = ".",
        const bool bDummy = false)
{
    typedef typename TGrid::BlockType B;

    int rank, nranks;

    // f_name is the base filename without file type extension
    std::ostringstream filename;
    filename << dump_path << "/" << f_name;

    MPI_Status status;
    MPI_File file_id;

    MPI_Comm comm = grid.getCartComm();
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    int coords[3];
    grid.peindex(coords);

    const unsigned int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
    const unsigned int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
    const unsigned int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
    static const unsigned int NCHANNELS = TStreamer::NCHANNELS;

    if (rank==0)
    {
        //Real memsize = (NX * NY * NZ * NCHANNELS * sizeof(Real))/(1024.*1024.*1024.);
        Real memsize = (NX * NY * NZ * sizeof(Real))/(1024.*1024.*1024.);
        std::cout << "Allocating " << memsize << " GB of BIN data per rank (" << memsize*nranks << " GB in total)" << std::endl;
    }
//	Real * array_all = new Real[NX * NY * NZ * NCHANNELS];
    Real * array_all = new Real[NX * NY * NZ];

    std::vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

    static const unsigned int sX = 0;
    static const unsigned int sY = 0;
    static const unsigned int sZ = 0;

    static const unsigned int eX = B::sizeX;
    static const unsigned int eY = B::sizeY;
    static const unsigned int eZ = B::sizeZ;

    int rc = MPI_File_open( MPI_COMM_SELF, const_cast<char*>( (filename.str()+".zbin").c_str() ), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_id );
    if (rc) {
        printf("Unable to create ZBIN file\n");
        exit(1);
    }


    long previous_offset = 0;
    header tag; 	// moved here

    for (unsigned int ichannel = 0; ichannel < NCHANNELS; ichannel++)
    {


#pragma omp parallel for
    for(unsigned int i=0; i<vInfo_local.size(); i++)
    {
        BlockInfo& info = vInfo_local[i];
        const unsigned int idx[3] = {info.index[0], info.index[1], info.index[2]};
        B & b = *(B*)info.ptrBlock;

        for(unsigned int ix=sX; ix<eX; ix++)
        {
            const unsigned int gx = idx[0]*B::sizeX + ix;
            for(unsigned int iy=sY; iy<eY; iy++)
            {
                const unsigned int gy = idx[1]*B::sizeY + iy;
                for(unsigned int iz=sZ; iz<eZ; iz++)
                {
                    const unsigned int gz = idx[2]*B::sizeZ + iz;

                    assert((gz + NZ * (gy + NY * gx)) < NX * NY * NZ);

                    Real * const ptr = array_all + (gz + NZ * (gy + NY * gx));

                    Real output;
                    TStreamer::operate(b, ix, iy, iz, &output, ichannel);	// point -> output,
                    ptr[0] = output;
                }
            }
        }
    }

//	long local_count = NX * NY * NZ * NCHANNELS;
    long local_count = NX * NY * NZ * 1;
    long local_bytes =  local_count * sizeof(Real);
    long offset; // global offset

    unsigned int max = local_bytes;
//	int layout[4] = {NCHANNELS, NX, NY, NZ};
    int layout[4] = {NX, NY, NZ, 1};
    long compressed_bytes = ZZcompress<typename TGrid::Real>((unsigned char *)array_all, local_bytes, layout, &max);	// "in place"
#if DBG
    printf("Writing %ld bytes of Compressed data (cr = %.2f)\n", compressed_bytes, local_bytes*1.0/compressed_bytes);
#endif
    MPI_Exscan( &compressed_bytes, &offset, 1, MPI_LONG, MPI_SUM, comm);

    if (rank == 0) offset = 0;

#if DBG
    printf("rank %d, offset = %ld, size = %ld\n", rank, offset, compressed_bytes); fflush(0);
#endif

//	header tag;
    tag.offset[ichannel] = offset + previous_offset;
    tag.size[ichannel] = compressed_bytes;
#if DBG
    printf("rank %d, offset = %ld, size = %ld\n", rank, tag.offset[ichannel], tag.size[ichannel]); fflush(0);
#endif
    previous_offset = (tag.offset[ichannel] + tag.size[ichannel]);
    MPI_Bcast(&previous_offset, 1, MPI_LONG, nranks-1, comm);

    long base = MAX_MPI_PROCS*sizeof(tag); 	// full Header

    MPI_File_write_at(file_id, base + tag.offset[ichannel], (char *)array_all, tag.size[ichannel], MPI_CHAR, &status);

    }	/* ichannel */

    MPI_File_write_at(file_id, rank*sizeof(tag), &tag, 2*8, MPI_LONG, &status);

    MPI_File_close(&file_id);
    delete [] array_all;
}

template<typename TStreamer, typename TGrid>
void ReadZBin_MPI(TGrid &grid, const std::string& f_name, const std::string& read_path=".")
{
    typedef typename TGrid::BlockType B;

    int rank, nranks;

    // f_name is the base filename without file type extension
    std::ostringstream filename;
    filename << read_path << "/" << f_name;

    MPI_Status status;
    MPI_File file_id;

    MPI_Comm comm = grid.getCartComm();
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    int coords[3];
    grid.peindex(coords);

    const int NX = grid.getResidentBlocksPerDimension(0)*B::sizeX;
    const int NY = grid.getResidentBlocksPerDimension(1)*B::sizeY;
    const int NZ = grid.getResidentBlocksPerDimension(2)*B::sizeZ;
    static const int NCHANNELS = TStreamer::NCHANNELS;

    Real * array_all = new Real[NX * NY * NZ * NCHANNELS];

    std::vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();

    static const int sX = 0;
    static const int sY = 0;
    static const int sZ = 0;

    const int eX = B::sizeX;
    const int eY = B::sizeY;
    const int eZ = B::sizeZ;

    int rc = MPI_File_open( MPI_COMM_SELF, const_cast<char*>( (filename.str()+".zbin").c_str() ), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_id );
    if (rc) {
        printf("Unable to read ZBIN file\n");
        exit(1);
    }

//	long local_count = NX * NY * NZ * NCHANNELS;
    long local_count = NX * NY * NZ * 1;
    long local_bytes = local_count * sizeof(Real);
    long offset;


    header tag;
    MPI_File_read_at(file_id, rank*sizeof(tag), &tag, 2*8, MPI_LONG, &status);

#if DBG
    printf("HEADER(%d):\n", rank);
    for (int i = 0; i < NCHANNELS; i++) {
        printf("channel %d: %ld %ld\n", i, tag.offset[i], tag.size[i]);
    }
#endif

    for (unsigned int ichannel = 0; ichannel < NCHANNELS; ichannel++)
    {

    //MPI_File_read_at(file_id, (rank*2+0)*sizeof(long), &offset     , 1, MPI_LONG, &status);
    //MPI_File_read_at(file_id, (rank*2+1)*sizeof(long), &compressed_bytes, 1, MPI_LONG, &status);
#if DBG
    printf("rank %d, offset = %ld, compr. size = %ld\n", rank, tag.offset[ichannel], tag.size[ichannel]); fflush(0);
#endif

    long compressed_bytes = tag.size[ichannel];
#if DBG
    printf("Reading %ld bytes of Compressed data (cr = %.2f)\n", compressed_bytes, local_bytes*1.0/compressed_bytes);
#endif
    unsigned char *tmp = (unsigned char *) malloc(compressed_bytes+4096);

    long base = MAX_MPI_PROCS*sizeof(tag);	// Header
    MPI_File_read_at(file_id, base + tag.offset[ichannel], (char *)tmp, compressed_bytes, MPI_CHAR, &status);

//	int layout[4] = {NCHANNELS, NX, NY, NZ};
    int layout[4] = {NX, NY, NZ, 1};
    size_t decompressed_bytes = ZZdecompress<typename TGrid::Real>(tmp, compressed_bytes, layout, (unsigned char *)array_all, local_bytes);
    free(tmp);
#if DBG
    printf("rank %d, offset = %ld, size = %ld (%ld)\n", rank, offset, decompressed_bytes, local_bytes); fflush(0);
#endif

    #pragma omp parallel for
    for(int i=0; i<vInfo_local.size(); i++)
    {
        BlockInfo& info = vInfo_local[i];
        const int idx[3] = {info.index[0], info.index[1], info.index[2]};
        B & b = *(B*)info.ptrBlock;

                for(int ix=sX; ix<eX; ix++)
          for(int iy=sY; iy<eY; iy++)
            for(int iz=sZ; iz<eZ; iz++)
              {
                    const int gx = idx[0]*B::sizeX + ix;
                    const int gy = idx[1]*B::sizeY + iy;
                    const int gz = idx[2]*B::sizeZ + iz;

                    //Real * const ptr_input = array_all + NCHANNELS*(gz + NZ * (gy + NY * gx));
                    Real * const ptr_input = array_all + (gz + NZ * (gy + NY * gx));

                    TStreamer::operate(b, *ptr_input, ix, iy, iz, ichannel);	// output -> point
                }
    }


    } /* ichannel */

    MPI_File_close(&file_id);
    delete [] array_all;
}

CUBISM_NAMESPACE_END
