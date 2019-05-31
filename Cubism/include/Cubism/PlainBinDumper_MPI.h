/*
 *  PlainBinDumper_MPI.h
 *  Cubism
 *
 *  Created by Panos Hadjidoukas on 9/30/16.
 *  Copyright 2016 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once
#define MAX_MPI_PROCS	(16*1024)	// header: 0.25MB

#include <iostream>
#include <string>
#include <sstream>
#include <mpi.h>

#include "Common.h"

CUBISM_NAMESPACE_BEGIN

typedef struct _header
{
    long offset;
    long size;
} header;

#define DBG 0

template <typename TReal>
void PlainDumpBin_MPI(MPI_Comm comm, TReal *buffer, long bytes, const std::string& f_name, const std::string& dump_path=".")
{
	int rank, nranks;
	MPI_Status status;
	MPI_File file_id;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    if (rank==0)
    {
        std::cout << "Writing BIN file for ICs\n";
    }

    std::ostringstream filename;
    filename << dump_path << "/" << f_name;

	int rc = MPI_File_open( MPI_COMM_SELF, const_cast<char*>( (filename.str()+".bin").c_str() ), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_id );
	if (rc) {
		printf("Unable to create BIN file\n");
		exit(1);
	}

    long offset; // global offset
    MPI_Exscan(&bytes, &offset, 1, MPI_LONG, MPI_SUM, comm);
    if (rank == 0) offset = 0;

#if DBG
    printf("rank %d, offset = %ld, size = %ld\n", rank, offset, bytes); fflush(0);
#endif

    header tag;
    tag.offset = offset;
    tag.size = bytes;

#if DBG
    printf("DUMP IC HEADER(%d): %ld %ld\n", rank, tag.offset, tag.size);
#endif

    long base = MAX_MPI_PROCS*sizeof(tag); 	// full size of Header

    MPI_File_write_at(file_id, rank*sizeof(tag), &tag, 2, MPI_LONG, &status);
    MPI_File_write_at(file_id, base + tag.offset, (char *)buffer, tag.size, MPI_CHAR, &status);

    MPI_File_close(&file_id);
}


template <typename TReal>
void PlainReadBin_MPI(MPI_Comm comm, TReal **buffer, long *bytes, const std::string& f_name, const std::string& read_path=".")
{
	int rank, nranks;
	MPI_Status status;
	MPI_File file_id;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    std::ostringstream filename;
    filename << read_path << "/" << f_name;

	int rc = MPI_File_open( MPI_COMM_SELF, const_cast<char*>( (filename.str()+".bin").c_str() ), MPI_MODE_RDONLY, MPI_INFO_NULL, &file_id );
	if (rc) {
		printf("Unable to read ZBIN file\n");
		exit(1);
	}

    header tag;
    tag.offset = -1;
    tag.size = -1;

    MPI_File_read_at(file_id, rank*sizeof(tag), &tag, 2, MPI_LONG, &status);

#if DBG
    printf("READ IC HEADER(%d): %ld %ld\n", rank, tag.offset, tag.size); fflush(0);
#endif

    long offset = tag.offset;
    long lbytes = tag.size;
#if DBG
    printf("Reading %ld bytes of binary data\n", lbytes);
#endif
    unsigned char *tmp = (unsigned char *) malloc(lbytes);

    long base = MAX_MPI_PROCS*sizeof(tag);	// Header
    MPI_File_read_at(file_id, base + offset, (char *)tmp, lbytes, MPI_CHAR, &status);

    *buffer = (TReal *)tmp;
    *bytes = lbytes;

    MPI_File_close(&file_id);
}

CUBISM_NAMESPACE_END
