#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <time.h>

#include <string>
#include <sstream>
#include <iomanip>


using namespace std;

#pragma once

struct MPI_Streamer
{
	MPI_Group orig_group;	// from MPI_COMM_WORLD
	int rank, size;			

	int mygroup, groupsize, mymaster;
	int first, last;

	MPI_Group new_group;	// groups with common output file
	MPI_Comm new_comm;
	int new_rank, new_size;

	int *counter_mem;
	MPI_Win win;


	string filename;
	char header[512];

	MPI_File fh;
	int mode;
	MPI_Offset base, offset;
	MPI_Offset old_offset;

	MPI_Status status[3];	
	MPI_Request req[3];
	int wbytes[3];


	MPI_Streamer(int _groupsize=-1)
	{
		groupsize = _groupsize;

		MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);

		if ((groupsize == -1)||(groupsize > size)) {
			groupsize = size;

			mygroup = rank / groupsize;
			mymaster = mygroup * groupsize; 
			first = mygroup * groupsize; 
			last = first + groupsize-1;
			if (last >= size) last = size-1;	// correct groupsize for the last group

			new_comm = MPI_COMM_WORLD;
			new_rank = rank;
			new_size = size;
		}
		else {
			// groupsize already set
		
			mygroup = rank / groupsize;
			mymaster = mygroup * groupsize; 
			first = mygroup * groupsize; 
			last = first + groupsize-1;
			if (last >= size) last = size-1;	// correct groupsize for the last group

#if 1
			new_comm = MPI_COMM_WORLD;
			new_rank = rank;
			new_size = size;
#else
			int range[1][3];
//		int ngroups = size / groupsize;	
//		printf("ngroups = %d\n", ngroups); fflush(0);

			range[0][0] = first;
			range[0][1] = last;
			range[0][2] = 1;
	
//			printf("[%d] first = %d last = %d\n", rank, first, last); fflush(0);
			MPI_Group_range_incl(orig_group, 1, range, &new_group);  
			MPI_Group_rank (new_group, &new_rank);
			MPI_Group_size (new_group, &new_size);
		 
			MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm); 
#endif
		}
		
//		if (new_rank == mymaster) {
		MPI_Alloc_mem(4, MPI_INFO_NULL, &counter_mem);
		MPI_Win_create(counter_mem, 1*sizeof(int), sizeof(int), MPI_INFO_NULL, new_comm, &win);
//		}
//		else {
//		MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win); 
//		}

	}

	~MPI_Streamer()
	{
		if (new_rank == mymaster) {
			printf("*counter_mem = %d\n", *counter_mem);
		}
		MPI_Free_mem(counter_mem);
		MPI_Win_free(&win);
	}

	void SetHeader(char * _header)
	{
		strcpy(header, _header);
	}

	void Init(string _filename)
	{
		*counter_mem = 0;	/* only for the master */
		old_offset = -1;

		filename = _filename;

		stringstream ss;
//		ss << "." << mygroup;
//		ss << "." << first << "_" << last << ".dat" ;
		ss << "." << setfill('0') << setw(5) << first << "_" << setfill('0') << setw(5) << last << ".dat" ;
		filename.append(ss.str());

		/* open file */
		mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;

		//MPI_File_open(MPI_COMM_WORLD, filename.c_str(), mode, MPI_INFO_NULL, &fh);
		MPI_File_open(MPI_COMM_SELF, filename.c_str(), mode, MPI_INFO_NULL, &fh);

		MPI_Barrier(MPI_COMM_WORLD);

		if (new_rank == mymaster) {
			MPI_Status tmp_status;
			printf("header:\n%s\n", header);
			MPI_File_write_at(fh, 0, (void *)header, 512, MPI_CHAR, &tmp_status);
			*counter_mem += 512;
			
			int firstoffset[groupsize];
			for (int i = 0; i < groupsize; i++) firstoffset[i] = -i;
			MPI_File_write_at(fh, 512, (void *)firstoffset, groupsize*sizeof(int), MPI_CHAR, &tmp_status);
			*counter_mem += groupsize*sizeof(int);
		}
	}

	void WaitResolve(char *buf)
	{
		/* Wait */
	}

	void Wait(char *buf)
	{
		MPI_Wait(&req[0], &status[0]);
		MPI_Wait(&req[1], &status[1]);
		MPI_Wait(&req[2], &status[2]);

//		if (req[0] != MPI_REQUEST_NULL) {
		MPI_Get_elements(&status[0], MPI_CHAR, &wbytes[0]);
//		}
//		else {
//			wbytes[0] = 0;
//		}
		MPI_Get_elements(&status[1], MPI_CHAR, &wbytes[1]);
		MPI_Get_elements(&status[2], MPI_CHAR, &wbytes[2]);
		printf( "TASK %d ====== number of bytes written = %d ======\n", rank, wbytes[0]+wbytes[1]+wbytes[2]);
	}
	
	void iWrite(char *buf, int nbytes)
	{
#define TAGSIZE (2*sizeof(int))
		static int tag[2];
		tag[0] = nbytes;
		tag[1] = -1;

		static int base_offset;
		base_offset = fetch_and_add(win, rank, size, nbytes+TAGSIZE);

		printf("Rank %d, base_offset %d old_offset = %d\n", rank, base_offset, old_offset); fflush(0);

		/*  write buffer to file */
		if (old_offset > -1) {
			MPI_File_iwrite_at(fh, old_offset+1*sizeof(int), (void *)&base_offset, sizeof(int), MPI_CHAR, &req[0]);
		}
		else {
			MPI_File_iwrite_at(fh, 512+new_rank*sizeof(int), (void *)&base_offset, sizeof(int), MPI_CHAR, &req[0]);
		}

		old_offset = base_offset;

		MPI_File_iwrite_at(fh, base_offset, (void *)tag, TAGSIZE, MPI_CHAR, &req[1]);
		MPI_File_iwrite_at(fh, base_offset+TAGSIZE, (void *)buf, nbytes, MPI_CHAR, &req[2]);

#undef TAGSIZE
	}


	void Write(char *buf, int nbytes)
	{
#define TAGSIZE (2*sizeof(int))
		int tag[2];
		tag[0] = nbytes;
		tag[1] = -1;

		int base_offset = fetch_and_add(win, rank, size, nbytes+TAGSIZE);

		printf("Rank %d, base_offset %d\n", rank, base_offset); fflush(0);

		/*  write buffer to file */
		if (old_offset > -1) {
			MPI_File_write_at(fh, old_offset+1*sizeof(int), (void *)&base_offset, sizeof(int), MPI_CHAR, &status[0]);
		}
		else {
			MPI_File_write_at(fh, 512+new_rank*sizeof(int), (void *)&base_offset, sizeof(int), MPI_CHAR, &status[0]);
		}
		MPI_Get_elements(&status[0], MPI_CHAR, &wbytes[0]);
		
		old_offset = base_offset;

		MPI_File_write_at(fh, base_offset, (void *)tag, TAGSIZE, MPI_CHAR, &status[1]);
		MPI_File_write_at(fh, base_offset+TAGSIZE, (void *)buf, nbytes, MPI_CHAR, &status[2]);

		/* print out number of bytes written */
		MPI_Get_elements(&status[1], MPI_CHAR, &wbytes[1]);
		MPI_Get_elements(&status[2], MPI_CHAR, &wbytes[2]);
		printf( "TASK %d ====== number of bytes written = %d ======\n", rank, wbytes[0]+wbytes[1]+wbytes[2]);
#undef TAGSIZE
	}


	void Flush()	/* Finalize */
	{
		/* close file */
		MPI_File_close(&fh);
	}

	int fetch_and_add(MPI_Win win, int rank, int nprocs, int num) 
	{
		int one = num;

		int value = 0;
	
		MPI_Win_lock(MPI_LOCK_EXCLUSIVE, mymaster, 0, win);
		MPI_Get(&value, 1, MPI_INT, mymaster, 0, 1, MPI_INT, win); 
		MPI_Accumulate(&one, 1, MPI_INT, mymaster, 0, 1, MPI_INT, MPI_SUM, win);
		MPI_Win_unlock(mymaster, win);

		return value;
	}


};

