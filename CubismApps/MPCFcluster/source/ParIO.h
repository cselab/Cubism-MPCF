#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#pragma once

#define MAX_REPORTFREQ	512

struct MPI_ParIO
{
	MPI_Group orig_group;	// from MPI_COMM_WORLD
	int rank, size;			

	MPI_Group new_group;	// groups for gathering
	MPI_Comm new_comm;
	int new_rank, new_size;

	MPI_Group master_group;	// group of masters
	MPI_Comm master_comm;
	int master_rank, master_size;

	int groupsize;
	int reportfreq;
	int layout;	// 0: chunked, 1: ordered

	int ngroups; // = size / groupsize;

	MPI_File outfile;
	MPI_Request request;	// for asychronous writes

	string filename;
	
	float numbers[MAX_REPORTFREQ];	// vector! 
	int numid;	// vector

	int mygstep;
	long jobid;
	
	int async_counter;
	int lastcall;

	float *all_numbers;
	
	MPI_Datatype filetype, buftype;	// for ordered


	MPI_ParIO(): mygstep(0), async_counter(0), lastcall(0), all_numbers(NULL)
	{
	}

	void Init(string name, int _groupsize, int _reportfreq, int _layout)
	{
		MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
		
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		
		filename = name;

#if 0
		if (rank == 0) {
			jobid = getpid();
		}
		MPI_Bcast(&jobid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
		
		std::stringstream ss;
		ss << "." << jobid;
		filename.append(ss.str());
#endif
		
#if 1
		if (rank == 0) {
			MPI_File testfile;

			int mode = MPI_MODE_RDONLY;
			int rc = MPI_File_open( MPI_COMM_SELF, (char *)filename.c_str(), mode, MPI_INFO_NULL, &testfile);
			if (!rc) {
 				printf("Warning: Hist file %s already exists and will be deleted (error code %d).\n", (char *)filename.c_str(), rc);fflush(stdout);
				rc = MPI_File_delete((char *)filename.c_str(), MPI_INFO_NULL);
				if (rc) {
					printf("Unable to delete file %s, error code %d\n", (char *)filename.c_str(), rc);fflush(stdout);
				}
				MPI_File_close(&testfile);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		
		groupsize = _groupsize;
		reportfreq = _reportfreq;
		layout = _layout;
		all_numbers = (float *)malloc(groupsize*reportfreq*sizeof(float));	// explicit allocation of data buffer
		
		numid = 0;
		
		int range[1][3];
//		int ngroups = size / groupsize;	
//		printf("ngroups = %d\n", ngroups); fflush(0);

		int mygroup = rank / groupsize;
		int first = mygroup * groupsize; 
		int last = first + groupsize-1;
		if (last >= size) last = size-1;	// correct groupsize for the last group
	
		range[0][0] = first;
		range[0][1] = last;
		range[0][2] = 1;
	
//		printf("[%d] first = %d last = %d\n", rank, first, last); fflush(0);
		MPI_Group_range_incl(orig_group, 1, range, &new_group);  
		MPI_Group_rank (new_group, &new_rank);
		MPI_Group_size (new_group, &new_size);
		 
		MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm); 

		int range_master[1][3];
	
		range_master[0][0] = 0;
		range_master[0][1] = size-1;
		range_master[0][2] = groupsize;
	
		master_comm = MPI_COMM_NULL;
		master_rank = -1;
		master_size = 0;
	
		MPI_Group_range_incl(orig_group, 1, range_master, &master_group);

		MPI_Comm_create(MPI_COMM_WORLD, master_group, &master_comm); 
		MPI_Group_rank (master_group, &master_rank); 
		MPI_Group_size (master_group, &master_size); 

		if (new_rank == 0) {	// parallel io
//			printf("I am here [master_comm = 0x%lx]!\n", master_comm);

			/* open file */
			int mode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
			//MPI_File_open( master_comm, (char *)filename.c_str(), mode, MPI_INFO_NULL, &outfile );
			MPI_File_open( MPI_COMM_SELF, (char *)filename.c_str(), mode, MPI_INFO_NULL, &outfile );

#if 1		/* ordered, single write call */
			if ((layout == 1) || (layout == 3)) {
				/* create and commit filetype */
				MPI_Type_vector(reportfreq, groupsize, size, MPI_FLOAT, &filetype);
				MPI_Type_commit(&filetype);

				/* create and commit buftype */              
				MPI_Type_contiguous(groupsize * reportfreq, MPI_FLOAT, &buftype);
				MPI_Type_commit(&buftype);
			}
#endif

			
		}
	}

	void Notify(float number)
	{
		numbers[numid] = number;
		numid++;
	}

	void _consolidate_chunked_async(int gstep)
	{
		async_counter++;

		//if (lastcall) {
		//	printf("chunked_async: lastcall = %d\n", lastcall); fflush(0);
		//}

		int step = (gstep - reportfreq)/reportfreq;	// normalized step
		
		numid = 0;	// reset the vector ;-)
		

		if (new_rank == 0) {
			MPI_Status status;
			if (async_counter > 1) MPI_Wait(&request, &status);
		}

		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers, reportfreq, MPI_FLOAT, 0, new_comm);	// 0: master of the group (new_rank == 0)

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;
		
			size_t base = step * size * reportfreq;
			size_t offset = master_rank * groupsize * reportfreq;
			
			MPI_File_iwrite_at(outfile, (base + offset)*sizeof(float), &all_numbers[0], reportfreq*groupsize, MPI_FLOAT, &request);
			if (lastcall) {
				MPI_Wait(&request, &status);
			}
		}
		
	}

	void _consolidate_chunked(int gstep)
	{
		int step = (gstep - reportfreq)/reportfreq;	// normalized step
		
		numid = 0;	// reset the vector ;-)
		
		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers, reportfreq, MPI_FLOAT, 0, new_comm);	// 0: master of the group (new_rank == 0)

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;
		
			size_t base = step * size * reportfreq;
			size_t offset = master_rank * groupsize * reportfreq;
		
			printf("[%d %d %d %d %d %ld %ld]\n", step, size, reportfreq, master_rank, groupsize, base, offset);
			
			MPI_File_write_at(outfile, (base + offset)*sizeof(float), &all_numbers[0], reportfreq*groupsize, MPI_FLOAT, &status);
		}
		
	}

	void _consolidate_ordered(int gstep)
	{
		int step = (gstep - reportfreq)/reportfreq;	// normalized step
		
		numid = 0;	// reset the vector ;-)
		
		float all_numbers_unordered[groupsize*reportfreq];	// explicit allocation / vector?
		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers_unordered, reportfreq, MPI_FLOAT, 0, new_comm);	// 0: master of the group (new_rank == 0)

		int k = 0;
		for (int i = 0; i < reportfreq; i++) {
			for (int j = i; j < groupsize*reportfreq; j+= reportfreq) {
				all_numbers[k] = all_numbers_unordered[j];
				k++;
			} 
		}

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;

			size_t base = step * (size * reportfreq) * sizeof(float);
			size_t offset  = master_rank* groupsize * sizeof(float);
			
			MPI_File_set_view(outfile, base+offset, MPI_CHAR, filetype, (char *)"native", MPI_INFO_NULL);
			MPI_File_write_at(outfile, 0, (void *)&all_numbers[0], 1, buftype, &status);
		
//			int nbytes;
//			MPI_Get_elements(&status, MPI_CHAR, &nbytes);
//			printf("====== number of bytes written = %d ======\n", nbytes); fflush(0);

		}
		
	}

	void _consolidate_ordered_async(int gstep)
	{
		async_counter++;

		if (lastcall) {
			printf("chunked_async: lastcall = %d\n", lastcall); fflush(0);
		}

		int step = (gstep - reportfreq)/reportfreq;	// normalized step
		
		if (new_rank == 0) {
			MPI_Status status;
			if (async_counter > 1) MPI_Wait(&request, &status);
		}

		numid = 0;	// reset the vector ;-)
		
		float all_numbers_unordered[groupsize*reportfreq];	// explicit allocation / vector?
		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers_unordered, reportfreq, MPI_FLOAT, 0, new_comm);	// 0: master of the group (new_rank == 0)

		int k = 0;
		for (int i = 0; i < reportfreq; i++) {
			for (int j = i; j < groupsize*reportfreq; j+= reportfreq) {
				all_numbers[k] = all_numbers_unordered[j];
				k++;
			} 
		}

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;

			size_t base = step * (size * reportfreq) * sizeof(float);
			size_t offset  = master_rank* groupsize * sizeof(float);
			
			MPI_File_set_view(outfile, base+offset, MPI_CHAR, filetype, (char *)"native", MPI_INFO_NULL);
			
			MPI_File_iwrite_at(outfile, 0, (void *)&all_numbers[0], 1, buftype, &request);
			if (lastcall) {
				MPI_Wait(&request, &status);
			}

		
//			int nbytes;
//			MPI_Get_elements(&status, MPI_CHAR, &nbytes);
//			printf("====== number of bytes written = %d ======\n", nbytes); fflush(0);

		}
		
	}

	void Consolidate(int gstep)
	{
		mygstep = gstep;

		if (layout == 0) {
			_consolidate_chunked_async(gstep);
		}
		else if (layout == 1) {
			_consolidate_ordered(gstep);
		}/*
		else if (layout == 2) {
			_consolidate_chunked_async(gstep);
		}
		else if (layout == 3) {
			_consolidate_ordered_async(gstep);
		}*/
	}

	void Finalize()
	{
		if (numid > 0) {
			lastcall = 1;
			for (int i = numid + 1; i < MAX_REPORTFREQ; i++) numbers[i] = -1;
			Consolidate(mygstep + reportfreq);
		}
		
		if (new_rank == 0) {	// parallel io
			MPI_File_close(&outfile);
		}

		free(all_numbers);
	}

};


#if 0
#include <stdio.h>
#include "ParIO.h"

using namespace std;

int main(int argc, char *argv[])
{
	int rank, size;             

	MPI_Init(&argc,&argv);     
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int groupsize = 8;
	int reportfreq = 5;
	int layout = 3;

	MPI_ParIO pio;

	pio.Init((char *)"output.bin", groupsize, reportfreq, layout);  // 4, 2, 0

	int step;

	for (step = 0; step < 15; step++) {
		if (step% 5 == 0 && step > 0) {
			pio.Consolidate(step);
		}

		//float number = rank*10.0 + rand()%10  + step*100;
		//float number = rank;
		float number = step;
		pio.Notify(number);
	}
	
	pio.Finalize();
	MPI_Finalize();
	
	return 0;
}
#endif
