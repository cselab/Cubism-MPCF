#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <string>
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


	MPI_ParIO()
	{
	}

	void Init(string name, int _groupsize, int _reportfreq, int _layout)
	{
		MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
		
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		
		filename = name;
		
		groupsize = _groupsize;
		reportfreq = _reportfreq;
		layout = _layout;
		
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
		}
	}

	void Notify(float number)
	{
		numbers[numid] = number;
		numid++;
	}

	void _consolidate_chunked(int gstep)
	{
		int step = (gstep - reportfreq)/reportfreq;	// normalized step
		
		numid = 0;	// reset the vector ;-)
		
		float all_numbers[groupsize*reportfreq];	// explicit allocation / vector?

		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers, reportfreq, MPI_FLOAT, 0, new_comm);	// 0: master of the group (new_rank == 0)

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;
		
			/* set file view */
			size_t base = step * size * reportfreq;
			size_t offset = master_rank * groupsize * reportfreq;
			
			MPI_File_write_at(outfile, (base + offset)*sizeof(float), &all_numbers[0], reportfreq*groupsize, MPI_FLOAT, &status);
		}
		
	}

	void _consolidate_ordered(int gstep)
	{
		int step = (gstep - reportfreq)/reportfreq;	// normalized step
		
		numid = 0;	// reset the vector ;-)
		
		float all_numbers[groupsize*reportfreq];	// explicit allocation / vector?
		MPI_Gather(&numbers[0], reportfreq, MPI_FLOAT, all_numbers, reportfreq, MPI_FLOAT, 0, new_comm);	// 0: master of the group (new_rank == 0)

		float all_numbers_ordered[groupsize*reportfreq];	// explicit allocation / vector?	 // todo: mpi_type and no buffering
		
		int k = 0;
		for (int i = 0; i < reportfreq; i++) {
			for (int j = i; j < groupsize*reportfreq; j+= reportfreq) {
				all_numbers_ordered[k] = all_numbers[j];
				k++;
			} 
		}

		if (new_rank == 0) {	// master: perform parallel io
			MPI_Status status;
		
			/* set file view */
			for (int k = 0, istep = gstep-reportfreq; istep < gstep; istep++,k++) {
				size_t base = istep * size;
				size_t offset = master_rank * groupsize;
			
				MPI_File_write_at(outfile, (base + offset)*sizeof(float), &all_numbers_ordered[k*groupsize], groupsize, MPI_FLOAT, &status);
			}
		}
		
	}


	void Consolidate(int gstep)
	{
		mygstep = gstep;

		if (layout == 0) {
			_consolidate_chunked(gstep);
		}
		else {
			_consolidate_ordered(gstep);
		}
	}

	void Finalize()
	{
		if (numid > 0) {
			for (int i = numid + 1; i < MAX_REPORTFREQ; i++) numbers[i] = -1;
			Consolidate(mygstep + reportfreq);
		}
		
		if (new_rank == 0) {	// parallel io
			MPI_File_close(&outfile);
		}
	}

};


#if 0
#include <mpi.h>
#include "ParIO.h" 

int main(int argc, char *argv[])
{
	int rank, size;

	MPI_Init(&argc,&argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size); 

	/* Extract the original group handle */ 
	/* Divide tasks into two distinct groups based upon rank */ 

	int groupsize = 4;
	int reportfreq = 2;
	int layout = 0;

	MPI_ParIO pio;

	pio.Init((char *)"output.bin", groupsize, reportfreq, layout);	// 4, 2, 0

	int step;
	
	for (step = 0; step < 10; step++) {
		if (step% 2 == 0 && step > 0) {
			pio.Consolidate(step);
		}

//		float number = rank*10.0 + rand()%10  + step*100;
//		float number = rank;
		float number = step;
		pio.Notify(step, number);
	}

	pio.Finalize();

	MPI_Finalize();
	
	return 0;
} 
#endif
