/*
 *  main.cpp
 *  MPCFnode/prototypes/HelloLibnuma
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

//in this file i try to gain knowledge about how to map thread ids to the core ids.

#include <cstdio>
#include <iostream>
#include <numa.h>
#include <omp.h>
#include <vector>
#include <set>

using namespace std;
void print_mask(bitmask* mask, const int nodes)
{
	for(int i=0; i<nodes; ++i)
		printf("%d ", numa_bitmask_isbitset(mask, i));	
}

int main()
{
	std::cout << "Hello I am going to play with libnuma.\n";
	
	
	if (numa_available() < 0)
	{
		std::cout << "OOops NUMA API not supported! Exiting now...\n";
		return -1;
	}
	else
		std::cout << "libnuma supported!\n";
	
	const int numanodes = numa_max_node()+1;
	std::cout << "There are " << numanodes << " numa nodes.\n";
	std::cout << "We can run on " << numa_num_thread_cpus() << " cpus\n";
	std::cout << "Number of configured cpus from libnuma: "<< numa_num_configured_cpus() << endl;
	std::cout << "We print now how CPU id is mapped to NUMA node:\n";
	const int NTHREADS = omp_get_max_threads();
	for(int t=0; t<NTHREADS; ++t)
		printf("%d -> %d\n", t, numa_node_of_cpu(t));
	
	std::cout << "We launch OpenMP threads, we identify on-the-fly to which numa-node they run:\n";
#pragma omp parallel
	{
		const int mynode = numa_node_of_cpu(omp_get_thread_num());		
		//numa_run_on_node(mynode);
		bitmask * mask = numa_get_run_node_mask();
#pragma omp critical
		{	
			std::cout << "omp thread nr. " << omp_get_thread_num() << " mask: ";
			print_mask(mask, numanodes) ;
			std::cout << " and i belong to numa node "<< mynode << "\n";
		}
	}
	
	std::cout << "We launch an openmp parallel for, we figure out which thread is responsible for each loop:\n";
	{
		const int N = 16*16*16;
		vector<int> values(N);
#pragma omp parallel for
		for(int i=0; i<N; ++i)
		{
			values[i] = omp_get_thread_num();
		}
		
		set<int> already_done;
		int oldID = -1;
		for(int i=0; i<N; ++i)
		{
			const int currID = values[i];
			if (oldID != currID)
				if (already_done.find(currID) != already_done.end())
				{
					std::cout << "Associated intervals are not contiguous!\n";
					break;
				}
				else
				{
					already_done.insert(currID);
					oldID = currID;
				} 
		}
		
		vector<int> count(NTHREADS);
		for(int i=0; i<N; ++i)
			count[values[i]]++;
		
		for(int t=0; t<NTHREADS; ++t)
			printf("%d:%d ", t, count[t]);
		cout << endl;
	}
	
	cout << "Final game: we control thread-data affinity with an openmp for loop\n";
	
	cout << "We launch an openmp parallel for, we figure out which thread is responsible for each loop:\n";
	{
		const int N = 16*16*16;
		vector<int> values(N);
		vector<bitmask*> assignednodes(NTHREADS);
#pragma omp parallel
		{
			const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
			const int mynode = omp_get_thread_num() / cores_per_node;
			numa_run_on_node(mynode);
			assignednodes[omp_get_thread_num()] = numa_get_run_node_mask();
#pragma omp for
			for(int i=0; i<N; ++i)
			{
				values[i] = omp_get_thread_num();
			}
		}
		
		set<int> already_done;
		int oldID = -1;
		for(int i=0; i<N; ++i)
		{
			const int currID = values[i];
			if (oldID != currID)
				if (already_done.find(currID) != already_done.end())
				{
					std::cout << "Associated intervals are not contiguous!\n";
					break;
				}
				else
				{
					already_done.insert(currID);
					oldID = currID;
				}
		}
		
		vector<int> count(NTHREADS);
		for(int i=0; i<N; ++i)
			count[values[i]]++;
		
		for(int t=0; t<NTHREADS; ++t)
		{
			printf("thread %d -> %d elements, with mask ", t, count[t]);
			print_mask(assignednodes[t], numanodes);
			cout << endl;
		}
		cout << endl;
	}
	
	return 0;
}
