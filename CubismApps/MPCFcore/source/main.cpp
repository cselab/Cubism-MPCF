/*
 *  main.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#if defined(_QPXEMU_) && _BLOCKSIZE_%4!=0
#error BLOCKSIZE NOT GOOD FOR QPX
#endif

#include <iostream>
#include <omp.h>

#ifdef _USE_HPM_
#include <mpi.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
//extern "C" void HPM_Init(void);
//extern "C" void HPM_Print(void);
#else
#define HPM_Start(x)
#define HPM_Stop(x)
#define MPI_Init(a,b)
#define MPI_Finalize()
#endif

#include <unistd.h>
#include <ArgumentParser.h>

#include "TestTypes.h"
#include "Test_Convection.h"
#include "Test_LocalKernel.h"

#include "Convection_CPP.h"
#include "Convection_CPP_omp.h"

#if defined(_QPX_) || defined(_QPXEMU_)
#include "Convection_QPX.h"
#include "Update_QPX.h"
#include "MaxSpeedOfSound_QPX.h"
#endif

#include "Update.h"
#include "MaxSpeedOfSound.h"

using namespace std;

struct TestInfo
{
	string name;
	bool profiling;
	double peakperf, peakbandwidth;
	double accuracythreshold;
	int nofblocks, noftimes;
	
	TestInfo(string name): name(name){}
};

template<typename Test , typename Kernel> void testing(Test test, Kernel kernel, const TestInfo info)
{
	if (info.name == "all")
		printKernelName(typeid(kernel).name());
	else
		printKernelName(info.name);
	
	printKernelName(info.name);
	
	if(!info.profiling)
	{
		test.accuracy(kernel, info.accuracythreshold, false);
		printf("done!"); exit(0);
		//test.performance(kernel, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);	
	}
	else {
		if (info.name == "Convection_CPP_omp")
			test.profile2(kernel, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);
		else
			test.profile(kernel, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);
	}
	
	printEndKernelTest();
}

template<typename Test , typename Kernel1, typename Kernel2> void comparing(Test test, Kernel1 kernel1, Kernel2 kernel2, const TestInfo info)
{
	if (info.name == "all")
		printKernelName(typeid(kernel1).name());
	else
		printKernelName(info.name);
	
	test.accuracy(kernel2, kernel1, info.accuracythreshold, false);
	test.performance(kernel2, kernel1, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);
	
	printEndKernelTest();
}

int main (int argc, const char ** argv) 
{
#ifdef _USE_HPM_
	MPI_Init(&argc, const_cast<char ***>(&argv));
#endif  	
	ArgumentParser parser(argc, argv);
	parser.loud();	
	//enable/disable the handling of denormalized numbers
#ifdef _QPXEMU_
	if (parser("-f2z").asBool(false))
#pragma omp parallel
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	}
#endif	
	
	//kernel name
	const string kernel = parser("-kernel").asString("all");
	//const string kernel = parser("-kernel").asString("Convection_CPP_omp");
	
	TestInfo info(kernel);
	
	//enable/disable the performance comparison betweens kernels and their baseline
	info.profiling = parser("-profile").asBool(false);
	
	//nominal peak performance per core
	info.peakperf = parser("-pp").asDouble(21);
	
	//nominal peak bandwidth per core
	info.peakbandwidth = parser("-pb").asDouble(4.5);
	
	//numerical discrepancies greater than this threshold will be reported
	info.accuracythreshold = parser("-accuracy").asDouble(1e-4);
	
	//memory footprint per thread, in terms of blocks.
	info.nofblocks = parser("-nblocks").asInt(50);
	
	//number of times that the kernel will be executed
	info.noftimes =  parser("-n").asInt(50);
	
	//C++ kernels
	{
		if (kernel == "Convection_CPP" || kernel == "all")
			testing(Test_Convection(), Convection_CPP(0, 1), info);		
		
		//this one does not pass the accuracy test
		//if (kernel == "Convection_CPP_omp" || kernel == "all")
		//	testing(Test_Convection(), Convection_CPP_omp(0, 1), info);		
		
		{
			Test_LocalKernel lt;
			
			if (kernel == "Update_CPP" || kernel == "all")
			{
				Update_CPP update_kernel;
				HPM_Start("Update_CPP");
				lt.profile_update(update_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
				HPM_Stop("Update_CPP");
			}
			
			if (kernel == "MaxSOS_CPP" || kernel == "all")
			{
				MaxSpeedOfSound_CPP maxsos_kernel;
				HPM_Start("MaxSOS_CPP");
				lt.profile_maxsos(maxsos_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
				HPM_Stop("MaxSOS_CPP");
			}
		}
	}
	
	//QPX kernels
#if defined(_QPX_) || defined(_QPXEMU_)
	{
		if (kernel == "Convection_QPX" || kernel == "all")
			testing(Test_Convection(), Convection_QPX(0, 1), info);

		Test_LocalKernel lt;

		if (kernel == "MaxSOS_QPX" || kernel == "all")
		{
			MaxSpeedOfSound_QPX maxsos_kernel;
			MaxSpeedOfSound_CPP refkernel;
			lt.accuracy(maxsos_kernel, refkernel, info.accuracythreshold);
			HPM_Start("MaxSOS_QPX");
			lt.profile_maxsos(maxsos_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
			HPM_Stop("MaxSOS_QPX");
		}
		
		if (kernel == "Update_QPX" || kernel == "all")
		{
			Update_CPP refkernel;
			Update_QPX update_kernel;
			lt.accuracy(update_kernel, refkernel, info.accuracythreshold);
			
			HPM_Start("Update_QPX");
			lt.profile_update(update_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
			HPM_Stop("Update_QPX");
		}
	}
#endif
		
#ifdef _USE_HPM_
	MPI_Finalize();
#endif
	
	return 0;
}
