/*
 *  main.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#if defined(_SSE_) && _BLOCKSIZE_%4!=0
#error BLOCKSIZE NOT GOOD FOR SSE
#elif defined(_AVX_) && (_BLOCKSIZE_%8!=0 || !defined(_SP_COMP_))
#error BLOCKSIZE NOT GOOD FOR AVX
#endif

#include <iostream>
#include <omp.h>

#ifdef _SSE_
#include <xmmintrin.h>
#endif
#include <unistd.h>
#include <ArgumentParser.h>

#include "TestTypes.h"
#include "Test_Convection.h"
#include "Test_LocalKernel.h"

#include "Convection_CPP.h"
#include "Convection_CPP_omp.h"
#include "Convection_QPX.h"

#include "Update.h"
#include "MaxSpeedOfSound.h"
#ifdef _SSE_
#include "Convection_SSE.h"
#endif
#ifdef _AVX_
#include "Convection_AVX.h"
#endif

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
		test.performance(kernel, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);	
	}
	else {
		//if (info.name == "Convection_CPP_omp")
		//	test.profile2(kernel, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);
		//else
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
	ArgumentParser parser(argc, argv);

	//enable/disable the handling of denormalized numbers
#ifdef _SSE_
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

                if (kernel == "Convection_QPX" || kernel == "all")
                        testing(Test_Convection(), Convection_QPX(0, 1), info);

		if (kernel == "Convection_CPP_omp" || kernel == "all")
			testing(Test_Convection(), Convection_CPP_omp(0, 1), info);		


		if (kernel == "Update_CPP" || kernel == "all")
		{
			Test_LocalKernel lt;
			Update_CPP update_kernel;
			lt.profile_update(update_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);

			//sleep(100);
		}

		if (kernel == "MaxSOS_CPP" || kernel == "all")
		{
			Test_LocalKernel lt;
			MaxSpeedOfSound_CPP maxsos_kernel;
			lt.profile_maxsos(maxsos_kernel, info.peakperf, info.peakbandwidth, info.nofblocks, info.noftimes);
		}

	}
	
	//SSE kernels
#if defined(_SSE_) && _ALIGNBYTES_ % 16 == 0
	{
		if (kernel == "Convection_SSE" || kernel == "all")
			testing(Test_Convection(), Convection_SSE(0, 1), info);
	}
#endif
	
	//AVX kernels
#if defined(_AVX_) && _ALIGNBYTES_ % 32 == 0
	{
		if (kernel == "Convection_AVX" || kernel == "all")
			testing(Test_Convection(), Convection_AVX(0, 1, 2.5, 2.1, 1, 0, 0), info);
		
	}
#endif
	
	return 0;
}
