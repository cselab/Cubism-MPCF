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
#include <xmmintrin.h>

#include <ArgumentParser.h>

#include "TestTypes.h"
#include "Test_Convection.h"
#include "Test_Diffusion.h"
#include "Test_SurfaceTension.h"
#include "Test_LocalKernel.h"

#include "Convection_CPP.h"
#include "Update.h"
#include "MaxSpeedOfSound.h"
#include "SurfaceTension_CPP.h"
#include "Diffusion_CPP.h"
#ifdef _SSE_
#include "Convection_SSE.h"
#include "Diffusion_SSE.h"
#include "SurfaceTension_SSE.h"
#endif
#ifdef _AVX_
#include "Convection_AVX.h"
#include "SurfaceTension_AVX.h"
#include "Diffusion_AVX.h"
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
	else
		test.profile(kernel, info.peakperf*1e9, info.peakbandwidth*1e9, info.nofblocks, info.noftimes, false);
	
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
	if (parser("-f2z").asBool(false))
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	
	//kernel name
	const string kernel = parser("-kernel").asString("all");
	
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
	info.nofblocks = parser("-nblocks").asInt(512);;
	
	//number of times that the kernel will be executed
	info.noftimes =  parser("-n").asInt(1000);
	
	//C++ kernels
	{
		if (kernel == "Convection_CPP" || kernel == "all")
			testing(Test_Convection(), Convection_CPP(0, 1, 2.5, 2.1, 1, 0, 0), info);
		
		if (kernel == "Diffusion_CPP" || kernel == "all")
			testing(Test_Diffusion(), Diffusion_CPP(1, 1, 2, 1/(2.1-1), 1/(2.1-1), 1, 1, 1), info);
		
		if (kernel == "SurfaceTension_CPP" || kernel == "all")
			testing(Test_SurfaceTension(), SurfaceTension_CPP(1, 1, 1/(2.5-1), 1/(2.1-1), 1, 1, 1), info);
	}
	
	//SSE kernels
#if defined(_SSE_) && _ALIGNBYTES_ % 16 == 0
	{
		if (kernel == "Convection_SSE" || kernel == "all")
			testing(Test_Convection(), Convection_SSE(0, 1, 2.5, 2.1, 1, 0, 0), info);
		
		if (kernel == "Diffusion_SSE" || kernel == "all")
			testing(Test_Diffusion(), Diffusion_SSE(1, 1, 2, 1/(2.1-1), 1/(2.1-1), 1, 1, 1), info);
				
		if (kernel == "SurfaceTension_SSE" || kernel == "all")
			testing(Test_SurfaceTension(), SurfaceTension_SSE(1, 1, 1/(2.5-1), 1/(2.1-1), 1, 1, 1), info);
		
		if (kernel == "MaxSpeedOfSound_SSE" || kernel == "all")
			comparing(Test_LocalKernel(), MaxSpeedOfSound_CPP(2.5, 2.1, 1, 0, 0), MaxSpeedOfSound_SSE(2.5, 2.1, 1, 0, 0), info);
		
		if (kernel == "Update_SSE" | kernel == "all")
			comparing(Test_LocalKernel(), Update_CPP(), Update_SSE(), info);
	}
#endif
	
	//AVX kernels
#if defined(_AVX_) && _ALIGNBYTES_ % 32 == 0
	{
		if (kernel == "Convection_AVX" || kernel == "all")
			testing(Test_Convection(), Convection_AVX(0, 1, 2.5, 2.1, 1, 0, 0), info);
		
		if (kernel == "Diffusion_AVX" || kernel == "all")
			testing(Test_Diffusion(), Diffusion_AVX(1, 1, 2, 1/(2.1-1), 1/(2.1-1), 1, 1, 1), info);

		if (kernel == "SurfaceTension_AVX" || kernel == "all")
			testing(Test_SurfaceTension(), SurfaceTension_AVX(1, 1, 1/(2.5-1), 1/(2.1-1), 1, 1, 1), info);			
	}
#endif
	
	return 0;
}
