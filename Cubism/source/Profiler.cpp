/*
 *  Profiler.cpp
 *  Cubism
 *
 *  Created by Diego Rossinelli on 9/13/08.
 *  Copyright 2008 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include "Profiler.h"

#include <tbb/tick_count.h>
using namespace tbb;

void ProfileAgent::_getTime(tick_count& time)
{
	time = tick_count::now();
}

double ProfileAgent::_getElapsedTime(const tick_count& tS, const tick_count& tE)
{
	return (tE - tS).seconds();
}



