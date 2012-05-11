/*
 *  Test_ShockTube.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Test_SteadyState.h"

class Test_ShockTube: public Test_SteadyState
{
	void _ic(FluidGrid& grid);	
public:
	
	Test_ShockTube(const int argc, const char ** argv): Test_SteadyState(argc, argv) { }
	void run();
	void setup();
};
