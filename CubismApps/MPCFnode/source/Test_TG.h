/*
 *  Test_TG.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Test_SteadyState.h"

class Test_TG: public Test_SteadyState
{
    friend class Test_TGMPI;
    
	Real radius, bubble_pos[3];
    
	void _ic(FluidGrid& grid);
	void _setup_constants();
    void _dumpStatistics(FluidGrid& grid, const int counter, const Real t, const Real dt);
    
public:	
	Test_TG(const int argc, const char ** argv): Test_SteadyState(argc, argv) { }
    
	void run();
	void setup();
};
