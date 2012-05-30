/*
 *  Test_SIC.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 5/30/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Test_ShockBubble.h"

class Test_SIC: public Test_ShockBubble
{
    friend class Test_SICMPI;
    
	//Real radius, bubble_pos[3];
    
	void _ic(FluidGrid& grid);
	//void _setup_constants();
    //void _dumpStatistics(FluidGrid& grid, const int counter, const Real t, const Real dt);
    
public:	
	Test_SIC(const int argc, const char ** argv): Test_ShockBubble(argc, argv) { }
    
	//void run();
	//void setup();
};
