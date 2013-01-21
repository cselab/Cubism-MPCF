/*
 *  Test_1D.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 1/18/13.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Test_ShockBubble.h"

class Test_1D: public Test_ShockBubble
{
    
    void _ic(FluidGrid& grid);
      
public:	
	Test_1D(const int argc, const char ** argv): Test_ShockBubble(argc, argv) { }
    
    void setup();
};
