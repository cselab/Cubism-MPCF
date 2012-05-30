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
    
	void _ic(FluidGrid& grid);
    
public:	
	Test_SIC(const int argc, const char ** argv): Test_ShockBubble(argc, argv) { }
};
