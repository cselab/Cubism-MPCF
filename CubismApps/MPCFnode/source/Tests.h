/*
 *  Tests.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/20/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

//#include "Test_ShockTube.h"
//typedef BlockLab<FluidBlock, std::allocator> Lab;

#include "Test_ShockBubble.h"
//typedef BlockLab<FluidBlock, std::allocator> Lab;
typedef BlockLabBubble<FluidBlock, std::allocator> Lab;

//#include "Test_SIC.h"
//maybe replace it with std::allocator
//typedef BlockLabCollapse<FluidBlock, tbb::scalable_allocator> Lab;
//typedef BlockLabCollapse<FluidBlock, std::allocator> Lab;

//typedef BlockLabBubbleYZSymmetric<FluidBlock, tbb::scalable_allocator> Lab;

//#include "Test_DoubleMachReflection.h"
//typedef BlockLabDMR<FluidBlock, tbb::scalable_allocator> Lab;

//#include "Test_CVT.h"
//typedef BlockLab_CVT<FluidBlock, tbb::scalable_allocator> Lab;
