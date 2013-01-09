/*
 *  Update.cpp
 *  MPCFcore
 *
 *  Created by Babak Hejazialhosseini  on 6/9/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cassert>
#include <iostream>
#include <cstdlib>

#include "common.h"
#include "Update.h"

void Update_CPP::compute(const Real * const src, Real * const dst, const int gptfloats) const
{
	assert(gptfloats >= 7);
	
	const int N = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_ * gptfloats;
	
	for(int i=0; i<N; i+=gptfloats)
        for(int comp = 0; comp < 7; comp++)
            dst[i+comp] += m_b * src[i+comp];
}
