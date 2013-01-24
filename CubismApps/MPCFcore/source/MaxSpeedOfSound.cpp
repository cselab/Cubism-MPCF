/*
 *  MaxSpeedOfSound.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "common.h"
#include "MaxSpeedOfSound.h"

Real MaxSpeedOfSound_CPP::compute(const Real * const src, const int gptfloats) const
{
	const int N=_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	Real sos = 0;
	
	for(int i=0; i<N; i+=gptfloats)
	{
		const Real r = src[i];
		const Real u = src[i+1];
		const Real v = src[i+2];
		const Real w = src[i+3];
		const Real e = src[i+4];
		const Real G = src[i+5];
  		const Real P = src[i+6];
        
		const Real p = (e - (u*u + v*v + w*w)*(0.5/r) - P)/G;
  		const Real c =  std::sqrt((1/G+1)*std::max((p+P/G/(1/G+1))/r, (Real)0));
		
		const Real cu = std::abs(u/r)+c;//std::max(std::abs(c + u/r), std::abs(c - u/r));
		const Real cv = std::abs(v/r)+c;//std::max(std::abs(c + v/r), std::abs(c - v/r));
		const Real cw = std::abs(w/r)+c;//std::max(std::abs(c + w/r), std::abs(c - w/r));
		
		sos = std::max(sos, std::max(std::max(cu, cv), cw));
	}
	
	return sos;
}
