/*
 *  MaxSpeedOfSound.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include "MaxSpeedOfSound.h"

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "common.h"

using namespace std;

Real MaxSpeedOfSound_CPP::compute(const Real * const src, const int gptfloats)
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

#ifdef _LIQUID_
  		const Real P = src[i+6];
#endif
        
#ifndef _LIQUID_
        const Real p = (e - (u*u + v*v + w*w)*(0.5/r))/G;
        const Real c =  mysqrt((1/G+1)*max(p/r, (Real)0));
#else
		const Real p = (e - (u*u + v*v + w*w)*(0.5/r) - P)/G;
  		const Real c =  mysqrt((1/G+1)*max((p+P/G/(1/G+1))/r, (Real)0));
#endif
		
		const Real cu = max(myabs(c + u/r), myabs(c - u/r));
		const Real cv = max(myabs(c + v/r), myabs(c - v/r));
		const Real cw = max(myabs(c + w/r), myabs(c - w/r));
		
		sos = max(sos, max( max( cu, cv), cw));
	}
	
	return sos;
}
