//
//  SurfaceTension.cpp
//  MPCFcore
//
//  Created by Babak Hejazialhosseini on 8/2/11.
//  Copyright 2011 ETH Zurich. All rights reserved.
//
#include <math.h>
#include <string.h>

#include "SurfaceTension_CPP.h"

void SurfaceTension_CPP::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{	
	InputSOA_ST& u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &l = ringls.ref();
	
	const Real ls_factor = G1 == G2 ? 1 : (Real)1/(G1 - G2);
	
	for(int sy=0; sy<_BLOCKSIZE_+2; sy++)
		for(int sx=0; sx<_BLOCKSIZE_+2; sx++)
		{
			AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));
			
			const int dx = sx-1;
			const int dy = sy-1;
			
			u.ref(dx, dy) = pt.u/pt.r;
			v.ref(dx, dy) = pt.v/pt.r;
			w.ref(dx, dy) = pt.w/pt.r;
			l.ref(dx, dy) = (Real)1-min(max((pt.l-G2)*ls_factor,(Real)0),(Real)1);
		}
}
