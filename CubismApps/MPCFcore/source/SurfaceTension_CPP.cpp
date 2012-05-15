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
	
	for(int sy=0; sy<_BLOCKSIZE_+2; sy++)
		for(int sx=0; sx<_BLOCKSIZE_+2; sx++)
		{
			AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));
			
			const int dx = sx-1;
			const int dy = sy-1;
			
			u.ref(dx, dy) = pt.u/pt.r;
			v.ref(dx, dy) = pt.v/pt.r;
			w.ref(dx, dy) = pt.w/pt.r;

			//project the levelset to volume fraction, same as FlowStep
			const Real x = min((Real)1, max((Real)-1, pt.l*(((Real)1)/smoothing_length)));
			const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
			const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
			l.ref(dx, dy) = 1 - (x<0 ? val_xneg : val_xpos);
		}
}
