/*
 *  MaxSpeedOfSound_QPX.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 7/2/13.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#ifndef _DIEGO_TRANSPOSE4
#define _DIEGO_TRANSPOSE4(a, b, c, d)\
{\
const vector4double v01L = vec_perm(a, b, vec_gpci(00415));\
const vector4double v01H = vec_perm(a, b, vec_gpci(02637));\
const vector4double v23L = vec_perm(c, d, vec_gpci(00415));\
const vector4double v23H = vec_perm(c, d, vec_gpci(02637));\
\
a = vec_perm(v01L, v23L, vec_gpci(00145));\
b = vec_perm(v01L, v23L, vec_gpci(02367));\
c = vec_perm(v01H, v23H, vec_gpci(00145));\
d = vec_perm(v01H, v23H, vec_gpci(02367));\
}
#endif

#include "MaxSpeedOfSound.h"

class MaxSpeedOfSound_QPX : public MaxSpeedOfSound_CPP
{	
	enum 
	{ 
		NPOINTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_,
		JUMP = sizeof(Real) * 4
	};
	
	inline vector4double _compute_p(const vector4double invr, const vector4double e, 
									const vector4double invG, const vector4double P,
									const vector4double speed2) const
	{
		const vector4double tmp = vec_sub(e, P);
		return vec_mul(invG, vec_madd(vec_mul(vec_splats(-0.5f), invr), speed2, tmp));
	}
	
	inline vector4double _compute_c(const vector4double invr, const vector4double p, 
									const vector4double invG, const vector4double P) const
	{
		vector4double tmp = vec_madd(invG, vec_add(p, P), p);
		return mysqrt<preclevel>(vec_mul(invr, tmp));
	}
	
	template<int GPTFLOATS> inline vector4double _sweep4(Real * const src) const
	{
		enum { PTJUMP = sizeof(Real) * GPTFLOATS };
		
		vector4double data0 = vec_lda(0L, src);			//rho
		vector4double data1 = vec_lda(PTJUMP, src);		//u
		vector4double data2 = vec_lda(PTJUMP * 2, src);	//v
		vector4double data3 = vec_lda(PTJUMP * 3, src);	//w
		
		_DIEGO_TRANSPOSE4(data0, data1, data2, data3);
		
		const vector4double invr = myreciprocal<preclevel>(data0);
		const vector4double speed2 = vec_madd(data1, data1, vec_madd(data2, data2, vec_mul(data3, data3)));
		const vector4double maxvel = mymax(vec_abs(data1), mymax(vec_abs(data2), vec_abs(data3))) ;
		
		data0 = vec_lda(JUMP, src);						//s
		data1 = vec_lda(PTJUMP + JUMP, src);			//G
		data2 = vec_lda(PTJUMP * 2 + JUMP, src);		//P
		data3 = vec_lda(PTJUMP * 3 + JUMP, src);		//mickey mouse
		
		_DIEGO_TRANSPOSE4(data0, data1, data2, data3);
		
		const vector4double invG = myreciprocal<preclevel>(data1);
		const vector4double p = _compute_p(invr, data0, invG, data2, speed2);
		const vector4double c = _compute_c(invr, p, invG, data2);
		
		return vec_madd(maxvel, invr, c);
	}
	
public:
	
	template<int GPTFLOATS> Real _compute(Real * const src) const
	{
		enum { NFLOATS = GPTFLOATS * NPOINTS };
		
		vector4double sos4 = vec_splats(0);
		
		for(int i=0; i < NFLOATS; i += 4 * GPTFLOATS)
			sos4 = mymax(sos4, _sweep4<GPTFLOATS>(src + i));
		
		sos4 = mymax(sos4, vec_perm(sos4, sos4, vec_gpci(2323)));
		
		/* why does this not work?
		 * sos4 = mymax(sos4, vec_perm(sos4, sos4, vec_gpci(1111)));
		 * return vec_extract(sos4, 0);
		 */
		
		return std::max(vec_extract(sos4, 0), vec_extract(sos4, 1));
	}
	
	Real compute(const Real * const _src, const int gptfloats) const
	{
		assert(gptfloats == 8 || gptfloats == 16);
		
		Real * const src = const_cast<Real *>(_src);	
		
		if (gptfloats == 8)
			return _compute<8>(src);
		else if (gptfloats == 16)
			return _compute<16>(src);
		else 
		{
			printf("ooops MaxSpeedOfSound_QPX::compute: gptfloats is not quite right. aborting.\n");
			abort();
			return -1;		
		}
	}
};
