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

#ifdef __INTEL_COMPILER
inline __m128 operator+(__m128 a, __m128 b){ return _mm_add_ps(a, b); }
inline __m128 operator|(__m128 a, __m128 b){ return _mm_or_ps(a, b); }
inline __m128 operator*(__m128 a, __m128 b){ return _mm_mul_ps(a, b); }
inline __m128 operator-(__m128 a,  __m128 b){ return _mm_sub_ps(a, b); }
inline __m128 operator/(__m128 a,  __m128 b){ return _mm_div_ps(a, b); }
#endif

using namespace std;

template<typename X> inline X mysqrt(X x){ abort(); return sqrt(x);}
template<>  inline float mysqrt<float>(float x){ return sqrtf(x);}
template<>  inline double mysqrt<double>(double x){ return sqrt(x);}

template<typename X> inline X myabs(X x){ abort(); return sqrt(x);}
template<>  inline float myabs<float>(float x){ return fabsf(x);}
template<>  inline double myabs<double>(double x){ return fabs(x);}

Real MaxSpeedOfSound_CPP::_getgamma(const Real phi)
{
	const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/m_smoothlength)));
	const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
	const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
	const Real hs = x<0 ? val_xneg : val_xpos;
	
	return m_gamma1*hs + m_gamma2*(((Real)1)-hs);
}

Real MaxSpeedOfSound_CPP::_getPC(const Real phi)
{
	const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/m_smoothlength)));
	const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
	const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
	const Real hs = x<0 ? val_xneg : val_xpos;
	
	return m_pc1*hs + m_pc2*(((Real)1)-hs);
}

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
		const Real l = src[i+5];
		
		const Real p = (e - (u*u + v*v + w*w)*(0.5/r))*(_getgamma(l)-((Real)1)) -_getgamma(l)*_getPC(l);
		const Real c =  mysqrt(_getgamma(l)* max((p+_getPC(l))*((Real)1/r), (Real)0));
		
		const Real cu = max(myabs(c + u/r), myabs(c - u/r));
		const Real cv = max(myabs(c + v/r), myabs(c - v/r));
		const Real cw = max(myabs(c + w/r), myabs(c - w/r));
		
		sos = max(sos, max( max( cu, cv), cw));
	}
	
	return sos;
}

#ifdef _SSE_
inline __m128 better_rcp(const __m128 a)
{
	const __m128 Ra0 = _mm_rcp_ps(a);
	return _mm_sub_ps(_mm_add_ps(Ra0, Ra0), _mm_mul_ps(_mm_mul_ps(Ra0, a), Ra0));
}

inline __m128 worse_sqrt(const __m128 a)
{
	const __m128 invz =  _mm_rsqrt_ps(a);
	const __m128 z = _mm_rcp_ps(invz);
	return z-(z*z-a)*invz*_mm_set_ps1(0.5f);
}

inline __m128 MaxSpeedOfSound_SSE::_heaviside(const __m128 phi, const __m128 inv_h, const __m128 one, const __m128 phalf, const __m128 mhalf) const
{
	const __m128 x = _mm_min_ps(one, _mm_max_ps(_mm_setzero_ps() - one, phi*inv_h));
	
	const __m128 val_xneg = (mhalf*x - one)*x + phalf;
	const __m128 val_xpos = (phalf*x - one)*x + phalf;
	
	const __m128 flag = _mm_cmplt_ps(x, _mm_setzero_ps());
	
	return _mm_or_ps(_mm_and_ps(flag, val_xneg),_mm_andnot_ps(flag, val_xpos));
}

inline __m128 MaxSpeedOfSound_SSE::_getgamma(const __m128 phi, const __m128 inv_smoothlength, 
											 const __m128 gamma1, const __m128 gamma2,
											 const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const
{
	const __m128 hs = _heaviside(phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	return (gamma1)*hs + (gamma2)*(F_1-hs); 
}

inline __m128 MaxSpeedOfSound_SSE::_getPC(const __m128 phi, const __m128 inv_smoothlength, 
										  const __m128 pc1, const __m128 pc2,
										  const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const
{
	const __m128 hs = _heaviside(phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	return (pc1)*hs + (pc2)*(F_1-hs); 
}

void MaxSpeedOfSound_SSE::_sse_convert(const float * const gptfirst, const int gptfloats, float * const r, 
									   float * const u, float * const v, float * const w, float * const p, float * const l)
{	
	const __m128 F_1_2 = _mm_set_ps1(0.5f);
	const __m128 M_1_2 = _mm_set_ps1(-0.5f);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/m_smoothlength);
	const __m128 g1 = _mm_set_ps1(m_gamma1);
	const __m128 g2 = _mm_set_ps1(m_gamma2);
	const __m128 pc1 = _mm_set_ps1(m_pc1);
	const __m128 pc2 = _mm_set_ps1(m_pc2);
	
#define ID (dx + _BLOCKSIZE_*dy)
	
	for(int dy=0; dy<_BLOCKSIZE_; dy++)
	{		
		for(int dx=0; dx<_BLOCKSIZE_; dx+=4)
		{
			__m128 dataA0 = _mm_loadu_ps(gptfirst + (ID + 0)*gptfloats);
			__m128 dataA1 = _mm_loadu_ps(gptfirst + (ID + 0)*gptfloats + 4);
			__m128 dataB0 = _mm_loadu_ps(gptfirst + (ID + 1)*gptfloats);
			__m128 dataB1 = _mm_loadu_ps(gptfirst + (ID + 1)*gptfloats + 4);
			__m128 dataC0 = _mm_loadu_ps(gptfirst + (ID + 2)*gptfloats);
			__m128 dataC1 = _mm_loadu_ps(gptfirst + (ID + 2)*gptfloats + 4);
			__m128 dataD0 = _mm_loadu_ps(gptfirst + (ID + 3)*gptfloats );
			__m128 dataD1 = _mm_loadu_ps(gptfirst + (ID + 3)*gptfloats + 4);
			
			_MM_TRANSPOSE4_PS(dataA0, dataB0, dataC0, dataD0);
			
			_mm_store_ps(r + ID, dataA0);
			
#ifdef _PREC_DIV_
			const __m128 inv_rho = F_1/dataA0;
#else
			const __m128 inv_rho = F_1*better_rcp(dataA0);
#endif
			_mm_store_ps(u + ID, dataB0*inv_rho);
			_mm_store_ps(v + ID, dataC0*inv_rho);
			_mm_store_ps(w + ID, dataD0*inv_rho);
			
			_MM_TRANSPOSE4_PS(dataA1, dataB1, dataC1, dataD1);
			
			_mm_store_ps(l + ID, dataB1);
			
			_mm_store_ps(p + ID,  
						 (dataA1 - (dataB0*dataB0 + dataC0*dataC0 + dataD0*dataD0)*(F_1_2*inv_rho))*
						 (_getgamma(dataB1, invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1) -
                         _getgamma(dataB1, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(dataB1, invsml, pc1, pc2, F_1, F_1_2, M_1_2));
		}
	}
	
#undef ID
}

void MaxSpeedOfSound_SSE::_sse_maxsos(const float * const r, const float * const  u, const float * const  v, 
									  const float * const w, const float * const  p, const float * const l, float * const sos)
{
	const __m128 F_1_2 = _mm_set_ps1(0.5);
	const __m128 M_1_2 = _mm_set_ps1(-0.5);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/m_smoothlength);
	const __m128 g1 = _mm_set_ps1(m_gamma1);
	const __m128 g2 = _mm_set_ps1(m_gamma2);
	const __m128 pc1 = _mm_set_ps1(m_pc1);
	const __m128 pc2 = _mm_set_ps1(m_pc2);
    
	const unsigned int absvalmask1 = 0x7fffffff;
	const float * const absmaskptr = (float *)&absvalmask1;
	const __m128 absmask = _mm_set1_ps(*absmaskptr);
	
	__m128 maxsos = _mm_setzero_ps();
	
	const int N = _BLOCKSIZE_*_BLOCKSIZE_;	
	for(int i=0; i<N; i+=4)
	{
#ifdef _PREC_DIV_
		const __m128 x = _getgamma(_mm_load_ps(l + i), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
		_mm_max_ps((_mm_load_ps(p + i) +  _getPC(_mm_load_ps(l + i), invsml, pc1, pc2, F_1, F_1_2, M_1_2))*
				   better_rcp(_mm_load_ps(r + i)), _mm_setzero_ps());
		
		const __m128 tmp = worse_sqrt(x);
		
		const __m128 c = _mm_and_ps(tmp, _mm_cmpgt_ps(x, _mm_setzero_ps()));
#else
		const __m128 c = _mm_sqrt_ps(_getgamma(_mm_load_ps(l + i), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
									 _mm_max_ps((_mm_load_ps(p + i)+_getPC(_mm_load_ps(l + i), invsml, pc1, pc2, F_1, F_1_2, M_1_2))/(_mm_load_ps(r + i)), _mm_setzero_ps()));
#endif
		
		const __m128 myu = _mm_load_ps(u + i);
		const __m128 myv = _mm_load_ps(v + i);
		const __m128 myw = _mm_load_ps(w + i);
		
#define cu  _mm_max_ps(_mm_and_ps(c + myu, absmask), _mm_and_ps(c - myu, absmask))
#define cv  _mm_max_ps(_mm_and_ps(c + myv, absmask), _mm_and_ps(c - myv, absmask))
#define cw  _mm_max_ps(_mm_and_ps(c + myw, absmask), _mm_and_ps(c - myw, absmask))
		
		maxsos = _mm_max_ps(_mm_max_ps(cu, cv), _mm_max_ps(cw, maxsos));
		
#undef cu
#undef cv
#undef cw
	}
	
	_mm_store_ps(sos, _mm_max_ps(_mm_load_ps(sos), maxsos));
}
#endif
