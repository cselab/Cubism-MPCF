/*
 *  FlowStep_SSE_diego.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/19/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <xmmintrin.h>
#include <emmintrin.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _mm_set_pd1
#define _mm_set_pd1(a) (_mm_set_pd((a),(a)))
#endif

#define _3ORPS_(a,b,c) _mm_or_ps(a, _mm_or_ps(b,c))
#define _3ORPD_(a,b,c) _mm_or_pd(a, _mm_or_pd(b,c))

#define WENO_CONSTANTS \
const __m128 one = _mm_set_ps1(1.); \
const __m128 F_4_3 = _mm_set_ps1(4./3.); \
const __m128 F_5_3 = _mm_set_ps1(5./3.); \
const __m128 F_10_3 = _mm_set_ps1(10./3.); \
const __m128 F_11_3 = _mm_set_ps1(11./3.); \
const __m128 F_13_3 = _mm_set_ps1(13./3.); \
const __m128 F_19_3 = _mm_set_ps1(19./3.);  \
const __m128 F_25_3 = _mm_set_ps1(25./3.); \
const __m128 F_31_3 = _mm_set_ps1(31./3.); \
const __m128 mywenoeps = _mm_set_ps1(WENOEPS); \
\
const __m128 M_1_6 = _mm_set_ps1(-1./6.); \
const __m128 F_1_3 = _mm_set_ps1(1./3.); \
const __m128 F_5_6 = _mm_set_ps1(5./6.); \
const __m128 M_7_6 = _mm_set_ps1(-7./6.); \
const __m128 F_11_6 = _mm_set_ps1(11./6.); \
\
const __m128 F_1_10 = _mm_set_ps1(1./10.); \
const __m128 F_6_10 = _mm_set_ps1(6./10.); \
const __m128 F_3_10 = _mm_set_ps1(3./10.);

#define WENO_CONSTANTS_DP \
const __m128d one = _mm_set_pd1(1.); \
const __m128d F_4_3 = _mm_set_pd1(4./3.); \
const __m128d F_5_3 = _mm_set_pd1(5./3.); \
const __m128d F_10_3 = _mm_set_pd1(10./3.); \
const __m128d F_11_3 = _mm_set_pd1(11./3.); \
const __m128d F_13_3 = _mm_set_pd1(13./3.); \
const __m128d F_19_3 = _mm_set_pd1(19./3.);  \
const __m128d F_25_3 = _mm_set_pd1(25./3.); \
const __m128d F_31_3 = _mm_set_pd1(31./3.); \
const __m128d mywenoeps = _mm_set_pd1(WENOEPS); \
\
const __m128d M_1_6 = _mm_set_pd1(-1./6.); \
const __m128d F_1_3 = _mm_set_pd1(1./3.); \
const __m128d F_5_6 = _mm_set_pd1(5./6.); \
const __m128d M_7_6 = _mm_set_pd1(-7./6.); \
const __m128d F_11_6 = _mm_set_pd1(11./6.); \
\
const __m128d F_1_10 = _mm_set_pd1(1./10.); \
const __m128d F_6_10 = _mm_set_pd1(6./10.); \
const __m128d F_3_10 = _mm_set_pd1(3./10.);

#ifndef _PREC_DIV_

#ifdef _PREC_DIV_NONE_

#define WENO_OMEGAS \
const __m128 alpha0 = F_1_10*better_rcp(is0*is0); \
const __m128 alpha1 = F_6_10*better_rcp(is1*is1); \
const __m128 alpha2 = F_3_10*better_rcp(is2*is2); \
const __m128 inv_alpha = better_rcp(alpha0+alpha1+alpha2); \
const __m128 omega0=alpha0 * inv_alpha; \
const __m128 omega1=alpha1 * inv_alpha; \
const __m128 omega2= one-omega0-omega1; 

#else

#define WENO_OMEGAS \
const __m128 alpha0 = F_1_10/(is0*is0);	\
const __m128 alpha1 = F_6_10/(is1*is1);	\
const __m128 alpha2 = F_3_10/(is2*is2);		   \
const __m128 inv_alpha = better_rcp(alpha0+alpha1+alpha2); \
const __m128 omega0=alpha0 * inv_alpha; \
const __m128 omega1=alpha1 * inv_alpha; \
const __m128 omega2= one-omega0-omega1; 
#endif

#else 

#define WENO_OMEGAS \
const __m128 alpha0 = F_1_10 / (is0*is0); \
const __m128 alpha1 = F_6_10 / (is1*is1); \
const __m128 alpha2 = F_3_10 / (is2*is2); \
const __m128 inv_alpha = one / (alpha0+alpha1+alpha2); \
const __m128 omega0 = alpha0 * inv_alpha; \
const __m128 omega1 = alpha1 * inv_alpha; \
const __m128 omega2 = one - omega0 - omega1; 

#endif

#define WENO_OMEGAS_DP \
const __m128d alpha0 = F_1_10 / (is0*is0); \
const __m128d alpha1 = F_6_10 / (is1*is1); \
const __m128d alpha2 = F_3_10 / (is2*is2); \
const __m128d inv_alpha = one /(alpha0+alpha1+alpha2); \
const __m128d omega0=alpha0 * inv_alpha; \
const __m128d omega1=alpha1 * inv_alpha; \
const __m128d omega2= one - omega0 - omega1; 

#define WENO_MINUS \
const __m128 is0 = a*(a*F_4_3  - b*F_19_3 + c*F_11_3)	+ b*(b*F_25_3 - c*F_31_3)	+ c*c*F_10_3	+ mywenoeps; \
const __m128 is1 = b*(b*F_4_3  - c*F_13_3 + d*F_5_3)	+ c*(c*F_13_3 - d*F_13_3)	+ d*d*F_4_3		+ mywenoeps; \
const __m128 is2 = c*(c*F_10_3 - d*F_31_3 + e*F_11_3) + d*(d*F_25_3  - e*F_19_3)	+ e*e*F_4_3		+ mywenoeps; \
WENO_OMEGAS \
const __m128 recdata = ( \
omega0 *(F_1_3*a + M_7_6*b + F_11_6*c)+   \
omega1 *(M_1_6*b + F_5_6*c + F_1_3*d) +  \
omega2 *(F_1_3*c + F_5_6*d + M_1_6*e)	);

#define WENO_MINUS_DP \
const __m128d is0 = a*(a*F_4_3  - b*F_19_3 + c*F_11_3)	+ b*(b*F_25_3 - c*F_31_3)	+ c*c*F_10_3	+ mywenoeps; \
const __m128d is1 = b*(b*F_4_3  - c*F_13_3 + d*F_5_3)	+ c*(c*F_13_3 - d*F_13_3)	+ d*d*F_4_3		+ mywenoeps; \
const __m128d is2 = c*(c*F_10_3 - d*F_31_3 + e*F_11_3) + d*(d*F_25_3  - e*F_19_3)	+ e*e*F_4_3		+ mywenoeps; \
WENO_OMEGAS_DP \
const __m128d recdata = ( \
omega0 *(F_1_3*a + M_7_6*b + F_11_6*c)+   \
omega1 *(M_1_6*b + F_5_6*c + F_1_3*d) +  \
omega2 *(F_1_3*c + F_5_6*d + M_1_6*e)	);

#define WENO_PLUS \
const __m128 is0 = (d*(d*F_10_3 - e*F_31_3 + f*F_11_3)	+ e*(e*F_25_3  - f*F_19_3)	+ f*f*F_4_3	)	+ mywenoeps; \
const __m128 is1 = (c*(c*F_4_3  - d*F_13_3 + e*F_5_3 )	+ d*(d*F_13_3 - e*F_13_3)	+ e*e*F_4_3	)	+ mywenoeps; \
const __m128 is2 = (b*(b*F_4_3  - c*F_19_3 + d*F_11_3) + c*(c*F_25_3  - d*F_31_3)	+ d*d*F_10_3)	+ mywenoeps; \
WENO_OMEGAS \
const __m128 recdata = ( \
omega0 *(F_1_3*f + M_7_6*e + F_11_6*d)+   \
omega1 *(M_1_6*e + F_5_6*d + F_1_3*c) +  \
omega2 *(F_1_3*d + F_5_6*c + M_1_6*b)	);

#define WENO_PLUS_DP \
const __m128d is0 = (d*(d*F_10_3 - e*F_31_3 + f*F_11_3)	+ e*(e*F_25_3  - f*F_19_3)	+ f*f*F_4_3	)	+ mywenoeps; \
const __m128d is1 = (c*(c*F_4_3  - d*F_13_3 + e*F_5_3 )	+ d*(d*F_13_3 - e*F_13_3)	+ e*e*F_4_3	)	+ mywenoeps; \
const __m128d is2 = (b*(b*F_4_3  - c*F_19_3 + d*F_11_3) + c*(c*F_25_3  - d*F_31_3)	+ d*d*F_10_3)	+ mywenoeps; \
WENO_OMEGAS_DP \
const __m128d recdata = ( \
omega0 *(F_1_3*f + M_7_6*e + F_11_6*d)+   \
omega1 *(M_1_6*e + F_5_6*d + F_1_3*c) +  \
omega2 *(F_1_3*d + F_5_6*c + M_1_6*b)	);

#include "FlowStep_SSE_diego.h"

inline __m128 FlowStep_SSE_diego::_heaviside(const __m128 phi, const __m128 inv_h, const __m128 one, const __m128 phalf, const __m128 mhalf) const
{
	const __m128 x = _mm_min_ps(one, _mm_max_ps(_mm_setzero_ps() - one, phi*inv_h));
	
	const __m128 val_xneg = (mhalf*x - one)*x + phalf;
	const __m128 val_xpos = (phalf*x - one)*x + phalf;
	
	const __m128 flag = _mm_cmplt_ps(x, _mm_setzero_ps());
	
	return _mm_or_ps(_mm_and_ps(flag, val_xneg),_mm_andnot_ps(flag, val_xpos));
}

inline __m128d FlowStep_SSE_diego::_heaviside(const __m128d phi, const __m128d inv_h, const __m128d one, const __m128d phalf, const __m128d mhalf) const
{
	const __m128d x = _mm_min_pd(one, _mm_max_pd(_mm_setzero_pd() - one, phi*inv_h));
	
	const __m128d val_xneg = (mhalf*x - one)*x + phalf;
	const __m128d val_xpos = (phalf*x - one)*x + phalf;
	
	const __m128d flag = _mm_cmplt_pd(x, _mm_setzero_pd());
	
	return _mm_or_pd(_mm_and_pd(flag, val_xneg),_mm_andnot_pd(flag, val_xpos));
}

inline __m128 FlowStep_SSE_diego::_getgamma(const __m128 phi, const __m128 inv_smoothlength, 
											const __m128 gamma1, const __m128 gamma2,
											const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const
{
	const __m128 hs = _heaviside(phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	return (gamma1)*hs + (gamma2)*(F_1-hs); 
}

inline __m128d FlowStep_SSE_diego::_getgamma(const __m128d phi, const __m128d inv_smoothlength, 
											 const __m128d gamma1, const __m128d gamma2,
											 const __m128d F_1, const __m128d F_1_2, const __m128d M_1_2) const
{
	const __m128d hs = _heaviside(phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	return (gamma1)*hs + (gamma2)*(F_1-hs); 
}

inline __m128 FlowStep_SSE_diego::_getPC(const __m128 phi, const __m128 inv_smoothlength, 
                                         const __m128 pc1, const __m128 pc2,
                                         const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const
{
	const __m128 hs = _heaviside(phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	return (pc1)*hs + (pc2)*(F_1-hs); 
}

inline __m128d FlowStep_SSE_diego::_getPC(const __m128d phi, const __m128d inv_smoothlength, 
										  const __m128d pc1, const __m128d pc2,
										  const __m128d F_1, const __m128d F_1_2, const __m128d M_1_2) const
{
	const __m128d hs = _heaviside(phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	return (pc1)*hs + (pc2)*(F_1-hs); 
}

void FlowStep_SSE_diego::_sse_convert_aligned(const float * const gptfirst, const int gptfloats, const int rowgpts,
											  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l)
{
	const __m128 F_1_2 = _mm_set_ps1(0.5f);
	const __m128 M_1_2 = _mm_set_ps1(-0.5f);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/m_smoothlength);
	const __m128 g1 = _mm_set_ps1(m_gamma1);
	const __m128 g2 = _mm_set_ps1(m_gamma2);
    const __m128 pc1 = _mm_set_ps1(m_pc1);
	const __m128 pc2 = _mm_set_ps1(m_pc2);
	
#define DESTID (dx + (InputSOA::PITCH)*dy)
	
	for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
	{
		const float * const in = gptfirst + dy*gptfloats*rowgpts -gptfloats;
		
		for(int dx=0; dx<_BLOCKSIZE_+8; dx+=4)
		{
			const int WID = (dx + (int)(dx==0))*gptfloats;
			const int CID = (dx+1)*gptfloats;
			const int EID = (dx+3 - (int)(dx==_BLOCKSIZE_+4))*gptfloats;
			
			__m128 dataA0 = _mm_load_ps(in + WID);
			__m128 dataA1 = _mm_load_ps(in + WID + 4);
			__m128 dataB0 = _mm_load_ps(in + CID );
			__m128 dataB1 = _mm_load_ps(in + CID + 4);
			__m128 dataC0 = _mm_load_ps(in + CID + gptfloats);
			__m128 dataC1 = _mm_load_ps(in + CID + gptfloats + 4);
			__m128 dataD0 = _mm_load_ps(in + EID );
			__m128 dataD1 = _mm_load_ps(in + EID + 4);
			
			_MM_TRANSPOSE4_PS(dataA0, dataB0, dataC0, dataD0);
			
			_mm_store_ps(rho + DESTID, dataA0);
#ifdef _PREC_DIV_
			const __m128 inv_rho = F_1/dataA0;			
#else
			const __m128 inv_rho = better_rcp(dataA0);
#endif
			_mm_store_ps(u + DESTID, dataB0*inv_rho);
			_mm_store_ps(v + DESTID, dataC0*inv_rho);
			_mm_store_ps(w + DESTID, dataD0*inv_rho);
			
			_MM_TRANSPOSE4_PS(dataA1, dataB1, dataC1, dataD1);
			
			_mm_store_ps(l + DESTID, dataB1);
			
			_mm_store_ps(p + DESTID,  
						 (dataA1 - (dataB0*dataB0 + dataC0*dataC0 + dataD0*dataD0)*(F_1_2*inv_rho))*
						 (_getgamma(dataB1, invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1) -
						 _getgamma(dataB1, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(dataB1, invsml, pc1, pc2, F_1, F_1_2, M_1_2));
		}
	}
	
#undef DESTID
}

void FlowStep_SSE_diego::_sse_convert(const double * const gptfirst, const int gptfloats, const int rowgpts,
									  double * const rho, double * const u, double * const v, double * const w, double * const p, double * const l)
{
	const __m128d F_1_2 = _mm_set_pd1(0.5);
	const __m128d M_1_2 = _mm_set_pd1(-0.5);
	const __m128d F_1 = _mm_set_pd1(1);
	const __m128d invsml = _mm_set_pd1(1./m_smoothlength);
	const __m128d g1 = _mm_set_pd1(m_gamma1);
	const __m128d g2 = _mm_set_pd1(m_gamma2);
	const __m128d pc1 = _mm_set_pd1(m_pc1);
	const __m128d pc2 = _mm_set_pd1(m_pc2);
	
	assert(InputSOA::PITCH == _BLOCKSIZE_+8); 
	
#define DESTID (dx + (_BLOCKSIZE_+8)*dy)
	
	for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
	{
		const double * const in = gptfirst + dy*gptfloats*rowgpts - gptfloats;
		
		for(int dx=0; dx<_BLOCKSIZE_+8; dx+=2)
		{
			const int WID = (dx + (int)(dx==0))*gptfloats;
			const int EID = (dx + 1 - (int)(dx==_BLOCKSIZE_+6))*gptfloats;
			
			const __m128d dataA0 = _mm_load_pd(in + WID);
			const __m128d dataA1 = _mm_load_pd(in + WID + 2);
			const __m128d dataA2 = _mm_load_pd(in + WID + 4);
			
			const __m128d dataB0 = _mm_load_pd(in + EID );
			const __m128d dataB1 = _mm_load_pd(in + EID + 2);
			const __m128d dataB2 = _mm_load_pd(in + EID + 4);
			
			const __m128d myrho = _mm_shuffle_pd(dataA0, dataB0, _MM_SHUFFLE2(0,0));
			_mm_store_pd(rho + DESTID, myrho);
			
			const __m128d inv_rho = F_1/myrho;
			
			const __m128d myu = _mm_shuffle_pd(dataA0, dataB0, _MM_SHUFFLE2(1,1));
			const __m128d myv = _mm_shuffle_pd(dataA1, dataB1, _MM_SHUFFLE2(0,0));
			const __m128d myw = _mm_shuffle_pd(dataA1, dataB1, _MM_SHUFFLE2(1,1));
			const __m128d mye = _mm_shuffle_pd(dataA2, dataB2, _MM_SHUFFLE2(0,0));
			const __m128d myl = _mm_shuffle_pd(dataA2, dataB2, _MM_SHUFFLE2(1,1));
			
			_mm_store_pd(u + DESTID, myu*inv_rho);
			_mm_store_pd(v + DESTID, myv*inv_rho);
			_mm_store_pd(w + DESTID, myw*inv_rho);
			
			_mm_store_pd(p + DESTID,  (mye - (myu*myu + myv*myv + myw*myw)*(F_1_2*inv_rho))*
						 (_getgamma(myl, invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1) -
						 _getgamma(myl, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(myl, invsml, pc1, pc2, F_1, F_1_2, M_1_2));
			
			_mm_store_pd(l + DESTID, myl);
		}
	}
#undef DESTID
}

void FlowStep_SSE_diego::_sse_convert(const float * const gptfirst, const int gptfloats, const int rowgpts,
									  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l)
{
	const __m128 F_1_2 = _mm_set_ps1(0.5f);
	const __m128 M_1_2 = _mm_set_ps1(-0.5f);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/m_smoothlength);
	const __m128 g1 = _mm_set_ps1(m_gamma1);
	const __m128 g2 = _mm_set_ps1(m_gamma2);
	const __m128 pc1 = _mm_set_ps1(m_pc1);
	const __m128 pc2 = _mm_set_ps1(m_pc2);
	
#define DESTID (dx + (InputSOA::PITCH)*dy)
	
	for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
	{
		const float * const in = gptfirst + dy*gptfloats*rowgpts -gptfloats;
		
		for(int dx=0; dx<_BLOCKSIZE_+8; dx+=4)
		{
			const int WID = (dx + (int)(dx==0))*gptfloats;
			const int CID = (dx+1)*gptfloats;
			const int EID = (dx+3 - (int)(dx==_BLOCKSIZE_+4))*gptfloats;
			
			__m128 dataA0 = _mm_loadu_ps(in + WID);
			__m128 dataA1 = _mm_loadu_ps(in + WID + 4);
			__m128 dataB0 = _mm_loadu_ps(in + CID );
			__m128 dataB1 = _mm_loadu_ps(in + CID + 4);
			__m128 dataC0 = _mm_loadu_ps(in + CID + gptfloats);
			__m128 dataC1 = _mm_loadu_ps(in + CID + gptfloats + 4);
			__m128 dataD0 = _mm_loadu_ps(in + EID );
			__m128 dataD1 = _mm_loadu_ps(in + EID + 4);
			
			_MM_TRANSPOSE4_PS(dataA0, dataB0, dataC0, dataD0);
			
			_mm_store_ps(rho + DESTID, dataA0);
#ifdef _PREC_DIV_
			const __m128 inv_rho = F_1/dataA0;			
#else
			const __m128 inv_rho = better_rcp(dataA0);
#endif
			_mm_store_ps(u + DESTID, dataB0*inv_rho);
			_mm_store_ps(v + DESTID, dataC0*inv_rho);
			_mm_store_ps(w + DESTID, dataD0*inv_rho);
			
			_MM_TRANSPOSE4_PS(dataA1, dataB1, dataC1, dataD1);
			
			_mm_store_ps(l + DESTID, dataB1);
			
			_mm_store_ps(p + DESTID,  
						 (dataA1 - (dataB0*dataB0 + dataC0*dataC0 + dataD0*dataD0)*(F_1_2*inv_rho))*
						 (_getgamma(dataB1, invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1) -
						 _getgamma(dataB1, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(dataB1, invsml, pc1, pc2, F_1, F_1_2, M_1_2));
		}
	}
	
#undef DESTID
}

void FlowStep_SSE_diego::_sse_convert_aligned(const double * const gptfirst, const int gptfloats, const int rowgpts,
											  double * const rho, double * const u, double * const v, double * const w, double * const p, double * const l)
{
	const __m128d F_1_2 = _mm_set_pd1(0.5);
	const __m128d M_1_2 = _mm_set_pd1(-0.5);
	const __m128d F_1 = _mm_set_pd1(1);
	const __m128d invsml = _mm_set_pd1(1./m_smoothlength);
	const __m128d g1 = _mm_set_pd1(m_gamma1);
	const __m128d g2 = _mm_set_pd1(m_gamma2);
	const __m128d pc1 = _mm_set_pd1(m_pc1);
	const __m128d pc2 = _mm_set_pd1(m_pc2);
	
	assert(InputSOA::PITCH == _BLOCKSIZE_+8); 
	
#define DESTID (dx + (_BLOCKSIZE_+8)*dy)
	
	for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
	{
		const double * const in = gptfirst + dy*gptfloats*rowgpts - gptfloats;
		
		for(int dx=0; dx<_BLOCKSIZE_+8; dx+=2)
		{
			const int WID = (dx + (int)(dx==0))*gptfloats;
			const int EID = (dx + 1 - (int)(dx==_BLOCKSIZE_+6))*gptfloats;
			
			const __m128d dataA0 = _mm_loadu_pd(in + WID);
			const __m128d dataA1 = _mm_loadu_pd(in + WID + 2);
			const __m128d dataA2 = _mm_loadu_pd(in + WID + 4);
			
			const __m128d dataB0 = _mm_loadu_pd(in + EID );
			const __m128d dataB1 = _mm_loadu_pd(in + EID + 2);
			const __m128d dataB2 = _mm_loadu_pd(in + EID + 4);
			
			const __m128d myrho = _mm_shuffle_pd(dataA0, dataB0, _MM_SHUFFLE2(0,0));
			_mm_store_pd(rho + DESTID, myrho);
			
			const __m128d inv_rho = F_1/myrho;
			
			const __m128d myu = _mm_shuffle_pd(dataA0, dataB0, _MM_SHUFFLE2(1,1));
			const __m128d myv = _mm_shuffle_pd(dataA1, dataB1, _MM_SHUFFLE2(0,0));
			const __m128d myw = _mm_shuffle_pd(dataA1, dataB1, _MM_SHUFFLE2(1,1));
			const __m128d mye = _mm_shuffle_pd(dataA2, dataB2, _MM_SHUFFLE2(0,0));
			const __m128d myl = _mm_shuffle_pd(dataA2, dataB2, _MM_SHUFFLE2(1,1));
			
			_mm_store_pd(u + DESTID, myu*inv_rho);
			_mm_store_pd(v + DESTID, myv*inv_rho);
			_mm_store_pd(w + DESTID, myw*inv_rho);
			
			_mm_store_pd(p + DESTID,  (mye - (myu*myu + myv*myv + myw*myw)*(F_1_2*inv_rho))*
						 (_getgamma(myl, invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1) -
						 _getgamma(myl, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(myl, invsml, pc1, pc2, F_1, F_1_2, M_1_2));
			
			_mm_store_pd(l + DESTID, myl);
		}
	}
#undef DESTID
}

void FlowStep_SSE_diego::_sse_xweno_minus(const float * const in, float * const out) const
{
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=4)
		{			
			const __m128 W = _mm_load_ps(&in[dx + SX*dy]); 
			const __m128 C = _mm_load_ps(&in[dx+4 + SX*dy]);
			const __m128 E = _mm_load_ps(&in[dx+8 + SX*dy]);
			
			const __m128 a = _mm_shuffle_ps(W, _mm_shuffle_ps(W,C, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1));
			const __m128 b = _mm_shuffle_ps(W, C, _MM_SHUFFLE(1,0,3,2));
			const __m128 c = _mm_shuffle_ps(_mm_shuffle_ps(W,C, _MM_SHUFFLE(0,0,3,3)), C, _MM_SHUFFLE(2,1,3,0));;
			const __m128 d = C;
			const __m128 e = _mm_shuffle_ps(C, _mm_shuffle_ps(C,E, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1));
			
			WENO_MINUS
			
			_mm_store_ps(&out[dx + TempSOA::PITCH*dy], recdata);
		}
}

void FlowStep_SSE_diego::_sse_xweno_minus(const double * const in, double * const out) const
{
	WENO_CONSTANTS_DP
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=2)
		{			
			const __m128d Y = _mm_load_pd(&in[dx + SX*dy]); 
			const __m128d W = _mm_load_pd(&in[dx+2 + SX*dy]); 
			const __m128d C = _mm_load_pd(&in[dx+4 + SX*dy]);
			const __m128d E = _mm_load_pd(&in[dx+6 + SX*dy]);
			
			const __m128d a = _mm_shuffle_pd(Y, W, _MM_SHUFFLE2(0,1));
			const __m128d b = W;
			const __m128d c = _mm_shuffle_pd(W, C, _MM_SHUFFLE2(0,1));
			const __m128d d = C;
			const __m128d e = _mm_shuffle_pd(C,E, _MM_SHUFFLE2(0,1));
			
			WENO_MINUS_DP
			
			_mm_store_pd(&out[dx + TempSOA::PITCH*dy], recdata);
		}
}

void FlowStep_SSE_diego::_sse_xweno_pluss(const float * const in, float * const out) const
{
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=4)
		{			
			const __m128 W = _mm_load_ps(&in[dx + SX*dy]); 
			const __m128 C = _mm_load_ps(&in[dx+4 + SX*dy]);
			const __m128 E = _mm_load_ps(&in[dx+8 + SX*dy]);
			
			const __m128 b = _mm_shuffle_ps(W, C, _MM_SHUFFLE(1,0,3,2));
			const __m128 c = _mm_shuffle_ps(_mm_shuffle_ps(W,C, _MM_SHUFFLE(0,0,3,3)), C, _MM_SHUFFLE(2,1,3,0));
			const __m128 d = C;
			const __m128 e = _mm_shuffle_ps(C, _mm_shuffle_ps(C,E, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1));
			const __m128 f = _mm_shuffle_ps(C, E, _MM_SHUFFLE(1,0,3,2));
			
			WENO_PLUS
			
			_mm_store_ps(&out[dx + TempSOA::PITCH*dy], recdata);
		}
}

void FlowStep_SSE_diego::_sse_xweno_pluss(const double * const in, double * const out) const
{
	WENO_CONSTANTS_DP
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=2)
		{			
			const __m128d W = _mm_load_pd(&in[dx+2 + SX*dy]); 
			const __m128d C = _mm_load_pd(&in[dx+4 + SX*dy]);
			const __m128d E = _mm_load_pd(&in[dx+6 + SX*dy]);
			
			const __m128d b = W;
			const __m128d c = _mm_shuffle_pd(W, C, _MM_SHUFFLE2(0,1));
			const __m128d d = C;
			const __m128d e = _mm_shuffle_pd(C,E, _MM_SHUFFLE2(0,1));
			const __m128d f = E;
			
			WENO_PLUS_DP
			
			_mm_store_pd(&out[dx + TempSOA::PITCH*dy], recdata);
		}
}

void FlowStep_SSE_diego::_sse_yweno_minus(const float * const in, float * const out)
{	
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy+=4)
	{
		const float * ptr = &in[dy];
		
		for(int dx=0; dx<TempSOA::NX; dx++)
		{			
			const __m128 a = _mm_load_ps(ptr + dx*SX);
			const __m128 b = _mm_load_ps(ptr + dx*SX + SX);
			const __m128 c = _mm_load_ps(ptr + dx*SX + 2*SX);
			const __m128 d = _mm_load_ps(ptr + dx*SX + 3*SX);
			const __m128 e = _mm_load_ps(ptr + dx*SX + 4*SX);
			
			WENO_MINUS
			
			_mm_store_ps(&tmp[dx][0], recdata);
		}
		
		for(int dx=0; dx<TempSOA::NX-1; dx+=4)
		{
			__m128 data0 = _mm_load_ps(&tmp[dx][0]);
			__m128 data1 = _mm_load_ps(&tmp[dx+1][0]);
			__m128 data2 = _mm_load_ps(&tmp[dx+2][0]);
			__m128 data3 = _mm_load_ps(&tmp[dx+3][0]);
			
			_MM_TRANSPOSE4_PS(data0, data1, data2, data3);
			
			_mm_store_ps(&out[dx + TempSOA::PITCH*dy], data0);
			_mm_store_ps(&out[dx + TempSOA::PITCH*(dy+1)], data1);
			_mm_store_ps(&out[dx + TempSOA::PITCH*(dy+2)], data2);
			_mm_store_ps(&out[dx + TempSOA::PITCH*(dy+3)], data3);
		}
		
		{
			out[TempSOA::NX-1 + TempSOA::PITCH*dy] = tmp[TempSOA::NX-1][0];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+1)] = tmp[TempSOA::NX-1][1];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+2)] = tmp[TempSOA::NX-1][2];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+3)] = tmp[TempSOA::NX-1][3];
		}
	}
}

void FlowStep_SSE_diego::_sse_yweno_minus(const double * const in, double * const out)
{	
	double __attribute__((aligned(16))) tmp[TempSOA::NX][_CPER16BYTES];
	
	WENO_CONSTANTS_DP
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy+=2)
	{
		const double * ptr = &in[dy];
		
		for(int dx=0; dx<TempSOA::NX; dx++)
		{			
			const __m128d a = _mm_load_pd(ptr + dx*SX);
			const __m128d b = _mm_load_pd(ptr + dx*SX + SX);
			const __m128d c = _mm_load_pd(ptr + dx*SX + 2*SX);
			const __m128d d = _mm_load_pd(ptr + dx*SX + 3*SX);
			const __m128d e = _mm_load_pd(ptr + dx*SX + 4*SX);
			
			WENO_MINUS_DP
			
			_mm_store_pd(&tmp[dx][0], recdata);
		}
		
		for(int dx=0; dx<TempSOA::NX-1; dx+=2)
		{
			const __m128d data0 = _mm_load_pd(&tmp[dx][0]);
			const __m128d data1 = _mm_load_pd(&tmp[dx+1][0]);
			
			_mm_store_pd(&out[dx + TempSOA::PITCH*dy], _mm_shuffle_pd(data0, data1, _MM_SHUFFLE2(0,0)));
			_mm_store_pd(&out[dx + TempSOA::PITCH*(dy+1)], _mm_shuffle_pd(data0, data1, _MM_SHUFFLE2(1,1)));
		}
		
		{
			out[TempSOA::NX-1 + TempSOA::PITCH*dy] = tmp[TempSOA::NX-1][0];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+1)] = tmp[TempSOA::NX-1][1];
		}
	}
}

void FlowStep_SSE_diego::_sse_yweno_pluss(const float * const in, float * const out)
{	
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy+=4)
	{
		const float * ptr = &in[dy];
		
		for(int dx=0; dx<TempSOA::NX; dx++)
		{			
			const __m128 b = _mm_load_ps(ptr + dx*SX + SX);
			const __m128 c = _mm_load_ps(ptr + dx*SX + 2*SX);
			const __m128 d = _mm_load_ps(ptr + dx*SX + 3*SX);
			const __m128 e = _mm_load_ps(ptr + dx*SX + 4*SX);
			const __m128 f = _mm_load_ps(ptr + dx*SX + 5*SX);
			
			WENO_PLUS
			
			_mm_store_ps(&tmp[dx][0], recdata);
		}
		
		for(int dx=0; dx<TempSOA::NX-1; dx+=4)
		{
			__m128 data0 = _mm_load_ps(&tmp[dx][0]);
			__m128 data1 = _mm_load_ps(&tmp[dx+1][0]);
			__m128 data2 = _mm_load_ps(&tmp[dx+2][0]);
			__m128 data3 = _mm_load_ps(&tmp[dx+3][0]);
			
			_MM_TRANSPOSE4_PS(data0, data1, data2, data3);
			
			_mm_store_ps(&out[dx + TempSOA::PITCH*dy], data0);
			_mm_store_ps(&out[dx + TempSOA::PITCH*(dy+1)], data1);
			_mm_store_ps(&out[dx + TempSOA::PITCH*(dy+2)], data2);
			_mm_store_ps(&out[dx + TempSOA::PITCH*(dy+3)], data3);
		}
		
		{
			out[TempSOA::NX-1 + TempSOA::PITCH*dy] = tmp[TempSOA::NX-1][0];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+1)] = tmp[TempSOA::NX-1][1];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+2)] = tmp[TempSOA::NX-1][2];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+3)] = tmp[TempSOA::NX-1][3];
		}
		
	}
}

void FlowStep_SSE_diego::_sse_yweno_pluss(const double * const in, double * const out)
{	
	double __attribute__((aligned(16))) tmp[TempSOA::NX][_CPER16BYTES];
	
	WENO_CONSTANTS_DP
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy+=2)
	{
		const double * ptr = &in[dy];
		
		for(int dx=0; dx<TempSOA::NX; dx++)
		{			
			const __m128d b = _mm_load_pd(ptr + dx*SX + SX);
			const __m128d c = _mm_load_pd(ptr + dx*SX + 2*SX);
			const __m128d d = _mm_load_pd(ptr + dx*SX + 3*SX);
			const __m128d e = _mm_load_pd(ptr + dx*SX + 4*SX);
			const __m128d f = _mm_load_pd(ptr + dx*SX + 5*SX);
			
			WENO_PLUS_DP
			
			_mm_store_pd(&tmp[dx][0], recdata);
		}
		
		for(int dx=0; dx<TempSOA::NX-1; dx+=2)
		{
			const __m128d data0 = _mm_load_pd(&tmp[dx][0]);
			const __m128d data1 = _mm_load_pd(&tmp[dx+1][0]);
			
			_mm_store_pd(&out[dx + TempSOA::PITCH*dy], _mm_shuffle_pd(data0, data1, _MM_SHUFFLE2(0,0)));
			_mm_store_pd(&out[dx + TempSOA::PITCH*(dy+1)], _mm_shuffle_pd(data0, data1, _MM_SHUFFLE2(1,1)));
		}
		
		{
			out[TempSOA::NX-1 + TempSOA::PITCH*dy] = tmp[TempSOA::NX-1][0];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+1)] = tmp[TempSOA::NX-1][1];
		}
	}
}

void FlowStep_SSE_diego::_sse_zweno_minus(const float * const a_, const float * const b_, 
										  const float * const c_, const float * const d_, 
										  const float * const e_ , float * const out) const
{
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=4)
		{			
			const __m128 a = _mm_load_ps(a_ + dx + SX*dy);
			const __m128 b = _mm_load_ps(b_ + dx + SX*dy);
			const __m128 c = _mm_load_ps(c_ + dx + SX*dy);
			const __m128 d = _mm_load_ps(d_ + dx + SX*dy);
			const __m128 e = _mm_load_ps(e_ + dx + SX*dy);
			
			WENO_MINUS
			
			_mm_store_ps(&out[dx + TempSOA::PITCH*dy], recdata);
		}	
}

void FlowStep_SSE_diego::_sse_zweno_minus(const double * const a_, const double * const b_, 
										  const double * const c_, const double * const d_, 
										  const double * const e_ , double * const out) const
{
	WENO_CONSTANTS_DP
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=2)
		{			
			const __m128d a = _mm_load_pd(a_ + dx + SX*dy);
			const __m128d b = _mm_load_pd(b_ + dx + SX*dy);
			const __m128d c = _mm_load_pd(c_ + dx + SX*dy);
			const __m128d d = _mm_load_pd(d_ + dx + SX*dy);
			const __m128d e = _mm_load_pd(e_ + dx + SX*dy);
			
			WENO_MINUS_DP
			
			_mm_store_pd(&out[dx + TempSOA::PITCH*dy], recdata);
		}	
}

void FlowStep_SSE_diego::_sse_zweno_pluss(const float * const b_, const float * const c_, 
										  const float * const d_, const float * const e_, 
										  const float * const f_ , float * const out) const
{
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=4)
		{			
			const __m128 b = _mm_load_ps(b_ + dx + SX*dy);
			const __m128 c = _mm_load_ps(c_ + dx + SX*dy);
			const __m128 d = _mm_load_ps(d_ + dx + SX*dy);
			const __m128 e = _mm_load_ps(e_ + dx + SX*dy);
			const __m128 f = _mm_load_ps(f_ + dx + SX*dy);
			
			WENO_PLUS
			
			_mm_store_ps(&out[dx + TempSOA::PITCH*dy], recdata);
		}	
}

void FlowStep_SSE_diego::_sse_zweno_pluss(const double * const b_, const double * const c_, 
										  const double * const d_, const double * const e_, 
										  const double * const f_ , double * const out) const
{
	WENO_CONSTANTS_DP
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=2)
		{			
			const __m128d b = _mm_load_pd(b_ + dx + SX*dy);
			const __m128d c = _mm_load_pd(c_ + dx + SX*dy);
			const __m128d d = _mm_load_pd(d_ + dx + SX*dy);
			const __m128d e = _mm_load_pd(e_ + dx + SX*dy);
			const __m128d f = _mm_load_pd(f_ + dx + SX*dy);
			
			WENO_PLUS_DP
			
			_mm_store_pd(&out[dx + TempSOA::PITCH*dy], recdata);
		}	
}

void FlowStep_SSE_diego::_sse_hlle_rho(const float * const rm, const float * const rp,
									   const float * const vm, const float * const vp,
									   const float * const am, const float * const ap,
									   float * const out)
{
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)	
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=4)
		{
			const __m128 flagminus	= _mm_cmpgt_ps(_mm_load_ps(am + ID),_mm_setzero_ps());
			const __m128 nonflagplus  = _mm_cmpge_ps(_mm_load_ps(ap + ID),_mm_setzero_ps()); 
			const __m128 flagother	= _mm_andnot_ps(flagminus, nonflagplus);
			
			const __m128 rminus = _mm_load_ps(rm + ID);
			const __m128 rpluss = _mm_load_ps(rp + ID);
			
			const __m128 fminus = _mm_load_ps(vm + ID)*rminus;
			const __m128 fpluss = _mm_load_ps(vp + ID)*rpluss;
			
#define aminus _mm_load_ps(am + ID)
#define apluss _mm_load_ps(ap + ID)
			
#ifdef _PREC_DIV_						
			const __m128 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(rpluss-rminus))/(apluss-aminus);
#else
			const __m128 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(rpluss-rminus))*better_rcp(apluss-aminus);
#endif
			_mm_store_ps(out + ID, _3ORPS_(_mm_and_ps(flagminus, fminus), _mm_andnot_ps(nonflagplus, fpluss), _mm_and_ps(flagother, fother)));
		}
#undef ID
#undef aminus
#undef apluss
}

void FlowStep_SSE_diego::_sse_hlle_rho(const double * const rm, const double * const rp,
									   const double * const vm, const double * const vp,
									   const double * const am, const double * const ap,
									   double * const out)
{
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)	
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=2)
		{
			const __m128d flagminus	= _mm_cmpgt_pd(_mm_load_pd(am + ID),_mm_setzero_pd());
			const __m128d nonflagplus  = _mm_cmpge_pd(_mm_load_pd(ap + ID),_mm_setzero_pd()); 
			const __m128d flagother	= _mm_andnot_pd(flagminus, nonflagplus);
			
			const __m128d rminus = _mm_load_pd(rm + ID);
			const __m128d rpluss = _mm_load_pd(rp + ID);
			
			const __m128d fminus = _mm_load_pd(vm + ID)*rminus;
			const __m128d fpluss = _mm_load_pd(vp + ID)*rpluss;
			
			const __m128d aminus = _mm_load_pd(am + ID);
			const __m128d apluss  = _mm_load_pd(ap + ID);
			
			const __m128d fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(rpluss-rminus))/(apluss-aminus);
			
			_mm_store_pd(out + ID, _3ORPD_(_mm_and_pd(flagminus , fminus), _mm_andnot_pd(nonflagplus, fpluss), _mm_and_pd(flagother, fother)));
		}
#undef ID
}

void FlowStep_SSE_diego::_sse_hlle_vel(const float * const rm, const float * const rp,
									   const float * const vm, const float * const vp,
									   const float * const vdm, const float * const vdp,
									   const float * const am, const float * const ap,
									   float * const out)
{
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)	
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=4)
		{
			const __m128 flagminus	= _mm_cmpgt_ps(_mm_load_ps(am + ID),_mm_setzero_ps());
			const __m128 nonflagplus  = _mm_cmpge_ps(_mm_load_ps(ap + ID),_mm_setzero_ps()); 
			const __m128 flagother	= _mm_andnot_ps(flagminus, nonflagplus);
			
			const __m128 uminus = _mm_load_ps(vm + ID)*_mm_load_ps(rm + ID);
			const __m128 upluss = _mm_load_ps(vp + ID)*_mm_load_ps(rp + ID);
			
			const __m128 fminus = _mm_load_ps(vdm + ID)*uminus;
			const __m128 fpluss = _mm_load_ps(vdp + ID)*upluss;
			
#define aminus _mm_load_ps(am + ID)
#define apluss _mm_load_ps(ap + ID)
			
#ifdef _PREC_DIV_						
			const __m128 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))/(apluss-aminus);
#else
			const __m128 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))*better_rcp(apluss-aminus);
#endif
			_mm_store_ps(out + ID, _3ORPS_(_mm_and_ps(flagminus, fminus),_mm_andnot_ps(nonflagplus, fpluss), _mm_and_ps(flagother, fother)));
		}
#undef ID
#undef aminus
#undef apluss
}

void FlowStep_SSE_diego::_sse_hlle_vel(const double * const rm, const double * const rp,
									   const double * const vm, const double * const vp,
									   const double * const vdm, const double * const vdp,
									   const double * const am, const double * const ap,
									   double * const out)
{
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)	
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=2)
		{
			const __m128d flagminus	= _mm_cmpgt_pd(_mm_load_pd(am + ID),_mm_setzero_pd());
			const __m128d nonflagplus  = _mm_cmpge_pd(_mm_load_pd(ap + ID),_mm_setzero_pd()); 
			const __m128d flagother	= _mm_andnot_pd(flagminus, nonflagplus);
			
			const __m128d uminus = _mm_load_pd(vm + ID)*_mm_load_pd(rm + ID);
			const __m128d upluss = _mm_load_pd(vp + ID)*_mm_load_pd(rp + ID);
			
			const __m128d fminus = _mm_load_pd(vdm + ID)*uminus;
			const __m128d fpluss = _mm_load_pd(vdp + ID)*upluss;
			
			const __m128d aminus = _mm_load_pd(am + ID);
			const __m128d apluss  = _mm_load_pd(ap + ID);
			
			const __m128d fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))/(apluss-aminus);
			
			_mm_store_pd(out + ID, _3ORPD_(_mm_and_pd(flagminus, fminus), _mm_andnot_pd(nonflagplus, fpluss), _mm_and_pd(flagother, fother)));
		}
#undef ID
}

void FlowStep_SSE_diego::_sse_hlle_pvel(const float * const rm, const float * const rp,
										const float * const vm, const float * const vp,
										const float * const pm, const float * const pp,
										const float * const am, const float * const ap,
										float * const out)
{
	static const int P = TempSOA::PITCH;
#define ID (ix + P*iy)
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=4)
		{			
#define myvminus _mm_load_ps(vm + ID)
#define myvpluss _mm_load_ps(vp + ID)
			
			const __m128 uminus = myvminus*_mm_load_ps(rm + ID);
			const __m128 upluss = myvpluss*_mm_load_ps(rp + ID);
			
			const __m128 fminus = myvminus*uminus + _mm_load_ps(pm + ID);
			const __m128 fpluss = myvpluss*upluss + _mm_load_ps(pp + ID);
			
#define aminus _mm_load_ps(am + ID)
#define apluss _mm_load_ps(ap + ID)
			
#define flagminus	_mm_cmpgt_ps(aminus,_mm_setzero_ps())
#define nonflagplus _mm_cmpge_ps(apluss,_mm_setzero_ps()) 
#define flagother	_mm_andnot_ps(flagminus, nonflagplus)		
			
#ifdef _PREC_DIV_						
#define fother (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))/(apluss-aminus)
#else
#define fother (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))*better_rcp(apluss-aminus)
#endif
			_mm_store_ps(out + ID, _3ORPS_(_mm_and_ps(flagminus, fminus), _mm_andnot_ps(nonflagplus, fpluss), _mm_and_ps(flagother, fother)));
		}
#undef myvminus
#undef myvpluss
#undef ID
#undef aminus
#undef apluss
#undef flagminus
#undef nonflagplus
#undef flagother
#undef fother
}

void FlowStep_SSE_diego::_sse_hlle_pvel(const double * const rm, const double * const rp,
										const double * const vm, const double * const vp,
										const double * const pm, const double * const pp,
										const double * const am, const double * const ap,
										double * const out)
{
	static const int P = TempSOA::PITCH;
#define ID (ix + P*iy)
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=2)
		{
			const __m128d flagminus	= _mm_cmpgt_pd(_mm_load_pd(am + ID),_mm_setzero_pd());
			const __m128d nonflagplus  = _mm_cmpge_pd(_mm_load_pd(ap + ID),_mm_setzero_pd()); 
			const __m128d flagother	= _mm_andnot_pd(flagminus, nonflagplus);
			
			const __m128d myvminus = _mm_load_pd(vm + ID);
			const __m128d myvpluss = _mm_load_pd(vp + ID);
			
			const __m128d uminus = myvminus*_mm_load_pd(rm + ID);
			const __m128d upluss = myvpluss*_mm_load_pd(rp + ID);
			
			const __m128d fminus = myvminus*uminus + _mm_load_pd(pm + ID);
			const __m128d fpluss = myvpluss*upluss + _mm_load_pd(pp + ID);
			
			const __m128d aminus = _mm_load_pd(am + ID);
			const __m128d apluss  = _mm_load_pd(ap + ID);
			
			const __m128d fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))/(apluss-aminus);
			
			_mm_store_pd(out + ID, _3ORPD_(_mm_and_pd(flagminus, fminus), _mm_andnot_pd(nonflagplus, fpluss), _mm_and_pd(flagother, fother)));
		}
#undef ID
}

void FlowStep_SSE_diego::_sse_hlle_e(const float * const rm, const float * const rp,
									 const float * const vdm, const float * const vdp,
									 const float * const v1m, const float * const v1p,
									 const float * const v2m, const float * const v2p,
									 const float * const pm, const float * const pp,
									 const float * const lm, const float * const lp, 
									 const float * const am, const float * const ap,
									 float * const out)
{
	static const int P = TempSOA::PITCH;
#define ID (ix + P*iy)
	
	const __m128 F_1_2 = _mm_set_ps1(0.5);
	const __m128 M_1_2 = _mm_set_ps1(-0.5);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/m_smoothlength);
	const __m128 g1 = _mm_set_ps1(m_gamma1);
	const __m128 g2 = _mm_set_ps1(m_gamma2);
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=4)
		{			
#define vdminus  _mm_load_ps(vdm + ID)
#define v1minus  _mm_load_ps(v1m + ID)
#define v2minus  _mm_load_ps(v2m + ID)
#define pminus  _mm_load_ps(pm + ID)
#ifdef _PREC_DIV_
			const __m128 eminus = pminus*(F_1/(_getgamma(_mm_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm_load_ps(rm + ID)*(vdminus*vdminus + v1minus*v1minus + v2minus*v2minus);
#else
			const __m128 eminus = pminus*(F_1 * better_rcp(_getgamma(_mm_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm_load_ps(rm + ID)*(vdminus*vdminus + v1minus*v1minus + v2minus*v2minus);
#endif
#define  vdplus  _mm_load_ps(vdp + ID)
#define  v1plus  _mm_load_ps(v1p + ID)
#define v2plus  _mm_load_ps(v2p + ID)
#define pplus  _mm_load_ps(pp + ID)
#ifdef _PREC_DIV_			
			const __m128 eplus = pplus*(F_1/(_getgamma(_mm_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm_load_ps(rp + ID)*(vdplus*vdplus + v1plus*v1plus + v2plus*v2plus);
#else
			const __m128 eplus = pplus*(F_1 * better_rcp(_getgamma(_mm_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm_load_ps(rp + ID)*(vdplus*vdplus + v1plus*v1plus + v2plus*v2plus);
#endif
			const __m128 fminus = vdminus*(pminus + eminus);
			const __m128 fpluss = vdplus *(pplus + eplus);
			
#define aminus _mm_load_ps(am + ID)
#define apluss _mm_load_ps(ap + ID)
			
#define flagminus	_mm_cmpgt_ps(aminus,_mm_setzero_ps())
#define nonflagplus _mm_cmpge_ps(apluss,_mm_setzero_ps()) 
#define flagother	_mm_andnot_ps(flagminus, nonflagplus)
			
#ifdef _PREC_DIV_						
#define fother	(apluss*fminus-aminus*fpluss+aminus*apluss*(eplus-eminus))/(apluss-aminus)
#else
#define fother	(apluss*fminus-aminus*fpluss+aminus*apluss*(eplus-eminus))*better_rcp(apluss-aminus)
#endif		
			_mm_store_ps(out + ID, _3ORPS_(_mm_and_ps(flagminus, fminus), _mm_andnot_ps(nonflagplus, fpluss), _mm_and_ps(flagother, fother)));
		}
#undef ID
#undef vdminus
#undef v1minus
#undef v2minus
#undef pminus
#undef vdplus
#undef v1plus
#undef v2plus
#undef pplus
#undef aminus
#undef apluss
#undef flagminus
#undef nonflagplus
#undef flagother
#undef fother
}

void FlowStep_SSE_diego::_sse_hlle_e(const double * const rm, const double * const rp,
									 const double * const vdm, const double * const vdp,
									 const double * const v1m, const double * const v1p,
									 const double * const v2m, const double * const v2p,
									 const double * const pm, const double * const pp,
									 const double * const lm, const double * const lp, 
									 const double * const am, const double * const ap,
									 double * const out)
{
	static const int P = TempSOA::PITCH;
#define ID (ix + P*iy)
	
	const __m128d F_1_2 = _mm_set_pd1(0.5);
	const __m128d M_1_2 = _mm_set_pd1(-0.5);
	const __m128d F_1 = _mm_set_pd1(1);
	const __m128d invsml = _mm_set_pd1(1/m_smoothlength);
	const __m128d g1 = _mm_set_pd1(m_gamma1);
	const __m128d g2 = _mm_set_pd1(m_gamma2);
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=2)
		{
			const __m128d flagminus	= _mm_cmpgt_pd(_mm_load_pd(am + ID),_mm_setzero_pd());
			const __m128d nonflagplus  = _mm_cmpge_pd(_mm_load_pd(ap + ID),_mm_setzero_pd()); 
			const __m128d flagother	= _mm_andnot_pd(flagminus, nonflagplus);
			
			const __m128d vdminus = _mm_load_pd(vdm + ID);
			const __m128d v1minus = _mm_load_pd(v1m + ID);
			const __m128d v2minus = _mm_load_pd(v2m + ID);
			const __m128d pminus = _mm_load_pd(pm + ID);
			
			const __m128d eminus = pminus*(F_1/(_getgamma(_mm_load_pd(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm_load_pd(rm + ID)*(vdminus*vdminus + v1minus*v1minus + v2minus*v2minus);
			
			const __m128d vdplus = _mm_load_pd(vdp + ID);
			const __m128d v1plus = _mm_load_pd(v1p + ID);
			const __m128d v2plus = _mm_load_pd(v2p + ID);
			const __m128d pplus = _mm_load_pd(pp + ID);
			
			const __m128d eplus = pplus*(F_1/(_getgamma(_mm_load_pd(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm_load_pd(rp + ID)*(vdplus*vdplus + v1plus*v1plus + v2plus*v2plus);
			
			const __m128d fminus = vdminus*(pminus + eminus);
			const __m128d fpluss = vdplus *(pplus + eplus);
			
			const __m128d aminus = _mm_load_pd(am + ID);
			const __m128d aplus  = _mm_load_pd(ap + ID);
			
			const __m128d fother = (aplus*fminus-aminus*fpluss+aminus*aplus*(eplus-eminus))/(aplus-aminus);
			
			_mm_store_pd(out + ID, _3ORPD_(_mm_and_pd(flagminus, fminus), _mm_andnot_pd(nonflagplus, fpluss), _mm_and_pd(flagother, fother)));
		}
#undef ID
}

void FlowStep_SSE_diego::_sse_char_vel(const float * const rm, const float * const rp, 
									   const float * const vm, const float * const vp,
									   const float * const pm, const float * const pp,
									   const float * const lm, const float * const lp, 
									   float * const outm, float * const outp)
{
	const __m128 F_1_2 = _mm_set_ps1(0.5);
	const __m128 M_1_2 = _mm_set_ps1(-0.5);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/m_smoothlength);
	const __m128 g1 = _mm_set_ps1(m_gamma1);
	const __m128 g2 = _mm_set_ps1(m_gamma2);
	const __m128 pc1 = _mm_set_ps1(m_pc1);
	const __m128 pc2 = _mm_set_ps1(m_pc2);
	
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=4)
		{
#ifdef _PREC_DIV_
			const __m128 cminus = _mm_sqrt_ps(_getgamma(_mm_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											  _mm_max_ps((_mm_load_ps(pm + ID)+_getPC(_mm_load_ps(lm + ID), invsml, pc1, pc2, F_1, F_1_2, M_1_2))/_mm_load_ps(rm + ID), _mm_setzero_ps()));
			const __m128 cplus = _mm_sqrt_ps(_getgamma(_mm_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											 _mm_max_ps((_mm_load_ps(pp + ID)+_getPC(_mm_load_ps(lp + ID), invsml, pc1, pc2, F_1, F_1_2, M_1_2))/_mm_load_ps(rp + ID), _mm_setzero_ps()));
#else
			const __m128 cminus = worse_sqrt(_getgamma(_mm_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											 _mm_max_ps((_mm_load_ps(pm + ID)+_getPC(_mm_load_ps(lm + ID), invsml, pc1, pc2, F_1, F_1_2, M_1_2))* better_rcp(_mm_load_ps(rm + ID)), _mm_setzero_ps()));
			const __m128 cplus = worse_sqrt(_getgamma(_mm_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											_mm_max_ps((_mm_load_ps(pp + ID)+_getPC(_mm_load_ps(lp + ID), invsml, pc1, pc2, F_1, F_1_2, M_1_2))* better_rcp(_mm_load_ps(rp + ID)), _mm_setzero_ps()));
#endif
			_mm_store_ps(outm + ID, _mm_min_ps(_mm_load_ps(vm + ID) - cminus, _mm_load_ps(vm + ID) - cplus));
			_mm_store_ps(outp + ID, _mm_max_ps(_mm_load_ps(vp + ID) + cminus, _mm_load_ps(vp + ID) + cplus));
		}
#undef ID
}

void FlowStep_SSE_diego::_sse_char_vel(const double * const rm, const double * const rp, 
									   const double * const vm, const double * const vp,
									   const double * const pm, const double * const pp,
									   const double * const lm, const double * const lp, 
									   double * const outm, double * const outp)
{
	const __m128d F_1_2 = _mm_set_pd1(0.5);
	const __m128d M_1_2 = _mm_set_pd1(-0.5);
	const __m128d F_1 = _mm_set_pd1(1);
	const __m128d invsml = _mm_set_pd1(1/m_smoothlength);
	const __m128d g1 = _mm_set_pd1(m_gamma1);
	const __m128d g2 = _mm_set_pd1(m_gamma2);
	const __m128d pc1 = _mm_set_pd1(m_pc1);
	const __m128d pc2 = _mm_set_pd1(m_pc2);
	
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=2)
		{
			const __m128d cminus = _mm_sqrt_pd(_getgamma(_mm_load_pd(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											   _mm_max_pd((_mm_load_pd(pm + ID)+_getPC(_mm_load_pd(lm + ID), invsml, pc1, pc2, F_1, F_1_2, M_1_2))/_mm_load_pd(rm + ID), _mm_setzero_pd()));
			const __m128d cplus = _mm_sqrt_pd(_getgamma(_mm_load_pd(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											  _mm_max_pd((_mm_load_pd(pp + ID)+_getPC(_mm_load_pd(lp + ID), invsml, pc1, pc2, F_1, F_1_2, M_1_2))/_mm_load_pd(rp + ID), _mm_setzero_pd()));
			_mm_store_pd(outm + ID, _mm_min_pd(_mm_load_pd(vm + ID) - cminus, _mm_load_pd(vm + ID) - cplus));
			_mm_store_pd(outp + ID, _mm_max_pd(_mm_load_pd(vp + ID) + cminus, _mm_load_pd(vp + ID) + cplus));
		}
#undef ID
}

void FlowStep_SSE_diego::_xflux(const int relid)
{	
	_xweno_minus(ringrho(relid), wenorho.ref(0));
	_xweno_pluss(ringrho(relid), wenorho.ref(1));
	_xweno_minus(ringu(relid), wenou.ref(0));
	_xweno_pluss(ringu(relid), wenou.ref(1));
	_xweno_minus(ringv(relid), wenov.ref(0));
	_xweno_pluss(ringv(relid), wenov.ref(1));
	_xweno_minus(ringw(relid), wenow.ref(0));
	_xweno_pluss(ringw(relid), wenow.ref(1));
	_xweno_minus(ringp(relid), wenop.ref(0));
	_xweno_pluss(ringp(relid), wenop.ref(1));
	_xweno_minus(ringls(relid), wenols.ref(0));
	_xweno_pluss(ringls(relid), wenols.ref(1));
	
	_xextraterm(wenou(0), wenou(1), wenols(0), wenols(1));	
	_char_vel(wenorho(0), wenorho(1), wenou(0), wenou(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel.ref(0), charvel.ref(1));
	
	_hlle_rho(wenorho(0), wenorho(1), wenou(0), wenou(1), charvel(0), charvel(1), fluxrho.ref());
	_hlle_pvel(wenorho(0), wenorho(1), wenou(0), wenou(1), wenop(0), wenop(1), charvel(0), charvel(1), fluxu.ref());
	_hlle_vel(wenorho(0), wenorho(1), wenov(0), wenov(1), wenou(0), wenou(1), charvel(0), charvel(1), fluxv.ref());
	_hlle_vel(wenorho(0), wenorho(1), wenow(0), wenow(1), wenou(0), wenou(1), charvel(0), charvel(1), fluxw.ref());
	_hlle_e(wenorho(0), wenorho(1), wenou(0), wenou(1), wenov(0), wenov(1), wenow(0), wenow(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel(0), charvel(1), fluxp.ref());
	_hlle_rho(wenols(0), wenols(1), wenou(0), wenou(1), charvel(0), charvel(1), fluxls.ref());
}

void FlowStep_SSE_diego::_yflux(const int relid)
{	
	_yweno_minus(ringrho(relid), wenorho.ref(0));
	_yweno_pluss(ringrho(relid), wenorho.ref(1));
	_yweno_minus(ringu(relid), wenou.ref(0));
	_yweno_pluss(ringu(relid), wenou.ref(1));
	_yweno_minus(ringv(relid), wenov.ref(0));
	_yweno_pluss(ringv(relid), wenov.ref(1));
	_yweno_minus(ringw(relid), wenow.ref(0));
	_yweno_pluss(ringw(relid), wenow.ref(1));
	_yweno_minus(ringp(relid), wenop.ref(0));
	_yweno_pluss(ringp(relid), wenop.ref(1));
	_yweno_minus(ringls(relid), wenols.ref(0));
	_yweno_pluss(ringls(relid), wenols.ref(1));
	
	_yextraterm(wenov(0), wenov(1), wenols(0), wenols(1));
	_char_vel(wenorho(0), wenorho(1), wenov(0), wenov(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel.ref(0), charvel.ref(1));
	
	_hlle_rho(wenorho(0), wenorho(1), wenov(0), wenov(1), charvel(0), charvel(1), fluxrho.ref());
	_hlle_vel(wenorho(0), wenorho(1), wenou(0), wenou(1), wenov(0), wenov(1), charvel(0), charvel(1), fluxu.ref());
	_hlle_pvel(wenorho(0), wenorho(1), wenov(0), wenov(1), wenop(0), wenop(1), charvel(0), charvel(1), fluxv.ref());
	_hlle_vel(wenorho(0), wenorho(1), wenow(0), wenow(1), wenov(0), wenov(1), charvel(0), charvel(1), fluxw.ref());
	_hlle_e(wenorho(0), wenorho(1), wenov(0), wenov(1), wenou(0), wenou(1), wenow(0), wenow(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel(0), charvel(1), fluxp.ref());
	_hlle_rho(wenols(0), wenols(1), wenov(0), wenov(1), charvel(0), charvel(1), fluxls.ref());
}

void FlowStep_SSE_diego::_zflux(const int relid)
{	
	_zweno_minus(relid, ringrho, wenorho.ref(0));
	_zweno_pluss(relid, ringrho, wenorho.ref(1));
	_zweno_minus(relid, ringu, wenou.ref(0));
	_zweno_pluss(relid, ringu, wenou.ref(1));
	_zweno_minus(relid, ringv, wenov.ref(0));
	_zweno_pluss(relid, ringv, wenov.ref(1));
	_zweno_minus(relid, ringw, wenow.ref(0));
	_zweno_pluss(relid, ringw, wenow.ref(1));
	_zweno_minus(relid, ringp, wenop.ref(0));
	_zweno_pluss(relid, ringp, wenop.ref(1));
	_zweno_minus(relid, ringls, wenols.ref(0));
	_zweno_pluss(relid, ringls, wenols.ref(1));
	
	_zextraterm(wenow(0), wenow(-1), wenols(0), wenols(-1));
	_char_vel(wenorho(0), wenorho(1), wenow(0), wenow(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel.ref(0), charvel.ref(1));
	
	_hlle_rho(wenorho(0), wenorho(1), wenow(0), wenow(1), charvel(0), charvel(1), fluxrho.ref());
	_hlle_vel(wenorho(0), wenorho(1), wenou(0), wenou(1), wenow(0), wenow(1), charvel(0), charvel(1), fluxu.ref());
	_hlle_vel(wenorho(0), wenorho(1), wenov(0), wenov(1), wenow(0), wenow(1), charvel(0), charvel(1), fluxv.ref());
	_hlle_pvel(wenorho(0), wenorho(1), wenow(0), wenow(1), wenop(0), wenop(1), charvel(0), charvel(1), fluxw.ref());
	_hlle_e(wenorho(0), wenorho(1), wenow(0), wenow(1), wenou(0), wenou(1), wenov(0), wenov(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel(0), charvel(1), fluxp.ref());
	_hlle_rho(wenols(0), wenols(1), wenow(0), wenow(1), charvel(0), charvel(1), fluxls.ref());
}

void FlowStep_SSE_diego::_sse_xrhsadd(const float * const f, float * const r)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
			_mm_store_ps(r + ix + OutputSOA::PITCH*iy,  
						 _mm_loadu_ps(f + ix + 1 + TempSOA::PITCH*iy) - _mm_load_ps(f + ix + TempSOA::PITCH*iy));
}

void FlowStep_SSE_diego::_sse_xrhsadd(const double * const f, double * const r)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=2)
			_mm_store_pd(r + ix + OutputSOA::PITCH*iy,  
						 _mm_loadu_pd(f + ix + 1 + TempSOA::PITCH*iy) - _mm_load_pd(f + ix + TempSOA::PITCH*iy));
}

void FlowStep_SSE_diego::_sse_yrhsadd(const float * const f, float * const r)
{
	static const int SP = TempSOA::PITCH;
	static const int DP = OutputSOA::PITCH;
	
	for(int iy=0; iy<OutputSOA::NY; iy+=4)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
		{
			__m128 rhs0 = _mm_loadu_ps(f + iy +  1 + SP*ix) - _mm_load_ps(f+ iy + SP*ix);
			__m128 rhs1 = _mm_loadu_ps(f + iy +  1 + SP*ix+SP) - _mm_load_ps(f+ iy + SP*ix+SP);
			__m128 rhs2 = _mm_loadu_ps(f + iy +  1 + SP*ix+2*SP) - _mm_load_ps(f+ iy + SP*ix+2*SP);
			__m128 rhs3 = _mm_loadu_ps(f + iy +  1 + SP*ix+3*SP) - _mm_load_ps(f+ iy + SP*ix+3*SP);
			
			_MM_TRANSPOSE4_PS(rhs0, rhs1, rhs2, rhs3);
			
			_mm_store_ps(r + ix + DP*iy, _mm_load_ps(r + ix + DP*iy) + rhs0);
			_mm_store_ps(r + ix + DP*iy+DP, _mm_load_ps(r + ix + DP*iy+DP) + rhs1);
			_mm_store_ps(r + ix + DP*iy+2*DP, _mm_load_ps(r + ix + DP*iy+2*DP) + rhs2);
			_mm_store_ps(r + ix + DP*iy+3*DP, _mm_load_ps(r + ix + DP*iy+3*DP) + rhs3);
		}
}

void FlowStep_SSE_diego::_sse_yrhsadd(const double * const f, double * const r)
{
	static const int SP = TempSOA::PITCH;
	static const int DP = OutputSOA::PITCH;
	
	for(int iy=0; iy<OutputSOA::NY; iy+=2)
		for(int ix=0; ix<OutputSOA::NX; ix+=2)
		{
			const __m128d rhs0 = _mm_loadu_pd(f + iy +  1 + SP*ix) - _mm_load_pd(f+ iy + SP*ix);
			const __m128d rhs1 = _mm_loadu_pd(f + iy +  1 + SP*ix+SP) - _mm_load_pd(f+ iy + SP*ix+SP);
			
			const __m128d output0 = _mm_shuffle_pd(rhs0, rhs1, _MM_SHUFFLE2(0,0));
			const __m128d output1 = _mm_shuffle_pd(rhs0, rhs1, _MM_SHUFFLE2(1,1));
			
			_mm_store_pd(r + ix + DP*iy, _mm_load_pd(r + ix + DP*iy) + output0);
			_mm_store_pd(r + ix + DP*iy+DP, _mm_load_pd(r + ix + DP*iy+DP) + output1);
		}
}

void FlowStep_SSE_diego::_sse_zrhsadd(const float * const fb, const float * const ff, float * const r)
{	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
			_mm_store_ps(r + ix + OutputSOA::PITCH*iy, _mm_load_ps(r + ix + OutputSOA::PITCH*iy) + 
						 _mm_load_ps(ff + ix + TempSOA::PITCH*iy) - _mm_load_ps(fb + ix + TempSOA::PITCH*iy));
}

void FlowStep_SSE_diego::_sse_zrhsadd(const double * const fb, const double * const ff, double * const r)
{	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=2)
			_mm_store_pd(r + ix + OutputSOA::PITCH*iy, _mm_load_pd(r + ix + OutputSOA::PITCH*iy) + 
						 _mm_load_pd(ff + ix + TempSOA::PITCH*iy) - _mm_load_pd(fb + ix + TempSOA::PITCH*iy));
}

void FlowStep_SSE_diego::compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
								 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
{
	for(int islice=0; islice<5; islice++)
	{
		_convert(srcfirst+islice*srcfloats*slicesrcs, srcfloats, rowsrcs);
		_next();
	}
	
	_convert(srcfirst + 5*srcfloats*slicesrcs, srcfloats, rowsrcs);
	_zflux(-2);
	_flux_next();
	
	for(int islice=0; islice<_BLOCKSIZE_; islice++)
	{
		_xflux(-2);
		_xrhs();
		
		_yflux(-2);
		_yrhs();
		
		_next();
		_convert(srcfirst + (islice+6)*srcfloats*slicesrcs, srcfloats, rowsrcs);
		
		_zflux(-2);
		_zrhs();
		
		_copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
		_flux_next();
	}
}
