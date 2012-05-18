/*
 *  Convection_SSEd.cpp
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

#define _3ORPD_(a,b,c) _mm_or_pd(a, _mm_or_pd(b,c))

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

#define WENO_OMEGAS_DP \
const __m128d alpha0 = F_1_10 / (is0*is0); \
const __m128d alpha1 = F_6_10 / (is1*is1); \
const __m128d alpha2 = F_3_10 / (is2*is2); \
const __m128d inv_alpha = one /(alpha0+alpha1+alpha2); \
const __m128d omega0=alpha0 * inv_alpha; \
const __m128d omega1=alpha1 * inv_alpha; \
const __m128d omega2= one - omega0 - omega1; 

#define WENO_MINUS_DP \
const __m128d is0 = a*(a*F_4_3  - b*F_19_3 + c*F_11_3)	+ b*(b*F_25_3 - c*F_31_3)	+ c*c*F_10_3	+ mywenoeps; \
const __m128d is1 = b*(b*F_4_3  - c*F_13_3 + d*F_5_3)	+ c*(c*F_13_3 - d*F_13_3)	+ d*d*F_4_3		+ mywenoeps; \
const __m128d is2 = c*(c*F_10_3 - d*F_31_3 + e*F_11_3) + d*(d*F_25_3  - e*F_19_3)	+ e*e*F_4_3		+ mywenoeps; \
WENO_OMEGAS_DP \
const __m128d recdata = ( \
omega0 *(F_1_3*a + M_7_6*b + F_11_6*c)+   \
omega1 *(M_1_6*b + F_5_6*c + F_1_3*d) +  \
omega2 *(F_1_3*c + F_5_6*d + M_1_6*e)	);

#define WENO_PLUS_DP \
const __m128d is0 = (d*(d*F_10_3 - e*F_31_3 + f*F_11_3)	+ e*(e*F_25_3  - f*F_19_3)	+ f*f*F_4_3	)	+ mywenoeps; \
const __m128d is1 = (c*(c*F_4_3  - d*F_13_3 + e*F_5_3 )	+ d*(d*F_13_3 - e*F_13_3)	+ e*e*F_4_3	)	+ mywenoeps; \
const __m128d is2 = (b*(b*F_4_3  - c*F_19_3 + d*F_11_3) + c*(c*F_25_3  - d*F_31_3)	+ d*d*F_10_3)	+ mywenoeps; \
WENO_OMEGAS_DP \
const __m128d recdata = ( \
omega0 *(F_1_3*f + M_7_6*e + F_11_6*d)+   \
omega1 *(M_1_6*e + F_5_6*d + F_1_3*c) +  \
omega2 *(F_1_3*d + F_5_6*c + M_1_6*b)	);

#include "Convection_SSE.h"

void Convection_SSE::_sse_convert(const double * const gptfirst, const int gptfloats, const int rowgpts,
								  double * const rho, double * const u, double * const v, double * const w, double * const p, double * const l)
{
	const __m128d F_1_2 = _mm_set_pd1(0.5);
	const __m128d M_1_2 = _mm_set_pd1(-0.5);
	const __m128d F_1 = _mm_set_pd1(1);
	const __m128d invsml = _mm_set_pd1(1./smoothlength);
	const __m128d g1 = _mm_set_pd1(gamma1);
	const __m128d g2 = _mm_set_pd1(gamma2);
	const __m128d _pc1 = _mm_set_pd1(pc1);
	const __m128d _pc2 = _mm_set_pd1(pc2);
	
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
						 _getgamma(myl, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(myl, invsml, _pc1, _pc2, F_1, F_1_2, M_1_2));
			
			_mm_store_pd(l + DESTID, myl);
		}
	}
#undef DESTID
}

void Convection_SSE::_sse_convert_aligned(const double * const gptfirst, const int gptfloats, const int rowgpts,
										  double * const rho, double * const u, double * const v, double * const w, double * const p, double * const l)
{
	const __m128d F_1_2 = _mm_set_pd1(0.5);
	const __m128d M_1_2 = _mm_set_pd1(-0.5);
	const __m128d F_1 = _mm_set_pd1(1);
	const __m128d invsml = _mm_set_pd1(1./smoothlength);
	const __m128d g1 = _mm_set_pd1(gamma1);
	const __m128d g2 = _mm_set_pd1(gamma2);
	const __m128d _pc1 = _mm_set_pd1(pc1);
	const __m128d _pc2 = _mm_set_pd1(pc2);
	
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
						 _getgamma(myl, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(myl, invsml, _pc1, _pc2, F_1, F_1_2, M_1_2));
			
			_mm_store_pd(l + DESTID, myl);
		}
	}
#undef DESTID
}

void Convection_SSE::_sse_xweno_minus(const double * const in, double * const out) const
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

void Convection_SSE::_sse_xweno_pluss(const double * const in, double * const out) const
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

void Convection_SSE::_sse_yweno_minus(const double * const in, double * const out)
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

void Convection_SSE::_sse_yweno_pluss(const double * const in, double * const out)
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

void Convection_SSE::_sse_zweno_minus(const double * const a_, const double * const b_, 
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

void Convection_SSE::_sse_zweno_pluss(const double * const b_, const double * const c_, 
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

void Convection_SSE::_sse_hlle_rho(const double * const rm, const double * const rp,
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

void Convection_SSE::_sse_hlle_vel(const double * const rm, const double * const rp,
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

void Convection_SSE::_sse_hlle_pvel(const double * const rm, const double * const rp,
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

void Convection_SSE::_sse_hlle_e(const double * const rm, const double * const rp,
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
	const __m128d invsml = _mm_set_pd1(1/smoothlength);
	const __m128d g1 = _mm_set_pd1(gamma1);
	const __m128d g2 = _mm_set_pd1(gamma2);
	
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

void Convection_SSE::_sse_char_vel(const double * const rm, const double * const rp, 
								   const double * const vm, const double * const vp,
								   const double * const pm, const double * const pp,
								   const double * const lm, const double * const lp, 
								   double * const outm, double * const outp)
{
	const __m128d F_1_2 = _mm_set_pd1(0.5);
	const __m128d M_1_2 = _mm_set_pd1(-0.5);
	const __m128d F_1 = _mm_set_pd1(1);
	const __m128d invsml = _mm_set_pd1(1/smoothlength);
	const __m128d g1 = _mm_set_pd1(gamma1);
	const __m128d g2 = _mm_set_pd1(gamma2);
	const __m128d _pc1 = _mm_set_pd1(pc1);
	const __m128d _pc2 = _mm_set_pd1(pc2);
	
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=2)
		{
			const __m128d cminus = _mm_sqrt_pd(_getgamma(_mm_load_pd(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											   _mm_max_pd((_mm_load_pd(pm + ID)+_getPC(_mm_load_pd(lm + ID), invsml, _pc1, _pc2, F_1, F_1_2, M_1_2))/_mm_load_pd(rm + ID), _mm_setzero_pd()));
			const __m128d cplus = _mm_sqrt_pd(_getgamma(_mm_load_pd(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											  _mm_max_pd((_mm_load_pd(pp + ID)+_getPC(_mm_load_pd(lp + ID), invsml, _pc1, _pc2, F_1, F_1_2, M_1_2))/_mm_load_pd(rp + ID), _mm_setzero_pd()));
			_mm_store_pd(outm + ID, _mm_min_pd(_mm_load_pd(vm + ID) - cminus, _mm_load_pd(vm + ID) - cplus));
			_mm_store_pd(outp + ID, _mm_max_pd(_mm_load_pd(vp + ID) + cminus, _mm_load_pd(vp + ID) + cplus));
		}
#undef ID
}

void Convection_SSE::_sse_xrhsadd(const double * const f, double * const r)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=2)
			_mm_store_pd(r + ix + OutputSOA::PITCH*iy,  
						 _mm_loadu_pd(f + ix + 1 + TempSOA::PITCH*iy) - _mm_load_pd(f + ix + TempSOA::PITCH*iy));
}

void Convection_SSE::_sse_yrhsadd(const double * const f, double * const r)
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

void Convection_SSE::_sse_zrhsadd(const double * const fb, const double * const ff, double * const r)
{	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=2)
			_mm_store_pd(r + ix + OutputSOA::PITCH*iy, _mm_load_pd(r + ix + OutputSOA::PITCH*iy) + 
						 _mm_load_pd(ff + ix + TempSOA::PITCH*iy) - _mm_load_pd(fb + ix + TempSOA::PITCH*iy));
}
