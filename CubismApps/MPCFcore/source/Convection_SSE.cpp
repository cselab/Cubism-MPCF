/*
 *  Convection_SSE.cpp
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

#define WENO_MINUS \
const __m128 is0 = a*(a*F_4_3  - b*F_19_3 + c*F_11_3)	+ b*(b*F_25_3 - c*F_31_3)	+ c*c*F_10_3	+ mywenoeps; \
const __m128 is1 = b*(b*F_4_3  - c*F_13_3 + d*F_5_3)	+ c*(c*F_13_3 - d*F_13_3)	+ d*d*F_4_3		+ mywenoeps; \
const __m128 is2 = c*(c*F_10_3 - d*F_31_3 + e*F_11_3) + d*(d*F_25_3  - e*F_19_3)	+ e*e*F_4_3		+ mywenoeps; \
WENO_OMEGAS \
const __m128 recdata = ( \
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

#include "Convection_SSE.h"

void Convection_SSE::_sse_convert_aligned(const float * const gptfirst, const int gptfloats, const int rowgpts,
											  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l)
{
	const __m128 F_1_2 = _mm_set_ps1(0.5f);
	const __m128 M_1_2 = _mm_set_ps1(-0.5f);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/smoothlength);
	const __m128 g1 = _mm_set_ps1(gamma1);
	const __m128 g2 = _mm_set_ps1(gamma2);
	const __m128 _pc1 = _mm_set_ps1(pc1);
	const __m128 _pc2 = _mm_set_ps1(pc2);
	
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
						 _getgamma(dataB1, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(dataB1, invsml, _pc1, _pc2, F_1, F_1_2, M_1_2));
		}
	}
	
#undef DESTID
}

void Convection_SSE::_sse_convert(const float * const gptfirst, const int gptfloats, const int rowgpts,
									  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l)
{
	const __m128 F_1_2 = _mm_set_ps1(0.5f);
	const __m128 M_1_2 = _mm_set_ps1(-0.5f);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/smoothlength);
	const __m128 g1 = _mm_set_ps1(gamma1);
	const __m128 g2 = _mm_set_ps1(gamma2);
	const __m128 _pc1 = _mm_set_ps1(pc1);
	const __m128 _pc2 = _mm_set_ps1(pc2);
	
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
						 _getgamma(dataB1, invsml, g1, g2, F_1, F_1_2, M_1_2)*_getPC(dataB1, invsml, _pc1, _pc2, F_1, F_1_2, M_1_2));
		}
	}
	
#undef DESTID
}

void Convection_SSE::_sse_xweno_minus(const float * const in, float * const out) const
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

void Convection_SSE::_sse_xweno_pluss(const float * const in, float * const out) const
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

void Convection_SSE::_sse_yweno_minus(const float * const in, float * const out)
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

void Convection_SSE::_sse_yweno_pluss(const float * const in, float * const out)
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

void Convection_SSE::_sse_zweno_minus(const float * const a_, const float * const b_, 
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

void Convection_SSE::_sse_zweno_pluss(const float * const b_, const float * const c_, 
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

void Convection_SSE::_sse_hlle_rho(const float * const rm, const float * const rp,
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

void Convection_SSE::_sse_hlle_vel(const float * const rm, const float * const rp,
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

void Convection_SSE::_sse_hlle_pvel(const float * const rm, const float * const rp,
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

void Convection_SSE::_sse_hlle_e(const float * const rm, const float * const rp,
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
	const __m128 invsml = _mm_set_ps1(1.f/smoothlength);
	const __m128 g1 = _mm_set_ps1(gamma1);
	const __m128 g2 = _mm_set_ps1(gamma2);
	
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

void Convection_SSE::_sse_char_vel(const float * const rm, const float * const rp, 
									   const float * const vm, const float * const vp,
									   const float * const pm, const float * const pp,
									   const float * const lm, const float * const lp, 
									   float * const outm, float * const outp)
{
	const __m128 F_1_2 = _mm_set_ps1(0.5);
	const __m128 M_1_2 = _mm_set_ps1(-0.5);
	const __m128 F_1 = _mm_set_ps1(1);
	const __m128 invsml = _mm_set_ps1(1.f/smoothlength);
	const __m128 g1 = _mm_set_ps1(gamma1);
	const __m128 g2 = _mm_set_ps1(gamma2);
	const __m128 _pc1 = _mm_set_ps1(pc1);
	const __m128 _pc2 = _mm_set_ps1(pc2);
	
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=4)
		{
#ifdef _PREC_DIV_
			const __m128 cminus = _mm_sqrt_ps(_getgamma(_mm_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											  _mm_max_ps((_mm_load_ps(pm + ID)+_getPC(_mm_load_ps(lm + ID), invsml, _pc1, _pc2, F_1, F_1_2, M_1_2))/_mm_load_ps(rm + ID), _mm_setzero_ps()));
			const __m128 cplus = _mm_sqrt_ps(_getgamma(_mm_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											 _mm_max_ps((_mm_load_ps(pp + ID)+_getPC(_mm_load_ps(lp + ID), invsml, _pc1, _pc2, F_1, F_1_2, M_1_2))/_mm_load_ps(rp + ID), _mm_setzero_ps()));
#else
			const __m128 cminus = worse_sqrt(_getgamma(_mm_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											 _mm_max_ps((_mm_load_ps(pm + ID)+_getPC(_mm_load_ps(lm + ID), invsml, _pc1, _pc2, F_1, F_1_2, M_1_2))* better_rcp(_mm_load_ps(rm + ID)), _mm_setzero_ps()));
			const __m128 cplus = worse_sqrt(_getgamma(_mm_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											_mm_max_ps((_mm_load_ps(pp + ID)+_getPC(_mm_load_ps(lp + ID), invsml, _pc1, _pc2, F_1, F_1_2, M_1_2))* better_rcp(_mm_load_ps(rp + ID)), _mm_setzero_ps()));
#endif
			_mm_store_ps(outm + ID, _mm_min_ps(_mm_load_ps(vm + ID) - cminus, _mm_load_ps(vm + ID) - cplus));
			_mm_store_ps(outp + ID, _mm_max_ps(_mm_load_ps(vp + ID) + cminus, _mm_load_ps(vp + ID) + cplus));
		}
#undef ID
}

void Convection_SSE::_sse_xrhsadd(const float * const f, float * const r)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
			_mm_store_ps(r + ix + OutputSOA::PITCH*iy,  
						 _mm_loadu_ps(f + ix + 1 + TempSOA::PITCH*iy) - _mm_load_ps(f + ix + TempSOA::PITCH*iy));
}

void Convection_SSE::_sse_yrhsadd(const float * const f, float * const r)
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

void Convection_SSE::_sse_zrhsadd(const float * const fb, const float * const ff, float * const r)
{	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=4)
			_mm_store_ps(r + ix + OutputSOA::PITCH*iy, _mm_load_ps(r + ix + OutputSOA::PITCH*iy) + 
						 _mm_load_ps(ff + ix + TempSOA::PITCH*iy) - _mm_load_ps(fb + ix + TempSOA::PITCH*iy));
}
