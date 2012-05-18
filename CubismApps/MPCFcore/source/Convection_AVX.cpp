/*
 *  FlowStep_AVX_diego.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/3/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include "Convection_AVX.h"

#define WENO_CONSTANTS \
const __m256 one = _mm256_set1_ps(1.); \
const __m256 F_4_3 = _mm256_set1_ps(4./3.); \
const __m256 F_5_3 = _mm256_set1_ps(5./3.); \
const __m256 F_10_3 = _mm256_set1_ps(10./3.); \
const __m256 F_11_3 = _mm256_set1_ps(11./3.); \
const __m256 F_13_3 = _mm256_set1_ps(13./3.); \
const __m256 F_19_3 = _mm256_set1_ps(19./3.);  \
const __m256 F_25_3 = _mm256_set1_ps(25./3.); \
const __m256 F_31_3 = _mm256_set1_ps(31./3.); \
const __m256 mywenoeps = _mm256_set1_ps(WENOEPS); \
\
const __m256 M_1_6 = _mm256_set1_ps(-1./6.); \
const __m256 F_1_3 = _mm256_set1_ps(1./3.); \
const __m256 F_5_6 = _mm256_set1_ps(5./6.); \
const __m256 M_7_6 = _mm256_set1_ps(-7./6.); \
const __m256 F_11_6 = _mm256_set1_ps(11./6.); \
\
const __m256 F_1_10 = _mm256_set1_ps(1./10.); \
const __m256 F_6_10 = _mm256_set1_ps(6./10.); \
const __m256 F_3_10 = _mm256_set1_ps(3./10.);

#ifndef _PREC_DIV_

#ifndef _PREC_DIV_NONE_

#define WENO_OMEGAS \
const __m256 alpha0 = F_1_10 / (is0*is0); \
const __m256 alpha1 = F_6_10 / (is1*is1); \
const __m256 alpha2 = F_3_10 / (is2*is2); \
const __m256 inv_alpha = better_rcp(alpha0+alpha1+alpha2); \
const __m256 omega0=alpha0 * inv_alpha; \
const __m256 omega1=alpha1 * inv_alpha; \
const __m256 omega2= one-omega0-omega1; 

#else

#define WENO_OMEGAS \
const __m256 alpha0 = F_1_10 * better_rcp(is0*is0); \
const __m256 alpha1 = F_6_10 * better_rcp(is1*is1); \
const __m256 alpha2 = F_3_10 * better_rcp(is2*is2); \
const __m256 inv_alpha = better_rcp(alpha0+alpha1+alpha2); \
const __m256 omega0=alpha0 * inv_alpha; \
const __m256 omega1=alpha1 * inv_alpha; \
const __m256 omega2= one-omega0-omega1; 

#endif

#else 

#define WENO_OMEGAS \
const __m256 alpha0 = F_1_10 / (is0*is0); \
const __m256 alpha1 = F_6_10 / (is1*is1); \
const __m256 alpha2 = F_3_10 / (is2*is2); \
const __m256 inv_alpha = one / (alpha0+alpha1+alpha2); \
const __m256 omega0 = alpha0 * inv_alpha; \
const __m256 omega1 = alpha1 * inv_alpha; \
const __m256 omega2 = one - omega0 - omega1; 

#endif

#define WENO_MINUS \
const __m256 is0 = a*(a*F_4_3  - b*F_19_3 + c*F_11_3)	+ b*(b*F_25_3 - c*F_31_3)	+ c*c*F_10_3	+ mywenoeps; \
const __m256 is1 = b*(b*F_4_3  - c*F_13_3 + d*F_5_3)	+ c*(c*F_13_3 - d*F_13_3)	+ d*d*F_4_3		+ mywenoeps; \
const __m256 is2 = c*(c*F_10_3 - d*F_31_3 + e*F_11_3) + d*(d*F_25_3  - e*F_19_3)	+ e*e*F_4_3		+ mywenoeps; \
WENO_OMEGAS \
const __m256 recdata = ( \
omega0 *(F_1_3*a + M_7_6*b + F_11_6*c)+   \
omega1 *(M_1_6*b + F_5_6*c + F_1_3*d) +  \
omega2 *(F_1_3*c + F_5_6*d + M_1_6*e)	);

#define WENO_PLUS \
const __m256 is0 = (d*(d*F_10_3 - e*F_31_3 + f*F_11_3)	+ e*(e*F_25_3  - f*F_19_3)	+ f*f*F_4_3	)	+ mywenoeps; \
const __m256 is1 = (c*(c*F_4_3  - d*F_13_3 + e*F_5_3 )	+ d*(d*F_13_3 - e*F_13_3)	+ e*e*F_4_3	)	+ mywenoeps; \
const __m256 is2 = (b*(b*F_4_3  - c*F_19_3 + d*F_11_3) + c*(c*F_25_3  - d*F_31_3)	+ d*d*F_10_3)	+ mywenoeps; \
WENO_OMEGAS \
const __m256 recdata = ( \
omega0 *(F_1_3*f + M_7_6*e + F_11_6*d)+   \
omega1 *(M_1_6*e + F_5_6*d + F_1_3*c) +  \
omega2 *(F_1_3*d + F_5_6*c + M_1_6*b)	);

inline void _MM_TRANSPOSE8_PS(__m256& row0, __m256& row1, __m256& row2, __m256& row3, __m256& row4, __m256& row5, __m256& row6, __m256& row7)
{
	__m256 __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;
	__m128 __f0, __f1, __f2, __f3;
	__t0 = _mm256_unpacklo_ps(row0,row1);
	__t1 = _mm256_unpacklo_ps(row2,row3);
	__t2 = _mm256_unpackhi_ps(row0,row1);
	__t3 = _mm256_unpackhi_ps(row2,row3);
	__t4 = _mm256_unpacklo_ps(row4,row5);
	__t5 = _mm256_unpacklo_ps(row6,row7);
	__t6 = _mm256_unpackhi_ps(row4,row5);
	__t7 = _mm256_unpackhi_ps(row6,row7);
	row0 = _mm256_shuffle_ps(__t0,__t1,_MM_SHUFFLE(1,0,1,0));
	row1 = _mm256_shuffle_ps(__t0,__t1,_MM_SHUFFLE(3,2,3,2));
	row2 = _mm256_shuffle_ps(__t2,__t3,_MM_SHUFFLE(1,0,1,0));
	row3 = _mm256_shuffle_ps(__t2,__t3,_MM_SHUFFLE(3,2,3,2));
	row4 = _mm256_shuffle_ps(__t4,__t5,_MM_SHUFFLE(1,0,1,0));
	row5 = _mm256_shuffle_ps(__t4,__t5,_MM_SHUFFLE(3,2,3,2));
	row6 = _mm256_shuffle_ps(__t6,__t7,_MM_SHUFFLE(1,0,1,0));
	row7 = _mm256_shuffle_ps(__t6,__t7,_MM_SHUFFLE(3,2,3,2));
	__f0 = _mm256_extractf128_ps(row0,1);
	__f1 = _mm256_extractf128_ps(row1,1);
	__f2 = _mm256_extractf128_ps(row2,1);
	__f3 = _mm256_extractf128_ps(row3,1);
	row0 = _mm256_insertf128_ps(row0, _mm256_extractf128_ps(row4,0), 1);
	row1 = _mm256_insertf128_ps(row1, _mm256_extractf128_ps(row5,0), 1);
	row2 = _mm256_insertf128_ps(row2, _mm256_extractf128_ps(row6,0), 1);
	row3 = _mm256_insertf128_ps(row3, _mm256_extractf128_ps(row7,0), 1);
	row4 = _mm256_insertf128_ps(row4, __f0, 0);
	row5 = _mm256_insertf128_ps(row5, __f1, 0);
	row6 = _mm256_insertf128_ps(row6, __f2, 0);
	row7 = _mm256_insertf128_ps(row7, __f3, 0);
}

inline void _MM_TRANSPOSE8x6_PS(__m256& row0, __m256& row1, __m256& row2, __m256& row3, __m256& row4, __m256& row5, __m256& row6, __m256& row7)
{
	__m256 __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;
	__m128 __f0, __f1;
	__t0 = _mm256_unpacklo_ps(row0,row1);
	__t1 = _mm256_unpacklo_ps(row2,row3);
	__t2 = _mm256_unpackhi_ps(row0,row1);
	__t3 = _mm256_unpackhi_ps(row2,row3);
	__t4 = _mm256_unpacklo_ps(row4,row5);
	__t5 = _mm256_unpacklo_ps(row6,row7);
	__t6 = _mm256_unpackhi_ps(row4,row5);
	__t7 = _mm256_unpackhi_ps(row6,row7);
	row0 = _mm256_shuffle_ps(__t0,__t1,_MM_SHUFFLE(1,0,1,0));
	row1 = _mm256_shuffle_ps(__t0,__t1,_MM_SHUFFLE(3,2,3,2));
	row2 = _mm256_shuffle_ps(__t2,__t3,_MM_SHUFFLE(1,0,1,0));
	row3 = _mm256_shuffle_ps(__t2,__t3,_MM_SHUFFLE(3,2,3,2));
	row4 = _mm256_shuffle_ps(__t4,__t5,_MM_SHUFFLE(1,0,1,0));
	row5 = _mm256_shuffle_ps(__t4,__t5,_MM_SHUFFLE(3,2,3,2));
	row6 = _mm256_shuffle_ps(__t6,__t7,_MM_SHUFFLE(1,0,1,0));
	row7 = _mm256_shuffle_ps(__t6,__t7,_MM_SHUFFLE(3,2,3,2));
	__f0 = _mm256_extractf128_ps(row0,1);
	__f1 = _mm256_extractf128_ps(row1,1);
	row0 = _mm256_insertf128_ps(row0, _mm256_extractf128_ps(row4,0), 1);
	row1 = _mm256_insertf128_ps(row1, _mm256_extractf128_ps(row5,0), 1);
	row2 = _mm256_insertf128_ps(row2, _mm256_extractf128_ps(row6,0), 1);
	row3 = _mm256_insertf128_ps(row3, _mm256_extractf128_ps(row7,0), 1);
	row4 = _mm256_insertf128_ps(row4, __f0, 0);
	row5 = _mm256_insertf128_ps(row5, __f1, 0);
}

void Convection_AVX::_avx_convert_aligned(const float * const gptfirst, const int gptfloats, const int rowgpts, const int slicegpts,
											  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l)
{	
	const __m256 F_1_2 = _mm256_set1_ps(0.5f);
	const __m256 M_1_2 = _mm256_set1_ps(-0.5f);
	const __m256 F_1 = _mm256_set1_ps(1);
	const __m256 invsml = _mm256_set1_ps(1.f/smoothlength);
	const __m256 g1 = _mm256_set1_ps(gamma1);
	const __m256 g2 = _mm256_set1_ps(gamma2);
	
#define DESTID (dx + (InputSOA::PITCH)*dy)
#define SRCID dx*gptfloats
	for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
	{
		const float * const in = gptfirst + dy*gptfloats*rowgpts -5*gptfloats;
		
		for(int dx=0; dx<_BLOCKSIZE_+16; dx+=8)
		{
			__m256 A,B,C,D,E,F,G,H;
			
			if (dx>0) A = _mm256_load_ps(in + SRCID);
			if (dx>0) B = _mm256_load_ps(in + SRCID + gptfloats);
			if (dx>0) C = _mm256_load_ps(in + SRCID + 2*gptfloats);
			if (dx>0 && dx<_BLOCKSIZE_+8) D = _mm256_load_ps(in + SRCID + 3*gptfloats);
			if (dx>0 && dx<_BLOCKSIZE_+8) E = _mm256_load_ps(in + SRCID + 4*gptfloats);
			if (dx<_BLOCKSIZE_+8) F = _mm256_load_ps(in + SRCID + 5*gptfloats);
			if (dx<_BLOCKSIZE_+8) G = _mm256_load_ps(in + SRCID + 6*gptfloats);
			if (dx<_BLOCKSIZE_+8) H = _mm256_load_ps(in + SRCID + 7*gptfloats);
			
			_MM_TRANSPOSE8x6_PS(A,B,C,D,E,F,G,H);
			
			_mm256_store_ps(rho + DESTID, A);
#ifdef _PREC_DIV_
			const __m256 inv_rho = F_1/A;			
#else
			const __m256 inv_rho = better_rcp(A);
#endif
			_mm256_store_ps(u + DESTID, B*inv_rho);
			_mm256_store_ps(v + DESTID, C*inv_rho);
			_mm256_store_ps(w + DESTID, D*inv_rho);
			
			_mm256_store_ps(l + DESTID, F);
			
			_mm256_store_ps(p + DESTID,  
							(E - (B*B + C*C + D*D)*(F_1_2*inv_rho))*
							(_getgamma(F, invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1));
		}
	}
	
#undef DESTID
#undef SRCID
}

void Convection_AVX::_avx_convert(const float * const gptfirst, const int gptfloats, const int rowgpts, const int slicegpts,
									  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l)
{
	const __m256 F_1_2 = _mm256_set1_ps(0.5f);
	const __m256 M_1_2 = _mm256_set1_ps(-0.5f);
	const __m256 F_1 = _mm256_set1_ps(1);
	const __m256 invsml = _mm256_set1_ps(1.f/smoothlength);
	const __m256 g1 = _mm256_set1_ps(gamma1);
	const __m256 g2 = _mm256_set1_ps(gamma2);
	
#define DESTID (dx + (InputSOA::PITCH)*dy)
#define SRCID dx*gptfloats
	for(int dy=0; dy<_BLOCKSIZE_+6; dy++)
	{
		const float * const in = gptfirst + dy*gptfloats*rowgpts -5*gptfloats;
		
		for(int dx=0; dx<_BLOCKSIZE_+16; dx+=8)
		{
			__m256 A,B,C,D,E,F,G,H;
			
			if (dx>0) A = _mm256_loadu_ps(in + SRCID);
			if (dx>0) B = _mm256_loadu_ps(in + SRCID + gptfloats);
			if (dx>0) C = _mm256_loadu_ps(in + SRCID + 2*gptfloats);
			if (dx>0 && dx<_BLOCKSIZE_+8) D = _mm256_loadu_ps(in + SRCID + 3*gptfloats);
			if (dx>0 && dx<_BLOCKSIZE_+8) E = _mm256_loadu_ps(in + SRCID + 4*gptfloats);
			if (dx<_BLOCKSIZE_+8) F = _mm256_loadu_ps(in + SRCID + 5*gptfloats);
			if (dx<_BLOCKSIZE_+8) G = _mm256_loadu_ps(in + SRCID + 6*gptfloats);
			if (dx<_BLOCKSIZE_+8) H = _mm256_loadu_ps(in + SRCID + 7*gptfloats);
			
			_MM_TRANSPOSE8x6_PS(A,B,C,D,E,F,G,H);
			
			_mm256_store_ps(rho + DESTID, A);
#ifdef _PREC_DIV_
			const __m256 inv_rho = F_1/A;			
#else
			const __m256 inv_rho = better_rcp(A);
#endif
			_mm256_store_ps(u + DESTID, B*inv_rho);
			_mm256_store_ps(v + DESTID, C*inv_rho);
			_mm256_store_ps(w + DESTID, D*inv_rho);
			
			_mm256_store_ps(l + DESTID, F);
			
			_mm256_store_ps(p + DESTID,  
							(E - (B*B + C*C + D*D)*(F_1_2*inv_rho))*
							(_getgamma(F, invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1));
		}
	}
	
#undef DESTID
#undef SRCID
}

void Convection_AVX::_avx_xweno_minus(const float * const in, float * const out) const
{	
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=8)
		{	
			
#define SL1(L,R) _mm256_shuffle_ps(L, _mm256_shuffle_ps(L,R, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1))
#define SL2(L,R) _mm256_shuffle_ps(L, R, _MM_SHUFFLE(1,0,3,2))
#define SL3(L,R) _mm256_shuffle_ps(_mm256_shuffle_ps(L,R, _MM_SHUFFLE(0,0,3,3)), R, _MM_SHUFFLE(2,1,3,0))
			
			//i rather prefer some access latency than register spills			
#define W _mm_load_ps(&in[dx + 4 + SX*dy])
#define C _mm256_load_ps(&in[dx + 8 + SX*dy])
#define E _mm_load_ps(&in[dx + 16 + SX*dy])
			
#define C0W _mm256_insertf128_ps(C, W, 1)
#define WC0 _mm256_permute2f128_ps(C0W, C0W, 1)
			
#define EC1 _mm256_insertf128_ps(C, E, 0)
#define C1E _mm256_permute2f128_ps(EC1, EC1, 1)
			
#define a  SL1(WC0, C)
#define b  SL2(WC0, C)
#define c  SL3(WC0, C)
#define d  C
#define e  SL1(C, C1E)
			
			WENO_MINUS
			
			_mm256_store_ps(&out[dx + TempSOA::PITCH*dy], recdata);
#undef a
#undef b
#undef c
#undef d
#undef e
			
#undef W
#undef C
#undef E
#undef C0W
#undef WC0
#undef EC1
#undef C1E
			
#undef SL1
#undef SL2
#undef SL3
			
		}
}

void Convection_AVX::_avx_xweno_pluss(const float * const in, float * const out) const
{
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=8)
		{	
			
#define SL1(L,R) _mm256_shuffle_ps(L, _mm256_shuffle_ps(L,R, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1))
#define SL2(L,R) _mm256_shuffle_ps(L, R, _MM_SHUFFLE(1,0,3,2))
#define SL3(L,R) _mm256_shuffle_ps(_mm256_shuffle_ps(L,R, _MM_SHUFFLE(0,0,3,3)), R, _MM_SHUFFLE(2,1,3,0))
			
			//i rather prefer some cache latency than register spills			
#define W _mm_load_ps(&in[dx + 4 + SX*dy])
#define C _mm256_load_ps(&in[dx + 8 + SX*dy])
#define E _mm_load_ps(&in[dx + 16 + SX*dy])
#define C0W _mm256_insertf128_ps(C, W, 1)
#define WC0 _mm256_permute2f128_ps(C0W, C0W, 1)
#define EC1 _mm256_insertf128_ps(C, E, 0)
#define C1E _mm256_permute2f128_ps(EC1, EC1, 1)
			
#define b  SL2(WC0, C)
#define c  SL3(WC0, C)
#define d  C
#define e  SL1(C, C1E)
#define f  SL2(C, C1E)
			
			WENO_PLUS
			
			_mm256_store_ps(&out[dx + TempSOA::PITCH*dy], recdata);
#undef b
#undef c
#undef d
#undef e
#undef f
			
#undef W
#undef C
#undef E
#undef C0W
#undef WC0
#undef EC1
#undef C1E
			
#undef SL1
#undef SL2
#undef SL3
			
		}
}

void Convection_AVX::_avx_yweno_minus(const float * const in, float * const out)
{	
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy+=8)
	{
		const float * ptr = &in[dy];
		
		for(int dx=0; dx<TempSOA::NX; dx++)
		{	
			//i rather prefer some cache latency than register spills			
#define a _mm256_load_ps(ptr + dx*SX)
#define b _mm256_load_ps(ptr + dx*SX + SX)
#define c _mm256_load_ps(ptr + dx*SX + 2*SX)
#define d _mm256_load_ps(ptr + dx*SX + 3*SX)
#define e _mm256_load_ps(ptr + dx*SX + 4*SX)
			
			WENO_MINUS
#undef a
#undef b
#undef c
#undef d
#undef e
			_mm256_store_ps(&tmp[dx][0], recdata);
		}
		
		for(int dx=0; dx<TempSOA::NX-1; dx+=8)
		{
			__m256 data0 = _mm256_load_ps(&tmp[dx][0]);
			__m256 data1 = _mm256_load_ps(&tmp[dx+1][0]);
			__m256 data2 = _mm256_load_ps(&tmp[dx+2][0]);
			__m256 data3 = _mm256_load_ps(&tmp[dx+3][0]);
			__m256 data4 = _mm256_load_ps(&tmp[dx+4][0]);
			__m256 data5 = _mm256_load_ps(&tmp[dx+5][0]);
			__m256 data6 = _mm256_load_ps(&tmp[dx+6][0]);
			__m256 data7 = _mm256_load_ps(&tmp[dx+7][0]);
			
			_MM_TRANSPOSE8_PS(data0, data1, data2, data3, data4, data5, data6, data7);
			
			_mm256_store_ps(&out[dx + TempSOA::PITCH*dy], data0);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+1)], data1);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+2)], data2);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+3)], data3);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+4)], data4);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+5)], data5);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+6)], data6);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+7)], data7);
		}
		
		{
			out[TempSOA::NX-1 + TempSOA::PITCH*dy] = tmp[TempSOA::NX-1][0];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+1)] = tmp[TempSOA::NX-1][1];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+2)] = tmp[TempSOA::NX-1][2];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+3)] = tmp[TempSOA::NX-1][3];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+4)] = tmp[TempSOA::NX-1][4];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+5)] = tmp[TempSOA::NX-1][5];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+6)] = tmp[TempSOA::NX-1][6];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+7)] = tmp[TempSOA::NX-1][7];
		}
	}
}

void Convection_AVX::_avx_yweno_pluss(const float * const in, float * const out)
{	
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy+=8)
	{
		const float * ptr = &in[dy];
		
		for(int dx=0; dx<TempSOA::NX; dx++)
		{
			//i rather prefer some cache latency than register spills
#define b _mm256_load_ps(ptr + dx*SX + SX)
#define c _mm256_load_ps(ptr + dx*SX + 2*SX)
#define d _mm256_load_ps(ptr + dx*SX + 3*SX)
#define e _mm256_load_ps(ptr + dx*SX + 4*SX)
#define f _mm256_load_ps(ptr + dx*SX + 5*SX)
			
			WENO_PLUS
#undef b
#undef c
#undef d
#undef e
#undef f
			_mm256_store_ps(&tmp[dx][0], recdata);
		}
		
		for(int dx=0; dx<TempSOA::NX-1; dx+=8)
		{
			__m256 data0 = _mm256_load_ps(&tmp[dx][0]);
			__m256 data1 = _mm256_load_ps(&tmp[dx+1][0]);
			__m256 data2 = _mm256_load_ps(&tmp[dx+2][0]);
			__m256 data3 = _mm256_load_ps(&tmp[dx+3][0]);
			__m256 data4 = _mm256_load_ps(&tmp[dx+4][0]);
			__m256 data5 = _mm256_load_ps(&tmp[dx+5][0]);
			__m256 data6 = _mm256_load_ps(&tmp[dx+6][0]);
			__m256 data7 = _mm256_load_ps(&tmp[dx+7][0]);
			
			_MM_TRANSPOSE8_PS(data0, data1, data2, data3, data4, data5, data6, data7);
			
			_mm256_store_ps(&out[dx + TempSOA::PITCH*dy], data0);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+1)], data1);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+2)], data2);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+3)], data3);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+4)], data4);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+5)], data5);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+6)], data6);
			_mm256_store_ps(&out[dx + TempSOA::PITCH*(dy+7)], data7);
		}
		
		{
			out[TempSOA::NX-1 + TempSOA::PITCH*dy] = tmp[TempSOA::NX-1][0];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+1)] = tmp[TempSOA::NX-1][1];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+2)] = tmp[TempSOA::NX-1][2];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+3)] = tmp[TempSOA::NX-1][3];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+4)] = tmp[TempSOA::NX-1][4];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+5)] = tmp[TempSOA::NX-1][5];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+6)] = tmp[TempSOA::NX-1][6];
			out[TempSOA::NX-1 + TempSOA::PITCH*(dy+7)] = tmp[TempSOA::NX-1][7];
		}
	}
}

void Convection_AVX::_avx_zweno_minus(const float * const a_, const float * const b_, 
										  const float * const c_, const float * const d_, 
										  const float * const e_ , float * const out) const
{
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=8)
		{	
			//i rather prefer some cache latency than register spills
#define a _mm256_load_ps(a_ + dx + SX*dy)
#define b _mm256_load_ps(b_ + dx + SX*dy)
#define c _mm256_load_ps(c_ + dx + SX*dy)
#define d _mm256_load_ps(d_ + dx + SX*dy)
#define e _mm256_load_ps(e_ + dx + SX*dy)
			
			WENO_MINUS
			
			_mm256_store_ps(&out[dx + TempSOA::PITCH*dy], recdata);
			
#undef a
#undef b
#undef c
#undef d
#undef e
			
		}	
}

void Convection_AVX::_avx_zweno_pluss(const float * const b_, const float * const c_, 
										  const float * const d_, const float * const e_, 
										  const float * const f_ , float * const out) const
{
	WENO_CONSTANTS
	
	static const int SX = InputSOA::PITCH;
	
	for(int dy=0; dy<TempSOA::NY; dy++)
		for(int dx=0; dx<TempSOA::NX; dx+=8)
		{		
			//i rather prefer some cache latency than register spills
#define b _mm256_load_ps(b_ + dx + SX*dy)
#define c _mm256_load_ps(c_ + dx + SX*dy)
#define d _mm256_load_ps(d_ + dx + SX*dy)
#define e _mm256_load_ps(e_ + dx + SX*dy)
#define f _mm256_load_ps(f_ + dx + SX*dy)
			
			WENO_PLUS
			
			_mm256_store_ps(&out[dx + TempSOA::PITCH*dy], recdata);
			
#undef b
#undef c
#undef d
#undef e
#undef f
			
		}	
}

void Convection_AVX::_avx_hlle_rho(const float * const rm, const float * const rp,
									   const float * const vm, const float * const vp,
									   const float * const am, const float * const ap,
									   float * const out)
{
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)	
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=8)
		{
#define aminus _mm256_load_ps(am + ID)
#define apluss  _mm256_load_ps(ap + ID)
#define flagminus _mm256_cmp_ps(aminus, _mm256_setzero_ps(), _CMP_NLE_US)
#define flagplus  _mm256_cmp_ps(apluss, _mm256_setzero_ps(), _CMP_LT_OS)
#define rminus _mm256_load_ps(rm + ID)
#define rpluss _mm256_load_ps(rp + ID)
			
			const __m256 fminus = _mm256_load_ps(vm + ID)*rminus;
			const __m256 fpluss = _mm256_load_ps(vp + ID)*rpluss;
			
#ifdef _PREC_DIV_						
			const __m256 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(rpluss-rminus))/(apluss-aminus);
#else
			const __m256 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(rpluss-rminus))*better_rcp(apluss-aminus);
#endif
			_mm256_store_ps(out + ID, _mm256_blendv_ps(_mm256_blendv_ps(fother, fpluss, flagplus), fminus, flagminus));
			
#undef rminus
#undef rpluss
#undef aminus
#undef apluss
#undef flagminus
#undef flagplus
		}
	
#undef ID
}

void Convection_AVX::_avx_hlle_vel(const float * const rm, const float * const rp,
									   const float * const vm, const float * const vp,
									   const float * const vdm, const float * const vdp,
									   const float * const am, const float * const ap,
									   float * const out)
{
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)	
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=8)
		{			
			const __m256 uminus = _mm256_load_ps(vm + ID)*_mm256_load_ps(rm + ID);
			const __m256 upluss = _mm256_load_ps(vp + ID)*_mm256_load_ps(rp + ID);
			
			const __m256 fminus = _mm256_load_ps(vdm + ID)*uminus;
			const __m256 fpluss = _mm256_load_ps(vdp + ID)*upluss;
			
#define aminus _mm256_load_ps(am + ID)
#define apluss  _mm256_load_ps(ap + ID)
#define flagminus _mm256_cmp_ps(aminus, _mm256_setzero_ps(), _CMP_NLE_US)
#define flagplus  _mm256_cmp_ps(apluss, _mm256_setzero_ps(), _CMP_LT_OS)
			
#ifdef _PREC_DIV_						
			const __m256 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))/(apluss-aminus);
#else
			const __m256 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))*better_rcp(apluss-aminus);
#endif
			_mm256_store_ps(out + ID, _mm256_blendv_ps(_mm256_blendv_ps(fother, fpluss, flagplus), fminus, flagminus));
			
#undef aminus
#undef apluss
#undef flagminus
#undef flagplus
		}
	
#undef ID
}

void Convection_AVX::_avx_hlle_pvel(const float * const rm, const float * const rp,
										const float * const vm, const float * const vp,
										const float * const pm, const float * const pp,
										const float * const am, const float * const ap,
										float * const out)
{
	static const int P = TempSOA::PITCH;
#define ID (ix + P*iy)
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=8)
		{
			const __m256 myvminus = _mm256_load_ps(vm + ID);
			const __m256 myvpluss = _mm256_load_ps(vp + ID);
			
			const __m256 uminus = myvminus*_mm256_load_ps(rm + ID);
			const __m256 upluss = myvpluss*_mm256_load_ps(rp + ID);
			
			const __m256 fminus = myvminus*uminus + _mm256_load_ps(pm + ID);
			const __m256 fpluss = myvpluss*upluss + _mm256_load_ps(pp + ID);
			
#define aminus _mm256_load_ps(am + ID)
#define apluss  _mm256_load_ps(ap + ID)
#define flagminus _mm256_cmp_ps(aminus, _mm256_setzero_ps(), _CMP_NLE_US)
#define flagplus  _mm256_cmp_ps(apluss, _mm256_setzero_ps(), _CMP_LT_OS)
			
#ifdef _PREC_DIV_						
			const __m256 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))/(apluss-aminus);
#else
			const __m256 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))*better_rcp(apluss-aminus);
#endif
			_mm256_store_ps(out + ID, _mm256_blendv_ps(_mm256_blendv_ps(fother, fpluss, flagplus), fminus, flagminus));
			
#undef aminus
#undef apluss
#undef flagminus
#undef flagplus
			
		}
	
#undef ID
}

void Convection_AVX::_avx_hlle_e(const float * const rm, const float * const rp,
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
	
	const __m256 F_1_2 = _mm256_set1_ps(0.5);
	const __m256 M_1_2 = _mm256_set1_ps(-0.5);
	const __m256 F_1 = _mm256_set1_ps(1);
	const __m256 invsml = _mm256_set1_ps(1.f/smoothlength);
	const __m256 g1 = _mm256_set1_ps(gamma1);
	const __m256 g2 = _mm256_set1_ps(gamma2);
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix+=8)
		{			
			//i rather prefer some cache latency than register spills
#define vdminus  _mm256_load_ps(vdm + ID)
#define v1minus  _mm256_load_ps(v1m + ID)
#define v2minus  _mm256_load_ps(v2m + ID)
#define pminus   _mm256_load_ps(pm + ID)
			
#ifdef _PREC_DIV_
			const __m256 eminus = pminus*(F_1/(_getgamma(_mm256_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm256_load_ps(rm + ID)*(vdminus*vdminus + v1minus*v1minus + v2minus*v2minus);
#else
			const __m256 eminus = pminus*(F_1 * better_rcp(_getgamma(_mm256_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm256_load_ps(rm + ID)*(vdminus*vdminus + v1minus*v1minus + v2minus*v2minus);
#endif
			//i rather prefer some cache latency than register spills
#define vdplus _mm256_load_ps(vdp + ID)
#define v1plus _mm256_load_ps(v1p + ID)
#define v2plus _mm256_load_ps(v2p + ID)
#define pplus  _mm256_load_ps(pp + ID)
			
#ifdef _PREC_DIV_			
			const __m256 eplus = pplus*(F_1/(_getgamma(_mm256_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm256_load_ps(rp + ID)*(vdplus*vdplus + v1plus*v1plus + v2plus*v2plus);
#else
			const __m256 eplus = pplus*(F_1 * better_rcp(_getgamma(_mm256_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)-F_1)) +
			F_1_2*_mm256_load_ps(rp + ID)*(vdplus*vdplus + v1plus*v1plus + v2plus*v2plus);
#endif
			const __m256 fminus = vdminus*(pminus + eminus);
			const __m256 fpluss = vdplus *(pplus + eplus);
			
#define aminus _mm256_load_ps(am + ID)
#define apluss  _mm256_load_ps(ap + ID)
#define flagminus _mm256_cmp_ps(aminus, _mm256_setzero_ps(), _CMP_NLE_US)
#define flagplus  _mm256_cmp_ps(apluss, _mm256_setzero_ps(), _CMP_LT_OS)
			
#ifdef _PREC_DIV_						
			const __m256 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(eplus-eminus))/(apluss-aminus);
#else
			const __m256 fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(eplus-eminus))*better_rcp(apluss-aminus);
#endif		
			_mm256_store_ps(out + ID, _mm256_blendv_ps(_mm256_blendv_ps(fother, fpluss, flagplus), fminus, flagminus));
			
#undef vdplus
#undef v1plus
#undef v2plus
#undef pplus
#undef vdplus
#undef v1plus
#undef v2plus
#undef pplus
#undef aminus
#undef apluss
#undef flagminus
#undef flagplus
			
		}
#undef ID
}

void Convection_AVX::_avx_char_vel(const float * const rm, const float * const rp, 
									   const float * const vm, const float * const vp,
									   const float * const pm, const float * const pp,
									   const float * const lm, const float * const lp, 
									   float * const outm, float * const outp)
{
	const __m256 F_1_2 = _mm256_set1_ps(0.5);
	const __m256 M_1_2 = _mm256_set1_ps(-0.5);
	const __m256 F_1 = _mm256_set1_ps(1);
	const __m256 invsml = better_rcp(_mm256_set1_ps(smoothlength));
	const __m256 g1 = _mm256_set1_ps(gamma1);
	const __m256 g2 = _mm256_set1_ps(gamma2);
	
	static const int P = TempSOA::PITCH;
	
#define ID (ix + P*iy)
	
	for(int iy=0; iy<TempSOA::NY; iy++)
#pragma ivdep
		for(int ix=0; ix<TempSOA::NX; ix+=8)
		{
#ifdef _PREC_DIV_
			const __m256 cminus = _mm256_sqrt_ps(_getgamma(_mm256_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
												 _mm256_max_ps(_mm256_load_ps(pm + ID)/_mm256_load_ps(rm + ID), _mm256_setzero_ps()));
			const __m256 cplus = _mm256_sqrt_ps(_getgamma(_mm256_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
												_mm256_max_ps(_mm256_load_ps(pp + ID)/_mm256_load_ps(rp + ID), _mm256_setzero_ps()));
#else
			const __m256 cminus = worse_sqrt(_getgamma(_mm256_load_ps(lm + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											 _mm256_max_ps(_mm256_load_ps(pm + ID)* better_rcp(_mm256_load_ps(rm + ID)), _mm256_setzero_ps()));
			const __m256 cplus = worse_sqrt(_getgamma(_mm256_load_ps(lp + ID), invsml, g1, g2, F_1, F_1_2, M_1_2)* 
											_mm256_max_ps(_mm256_load_ps(pp + ID)* better_rcp(_mm256_load_ps(rp + ID)), _mm256_setzero_ps()));
#endif
			_mm256_store_ps(outm + ID, _mm256_min_ps(_mm256_load_ps(vm + ID) - cminus, _mm256_load_ps(vm + ID) - cplus));
			_mm256_store_ps(outp + ID, _mm256_max_ps(_mm256_load_ps(vp + ID) + cminus, _mm256_load_ps(vp + ID) + cplus));
		}
#undef ID
}

void Convection_AVX::_avx_xrhsadd(const float * const f, float * const r)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=8)
			_mm256_store_ps(r + ix + OutputSOA::PITCH*iy,  
							_mm256_loadu_ps(f + ix + 1 + TempSOA::PITCH*iy) - _mm256_load_ps(f + ix + TempSOA::PITCH*iy));
}

void Convection_AVX::_avx_yrhsadd(const float * const f, float * const r)
{
	static const int SP = TempSOA::PITCH;
	static const int DP = OutputSOA::PITCH;
	
	for(int iy=0; iy<OutputSOA::NY; iy+=8)
		for(int ix=0; ix<OutputSOA::NX; ix+=8)
		{
			__m256 rhs0 = _mm256_loadu_ps(f + iy +  1 + SP*ix) - _mm256_load_ps(f+ iy + SP*ix);
			__m256 rhs1 = _mm256_loadu_ps(f + iy +  1 + SP*ix+SP) - _mm256_load_ps(f+ iy + SP*ix+SP);
			__m256 rhs2 = _mm256_loadu_ps(f + iy +  1 + SP*ix+2*SP) - _mm256_load_ps(f+ iy + SP*ix+2*SP);
			__m256 rhs3 = _mm256_loadu_ps(f + iy +  1 + SP*ix+3*SP) - _mm256_load_ps(f+ iy + SP*ix+3*SP);
			__m256 rhs4 = _mm256_loadu_ps(f + iy +  1 + SP*ix+4*SP) - _mm256_load_ps(f+ iy + SP*ix+4*SP);
			__m256 rhs5 = _mm256_loadu_ps(f + iy +  1 + SP*ix+5*SP) - _mm256_load_ps(f+ iy + SP*ix+5*SP);
			__m256 rhs6 = _mm256_loadu_ps(f + iy +  1 + SP*ix+6*SP) - _mm256_load_ps(f+ iy + SP*ix+6*SP);
			__m256 rhs7 = _mm256_loadu_ps(f + iy +  1 + SP*ix+7*SP) - _mm256_load_ps(f+ iy + SP*ix+7*SP);
			
			_MM_TRANSPOSE8_PS(rhs0, rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7);
			
			_mm256_store_ps(r + ix + DP*iy, _mm256_load_ps(r + ix + DP*iy) + rhs0);
			_mm256_store_ps(r + ix + DP*iy+DP, _mm256_load_ps(r + ix + DP*iy+DP) + rhs1);
			_mm256_store_ps(r + ix + DP*iy+2*DP, _mm256_load_ps(r + ix + DP*iy+2*DP) + rhs2);
			_mm256_store_ps(r + ix + DP*iy+3*DP, _mm256_load_ps(r + ix + DP*iy+3*DP) + rhs3);
			_mm256_store_ps(r + ix + DP*iy+4*DP, _mm256_load_ps(r + ix + DP*iy+4*DP) + rhs4);
			_mm256_store_ps(r + ix + DP*iy+5*DP, _mm256_load_ps(r + ix + DP*iy+5*DP) + rhs5);
			_mm256_store_ps(r + ix + DP*iy+6*DP, _mm256_load_ps(r + ix + DP*iy+6*DP) + rhs6);
			_mm256_store_ps(r + ix + DP*iy+7*DP, _mm256_load_ps(r + ix + DP*iy+7*DP) + rhs7);
		}
}

void Convection_AVX::_avx_zrhsadd(const float * const fb, const float * const ff, float * const r)
{	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix+=8)
			_mm256_store_ps(r + ix + OutputSOA::PITCH*iy, _mm256_load_ps(r + ix + OutputSOA::PITCH*iy) + 
							_mm256_load_ps(ff + ix + TempSOA::PITCH*iy) - _mm256_load_ps(fb + ix + TempSOA::PITCH*iy));
}
