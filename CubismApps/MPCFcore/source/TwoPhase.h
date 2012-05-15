/*
 *  TwoPhase.h
 *  
 *
 *  Created by Diego Rossinelli on 5/14/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

//C++ related functions
template<typename X> inline X mysqrt(X x){ abort(); return sqrt(x);}
template<>  inline float mysqrt<float>(float x){ return sqrtf(x);}
template<>  inline double mysqrt<double>(double x){ return sqrt(x);}

//SSE-related functions
#ifdef __INTEL_COMPILER
inline __m128 operator+(__m128 a, __m128 b){ return _mm_add_ps(a, b); }
inline __m128 operator|(__m128 a, __m128 b){ return _mm_or_ps(a, b); }
inline __m128 operator*(__m128 a, __m128 b){ return _mm_mul_ps(a, b); }
inline __m128 operator-(__m128 a,  __m128 b){ return _mm_sub_ps(a, b); }
inline __m128 operator/(__m128 a, __m128 b){ return _mm_div_ps(a, b); }
inline __m128d operator+(__m128d a, __m128d b){ return _mm_add_pd(a, b); }
inline __m128d operator*(__m128d a, __m128d b){ return _mm_mul_pd(a, b); }
inline __m128d operator-(__m128d a, __m128d b){ return _mm_sub_pd(a, b); }
inline __m128d operator/(__m128d a, __m128d b){ return _mm_div_pd(a, b); }
inline __m128d operator|(__m128d a, __m128d b){ return _mm_or_pd(a, b); }
#endif

inline __m128 better_rcp(const __m128 a)
{
	const __m128 Ra0 = _mm_rcp_ps(a);
	return _mm_sub_ps(_mm_add_ps(Ra0, Ra0), _mm_mul_ps(_mm_mul_ps(Ra0, a), Ra0));
}

inline __m128 better_rcp(const __m128 a)
{
	const __m128 Ra0 = _mm_rcp_ps(a);
	return _mm_sub_ps(_mm_add_ps(Ra0, Ra0), _mm_mul_ps(_mm_mul_ps(Ra0, a), Ra0));
}

inline __m128 worse_sqrt(const __m128 a)
{	
 	const __m128 invz =  _mm_rsqrt_ps(a);
	const __m128 z = _mm_rcp_ps(invz);
	const __m128 tmp =  z-(z*z-a)*invz*_mm_set_ps1(0.5f);
	return  _mm_and_ps(tmp, _mm_cmpgt_ps(a, _mm_setzero_ps()));
}

inline __m128 _heaviside(const __m128 phi, const __m128 inv_h, const __m128 one, const __m128 phalf, const __m128 mhalf) const
{
	const __m128 x = _mm_min_ps(one, _mm_max_ps(_mm_setzero_ps() - one, phi*inv_h));
	
	const __m128 val_xneg = (mhalf*x - one)*x + phalf;
	const __m128 val_xpos = (phalf*x - one)*x + phalf;
	
	const __m128 flag = _mm_cmplt_ps(x, _mm_setzero_ps());
	
	return _mm_or_ps(_mm_and_ps(flag, val_xneg),_mm_andnot_ps(flag, val_xpos));
}

inline __m128d _heaviside(const __m128d phi, const __m128d inv_h, const __m128d one, const __m128d phalf, const __m128d mhalf) const
{
	const __m128d x = _mm_min_pd(one, _mm_max_pd(_mm_setzero_pd() - one, phi*inv_h));
	
	const __m128d val_xneg = (mhalf*x - one)*x + phalf;
	const __m128d val_xpos = (phalf*x - one)*x + phalf;
	
	const __m128d flag = _mm_cmplt_pd(x, _mm_setzero_pd());
	
	return _mm_or_pd(_mm_and_pd(flag, val_xneg),_mm_andnot_pd(flag, val_xpos));
}

inline __m128 _getgamma(const __m128 phi, const __m128 inv_smoothlength, 
											const __m128 gamma1, const __m128 gamma2,
											const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const
{
	const __m128 hs = _heaviside(phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	return (gamma1)*hs + (gamma2)*(F_1-hs); 
}

inline __m128d _getgamma(const __m128d phi, const __m128d inv_smoothlength, 
											 const __m128d gamma1, const __m128d gamma2,
											 const __m128d F_1, const __m128d F_1_2, const __m128d M_1_2) const
{
	const __m128d hs = _heaviside(phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	return (gamma1)*hs + (gamma2)*(F_1-hs); 
}

//AVX-related functions
#define _MM_TRANSPOSE8_PS(row0, row1, row2, row3, row4, row5, row6, row7)\
{\
__m256 __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;\
__m128 __f0, __f1, __f2, __f3;\
__t0 = _mm256_unpacklo_ps(row0,row1);\
__t1 = _mm256_unpacklo_ps(row2,row3);\
__t2 = _mm256_unpackhi_ps(row0,row1);\
__t3 = _mm256_unpackhi_ps(row2,row3);\
__t4 = _mm256_unpacklo_ps(row4,row5);\
__t5 = _mm256_unpacklo_ps(row6,row7);\
__t6 = _mm256_unpackhi_ps(row4,row5);\
__t7 = _mm256_unpackhi_ps(row6,row7);\
row0 = _mm256_shuffle_ps(__t0,__t1,_MM_SHUFFLE(1,0,1,0));\
row1 = _mm256_shuffle_ps(__t0,__t1,_MM_SHUFFLE(3,2,3,2));\
row2 = _mm256_shuffle_ps(__t2,__t3,_MM_SHUFFLE(1,0,1,0));\
row3 = _mm256_shuffle_ps(__t2,__t3,_MM_SHUFFLE(3,2,3,2));\
row4 = _mm256_shuffle_ps(__t4,__t5,_MM_SHUFFLE(1,0,1,0));\
row5 = _mm256_shuffle_ps(__t4,__t5,_MM_SHUFFLE(3,2,3,2));\
row6 = _mm256_shuffle_ps(__t6,__t7,_MM_SHUFFLE(1,0,1,0));\
row7 = _mm256_shuffle_ps(__t6,__t7,_MM_SHUFFLE(3,2,3,2));\
__f0 = _mm256_extractf128_ps(row0,1);\
__f1 = _mm256_extractf128_ps(row1,1);\
__f2 = _mm256_extractf128_ps(row2,1);\
__f3 = _mm256_extractf128_ps(row3,1);\
row0 = _mm256_insertf128_ps(row0, _mm256_extractf128_ps(row4,0), 1);\
row1 = _mm256_insertf128_ps(row1, _mm256_extractf128_ps(row5,0), 1);\
row2 = _mm256_insertf128_ps(row2, _mm256_extractf128_ps(row6,0), 1);\
row3 = _mm256_insertf128_ps(row3, _mm256_extractf128_ps(row7,0), 1);\
row4 = _mm256_insertf128_ps(row4, __f0, 0);\
row5 = _mm256_insertf128_ps(row5, __f1, 0);\
row6 = _mm256_insertf128_ps(row6, __f2, 0);\
row7 = _mm256_insertf128_ps(row7, __f3, 0);\
}

#define _MM_TRANSPOSE8x6_PS(row0, row1, row2, row3, row4, row5, row6, row7)\
{\
__m256 __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;\
__m128 __f0, __f1;\
__t0 = _mm256_unpacklo_ps(row0,row1);\
__t1 = _mm256_unpacklo_ps(row2,row3);\
__t2 = _mm256_unpackhi_ps(row0,row1);\
__t3 = _mm256_unpackhi_ps(row2,row3);\
__t4 = _mm256_unpacklo_ps(row4,row5);\
__t5 = _mm256_unpacklo_ps(row6,row7);\
__t6 = _mm256_unpackhi_ps(row4,row5);\
__t7 = _mm256_unpackhi_ps(row6,row7);\
row0 = _mm256_shuffle_ps(__t0,__t1,_MM_SHUFFLE(1,0,1,0));\
row1 = _mm256_shuffle_ps(__t0,__t1,_MM_SHUFFLE(3,2,3,2));\
row2 = _mm256_shuffle_ps(__t2,__t3,_MM_SHUFFLE(1,0,1,0));\
row3 = _mm256_shuffle_ps(__t2,__t3,_MM_SHUFFLE(3,2,3,2));\
row4 = _mm256_shuffle_ps(__t4,__t5,_MM_SHUFFLE(1,0,1,0));\
row5 = _mm256_shuffle_ps(__t4,__t5,_MM_SHUFFLE(3,2,3,2));\
row6 = _mm256_shuffle_ps(__t6,__t7,_MM_SHUFFLE(1,0,1,0));\
row7 = _mm256_shuffle_ps(__t6,__t7,_MM_SHUFFLE(3,2,3,2));\
__f0 = _mm256_extractf128_ps(row0,1);\
__f1 = _mm256_extractf128_ps(row1,1);\
row0 = _mm256_insertf128_ps(row0, _mm256_extractf128_ps(row4,0), 1);\
row1 = _mm256_insertf128_ps(row1, _mm256_extractf128_ps(row5,0), 1);\
row2 = _mm256_insertf128_ps(row2, _mm256_extractf128_ps(row6,0), 1);\
row3 = _mm256_insertf128_ps(row3, _mm256_extractf128_ps(row7,0), 1);\
row4 = _mm256_insertf128_ps(row4, __f0, 0);\
row5 = _mm256_insertf128_ps(row5, __f1, 0);\
}

inline __m256 better_rcp(const __m256 a)
{
	const __m256 Ra0 = _mm256_rcp_ps(a);
	return _mm256_sub_ps(_mm256_add_ps(Ra0, Ra0), _mm256_mul_ps(_mm256_mul_ps(Ra0, a), Ra0));
}

inline __m256 worse_sqrt(const __m256 a)
{
	const __m256 invz =  _mm256_rsqrt_ps(a);
	const __m256 z = _mm256_rcp_ps(invz);
	const __m256 tmp = z-(z*z-a)*invz*_mm256_set1_ps(0.5f);
	
	return  _mm256_and_ps(tmp, _mm256_cmp_ps(a, _mm256_setzero_ps(), _CMP_GT_OS));
}