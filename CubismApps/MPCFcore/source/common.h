/*
 *  TwoPhase.h
 *  
 *
 *  Created by Diego Rossinelli on 5/14/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

inline void printKernelName(string kernelname)
{
	cout << endl;
	for (int i=0; i<100; i++)
		cout << "=";
	cout << endl << endl;
	cout << "KERNEL - " << kernelname.c_str() << endl << endl;
}

inline void printEndKernelTest()
{
	cout << endl;
	for (int i=0; i<100; i++)
		cout << "=";
	cout << endl << endl;
}

inline void printAccuracyTitle()
{
	string acc = " ACCURACY ";
	int l = (80-acc.size())/2;
	cout << "\t";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << acc.c_str();
	cout << "";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << endl;
}

inline void printPerformanceTitle()
{
	string perf = " PERFORMANCE ";
	int l = (80-perf.size())/2;
	cout << "\t";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << perf.c_str();
	cout << "";
	for (int i=0; i<80-l-(int)perf.size(); i++)
		cout << "-";
	cout << endl;
}

inline void printEndLine()
{
	cout << "\t";
	for (int i=0; i<80; i++)
		cout << "-";
	cout << endl;
}

inline void awkAcc(string kernelname,
				   double linf0, double linf1, double linf2, double linf3, double linf4, double linf5,
				   double l10, double l11, double l12, double l13, double l14, double l15)
{
	// should also output name
	cout << "awkAccLinf\t" <<  kernelname.c_str();
	cout << "\t" << setprecision(4) << linf0;
	cout << "\t" << setprecision(4) << linf1;
	cout << "\t" << setprecision(4) << linf2;
	cout << "\t" << setprecision(4) << linf3;
	cout << "\t" << setprecision(4) << linf4;
	cout << "\t" << setprecision(4) << linf5;
	cout << endl;
	
	cout << "awkAccL1\t" <<  kernelname.c_str();
	cout << "\t" << setprecision(4) << l10;
	cout << "\t" << setprecision(4) << l11;
	cout << "\t" << setprecision(4) << l12;
	cout << "\t" << setprecision(4) << l13;
	cout << "\t" << setprecision(4) << l14;
	cout << "\t" << setprecision(4) << l15;
	cout << endl;
}

//C++ related functions
template<typename X> inline X mysqrt(X x){ abort(); return sqrt(x);}
template<>  inline float mysqrt<float>(float x){ return sqrtf(x);}
template<>  inline double mysqrt<double>(double x){ return sqrt(x);}

template<typename X> inline X myabs(X x){ abort(); return sqrt(x);}
template<>  inline float myabs<float>(float x){ return fabsf(x);}
template<>  inline double myabs<double>(double x){ return fabs(x);}


//SSE-related functions
#include <xmmintrin.h>

#ifdef __INTEL_COMPILER
inline __m128 operator+(__m128 a, __m128 b){ return _mm_add_ps(a, b); }
inline __m128 operator&(__m128 a, __m128 b){ return _mm_and_ps(a, b); }
inline __m128 operator|(__m128 a, __m128 b){ return _mm_or_ps(a, b); }
inline __m128 operator*(__m128 a, __m128 b){ return _mm_mul_ps(a, b); }
inline __m128 operator-(__m128 a,  __m128 b){ return _mm_sub_ps(a, b); }
inline __m128 operator/(__m128 a, __m128 b){ return _mm_div_ps(a, b); }
inline __m128d operator+(__m128d a, __m128d b){ return _mm_add_pd(a, b); }
inline __m128d operator*(__m128d a, __m128d b){ return _mm_mul_pd(a, b); }
inline __m128d operator-(__m128d a, __m128d b){ return _mm_sub_pd(a, b); }
inline __m128d operator/(__m128d a, __m128d b){ return _mm_div_pd(a, b); }
inline __m128d operator&(__m128d a, __m128d b){ return _mm_and_pd(a, b); }
inline __m128d operator|(__m128d a, __m128d b){ return _mm_or_pd(a, b); }
#endif


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

inline __m128 worse_sqrt_unsafe(const __m128 a)
{
	const __m128 invz =  _mm_rsqrt_ps(a);
	const __m128 z = _mm_rcp_ps(invz);
	return z-(z*z-a)*invz*_mm_set_ps1(0.5f);
}

inline __m128 myrcp(const __m128 a) 
{
#ifdef _PREC_DIV_
	return _mm_set1_ps(1.f)/a;
#else
	return better_rcp(a);
#endif
}

inline __m128 mysqrt(const __m128 a) 
{
#ifdef _PREC_DIV_
	return _mm_sqrt_ps(a);
#else
	const __m128 invz =  _mm_rsqrt_ps(a);
	const __m128 z = _mm_rcp_ps(invz);
	return z-(z*z-a)*invz*_mm_set_ps1(0.5f);
#endif
}

inline __m128 myrsqrt(const __m128 v) 
{
#ifdef _PREC_DIV_
	return _mm_div_ps(_mm_set1_ps(1.), _mm_sqrt_ps(v));
#else
	const __m128 approx = _mm_rsqrt_ps( v );
	const __m128 muls = _mm_mul_ps(_mm_mul_ps(v, approx), approx);
	return _mm_mul_ps(_mm_mul_ps(_mm_set1_ps(0.5f), approx), _mm_sub_ps(_mm_set1_ps(3.f), muls) );
#endif
}

//AVX-related functions
#ifdef _AVX_
#include <immintrin.h>

#if defined(__INTEL_COMPILER)
inline __m256 operator+(__m256 a, __m256 b){ return _mm256_add_ps(a, b); }
inline __m256 operator&(__m256 a, __m256 b){ return _mm256_and_ps(a, b); }
inline __m256 operator|(__m256 a, __m256 b){ return _mm256_or_ps(a, b); }
inline __m256 operator*(__m256 a, __m256 b){ return _mm256_mul_ps(a, b); }
inline __m256 operator-(__m256 a, __m256 b){ return _mm256_sub_ps(a, b); }
inline __m256 operator/(__m256 a, __m256 b){ return _mm256_div_ps(a, b); }
inline __m256d operator+(__m256d a, __m256d b){ return _mm256_add_pd(a, b); }
inline __m256d operator*(__m256d a, __m256d b){ return _mm256_mul_pd(a, b); }
inline __m256d operator-(__m256d a, __m256d b){ return _mm256_sub_pd(a, b); }
inline __m256d operator/(__m256d a, __m256d b){ return _mm256_div_pd(a, b); }
inline __m256d operator&(__m256d a, __m256d b){ return _mm256_and_pd(a, b); }
inline __m256d operator|(__m256d a, __m256d b){ return _mm256_or_pd(a, b); }
#endif

inline __m256 heaviside(const __m256 phi, const __m256 inv_h, const __m256 one, const __m256 phalf, const __m256 mhalf) const
{
	const __m256 x = _mm256_min_ps(one, _mm256_max_ps(_mm256_setzero_ps() - one, phi*inv_h));
	
	return (_mm256_blendv_ps(phalf, mhalf,  _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OS))*x - one)*x + phalf;	
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

inline __m256 myrsqrt(const __m256 v) const
{
#ifdef _PREC_DIV_
	return _mm256_div_ps(_mm256_set1_ps(1.), _mm_sqrt_ps(v));
#else
	const __m256 approx = _mm256_rsqrt_ps( v );
	const __m256 muls = _mm256_mul_ps(_mm256_mul_ps(v, approx), approx);
	return _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.5f), approx), _mm256_sub_ps(_mm256_set1_ps(3.f), muls) );
#endif
}

#endif
