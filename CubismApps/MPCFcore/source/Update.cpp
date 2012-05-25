/*
 *  Update.cpp
 *  MPCFcore
 *
 *  Created by Babak Hejazialhosseini  on 6/9/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cassert>
#include <iostream>
#include <cstdlib>

#include "common.h"
#include "Update.h"

void Update_CPP::compute(const Real * const src, Real * const dst, const int gptfloats)
{
	assert(gptfloats >= 6);
	
	const int N=_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	
	for(int i=0; i<N; i+=gptfloats)
		for(int comp=0;comp<6;comp++)
			dst[i+comp] += m_b*src[i+comp];
}

#ifdef _SSE_
void Update_SSE::_sse_update_aligned(const float * const src, float * const dst, const int gptfloats)
{
	const int N=_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	const __m128 b = _mm_set1_ps(m_b);
	const __m128 b2 = _mm_set_ps(0, 0, m_b, m_b);
	
	for(int i=0; i<N; i+=gptfloats)
	{
		_mm_store_ps(dst+i,   _mm_load_ps(dst+i) + b*_mm_load_ps(src+i));
		_mm_store_ps(dst+i+4, _mm_load_ps(dst+i+4) + b2*_mm_load_ps(src+i+4));
	}	
}

void Update_SSE::_sse_update6(const float * const src, float * const dst, const int gptfloats)
{
	const int N=_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	const __m128 b = _mm_set1_ps(m_b);
	const __m128 b2 = _mm_set_ps(0, 0, m_b, m_b);
	
	assert(N % 2 == 0);
	
	for(int i=0; i<N; i+=gptfloats*2)
	{
		_mm_store_ps(dst+i,   _mm_load_ps(dst+i) + b*_mm_load_ps(src+i));
		_mm_store_ps(dst+i+4, _mm_load_ps(dst+i+4) + b2*_mm_load_ps(src+i+4));
		_mm_storeu_ps(dst+gptfloats+i,   _mm_loadu_ps(dst+gptfloats+i) + b*_mm_loadu_ps(src+gptfloats+i));
		_mm_storeu_ps(dst+gptfloats+i+4, _mm_loadu_ps(dst+gptfloats+i+4) + b2*_mm_loadu_ps(src+gptfloats+i+4));
	}	
}

void Update_SSE::_sse_update(const float * const src, float * const dst, const int gptfloats)
{
	const int Nm1=(_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_-1)*gptfloats;
	const __m128 b = _mm_set1_ps(m_b);
	const __m128 b2 = _mm_set_ps(0, 0, m_b, m_b);
	
	for(int i=0; i<Nm1; i+=gptfloats)
	{
		_mm_storeu_ps(dst+i,   _mm_loadu_ps(dst+i) + b*_mm_loadu_ps(src+i));
		_mm_storeu_ps(dst+i+4, _mm_loadu_ps(dst+i+4) + b2*_mm_loadu_ps(src+i+4));
	}
	
	{
		const int i = Nm1;
		_mm_storeu_ps(dst+i,   _mm_loadu_ps(dst+i) + b*_mm_loadu_ps(src+i));
		
		dst[i+4]+=m_b*src[i+4];
		dst[i+5]+=m_b*src[i+5];		
	}	
}

void Update_SSE::_sse_updateLD(const float * const src, float * const dst, const int gptfloats)
{
	const int Nfloats = _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	const int N = (Nfloats/4)*4;
	
	const __m128 b = _mm_set1_ps(m_b);
	
	for(int i=0; i<N; i+=4)
		_mm_store_ps(dst+i, _mm_load_ps(dst+i) + b*_mm_load_ps(src+i));
	
	if (Nfloats % 4 > 0)
		_mm_store_ps(dst+N, _mm_load_ps(dst+N) + _mm_set_ps(0, 0, m_b, m_b)*_mm_load_ps(src+N));
}
#endif

#ifdef _AVX_
void Update_AVX::_avx_sparse(const float * const src, float * const dst, const int gptfloats)
{
	const int Nfloats = _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	
	const __m256 b = _mm256_set_ps(0, 0, m_b, m_b, m_b, m_b, m_b, m_b);
	
	for(int i=0; i<Nfloats; i+=gptfloats)
		_mm256_storeu_ps(dst+i, _mm256_loadu_ps(dst+i) + b*_mm256_loadu_ps(src+i));
}

void Update_AVX::_avx_lockdown(const float * const src, float * const dst, const int gptfloats)
{
	const int Nfloats = _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	const int N = (Nfloats/8)*8;
	
	const __m256 b = _mm256_set1_ps(m_b);
	
	for(int i=0; i<N; i+=8)
		_mm256_storeu_ps(dst+i, _mm256_loadu_ps(dst+i) + b*_mm256_loadu_ps(src+i));
	
	//diego dai salta gio'! vo a pe'
	for(int i=N; i<Nfloats; i++) dst[i] += m_b*src[i];
}

void Update_AVX::_avx_lockdown32(const float * const src, float * const dst, const int gptfloats)
{
	const int Nfloats = _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*gptfloats;
	const int N = (Nfloats/8)*8;
	const int R = Nfloats % 8;
	
	const __m256 b = _mm256_set1_ps(m_b);
	
	for(int i=0; i<N; i+=8)
		_mm256_store_ps(dst+i, _mm256_load_ps(dst+i) + b*_mm256_load_ps(src+i));
	
	if (R > 0)
	{
		const __m256i mask = _mm256_castps_si256(_mm256_cmp_ps(_mm256_set_ps(7,6,5,4,3,2,1,0), _mm256_set1_ps(R), _CMP_LT_OS));
		_mm256_maskstore_ps(dst+N, mask, _mm256_maskload_ps(dst+N, mask) + b*_mm256_maskload_ps(src+N, mask));
	}
}
#endif
