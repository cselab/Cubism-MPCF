/*
 *  Update.h
 *  MPCFcore
 *
 *  Created by Babak Hejazialhosseini  on 6/9/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "common.h"

class Update_CPP
{
protected:
	
	Real m_b;
	
	bool _is_aligned(const void * const ptr, unsigned int alignment)
	{
		return ((size_t)ptr) % alignment == 0;
	}
	
public:
	
	Update_CPP(Real b=1): m_b(b) {}
	
	void compute(const Real * const src, Real * const dst, const int gptfloats);
	
	static void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const size_t NCORES, const size_t NT, const size_t NBLOCKS, const float MEASUREDTIME, const bool bAwk=false)
	{
		const float PEAKPERF = PEAKPERF_CORE*NCORES;
		
		//FLOP estimation
		const double GFLOPUPDATE = 6.*2.e-9*NBLOCKS*1.*(_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_);
		
		//FLOP/s estimation
		const float OIUpdate = 2.f/(3*sizeof(Real));
		const float AIUpdate = 2.f/(4*sizeof(Real));
		const double EPERFUPDATE   = min(OIUpdate*   PEAKBAND, PEAKPERF);
		const double EPERFUPDATEAI   = min(AIUpdate*   PEAKBAND, PEAKPERF);
		
		//execution time estimation
		const double TUPDATE = 1.e9*GFLOPUPDATE/EPERFUPDATE;
		const double TUPDATEAI = 1.e9*GFLOPUPDATE/EPERFUPDATEAI;
		
		printPerformanceTitle();
		printf("\tUP TOTAL GFLOPS: %.2f\n", GFLOPUPDATE);
		printf("\tUP ASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tUP RIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tUP THIS ONE IS %.2f GFLOP/s,\t\"PER Update\" %.2f FLOP/B\n", GFLOPUPDATE/MEASUREDTIME, OIUpdate);
		printf("\tUP TIME PER BLOCK: %.5f ms (expected %.5f ms)\n",  1e3*MEASUREDTIME/NBLOCKS, 1e3*TUPDATE/NBLOCKS);
		printf("\tUP Expected Performance is: %.2f GFLOP/s [AI], %.2f GFLOP/s [OI]\n", GFLOPUPDATE/TUPDATEAI, GFLOPUPDATE/TUPDATE);
		printf("\tUP EFFICIENCY: %.2f%% [AI] - %.2f%% [OI], HW-UTILIZATION: %.2f%%\n", 100.*TUPDATEAI/MEASUREDTIME, 100.*TUPDATE/MEASUREDTIME, 100*(GFLOPUPDATE/MEASUREDTIME*1e9)/PEAKPERF);
		printEndLine();
	}
};

class Update_SSE : public Update_CPP
{	
protected:
	void _sse_update_aligned(const float * const src, float * const dst, const int gptfloats);
	void _sse_update6(const float * const src, float * const dst, const int gptfloats);
	void _sse_update(const float * const src, float * const dst, const int gptfloats);
	void _sse_updateLD(const float * const src, float * const dst, const int gptfloats);
	
public:
	
	Update_SSE(Real b=1): Update_CPP(b) { }
	
	void compute(const Real * const src, Real * const dst, const int gptfloats)
	{
#if defined(_SP_COMP_)
		assert(gptfloats >= 6);
		
		const bool bAligned = _is_aligned(src, 16) && _is_aligned(dst, 16);
		const bool b4Multiple = gptfloats % 4 == 0;
		
		if (bAligned && gptfloats == 6)
			_sse_updateLD(src, dst, gptfloats);
		else if (bAligned && b4Multiple)
			_sse_update_aligned(src, dst, gptfloats);
		else if (bAligned && gptfloats % 6 == 0)
			_sse_update6(src, dst, gptfloats);
		else
			_sse_update(src, dst, gptfloats);
#else
		Update_CPP::compute(src, dst, gptfloats);
#endif
	}
};

#ifdef _AVX_
class Update_AVX : public Update_CPP
{	
protected:
	void _avx_sparse(const float * const src, float * const dst, const int gptfloats);
	void _avx_lockdown(const float * const src, float * const dst, const int gptfloats);
	void _avx_lockdown32(const float * const src, float * const dst, const int gptfloats);
	
public:
	
	Update_AVX(Real b=1): Update_CPP(b) { }
	
	void compute(const float * const src, float * const dst, const int gptfloats)
	{
		assert(gptfloats >= 6);
		
		const bool aligned32B = _is_aligned(src, 32) && _is_aligned(dst, 32);
		
		const bool unitstride = gptfloats == 6;
		
		if (aligned32B && unitstride)
			_avx_lockdown32(src, dst, gptfloats);
		else if (unitstride)
			_avx_lockdown(src, dst, gptfloats);	
		else
			_avx_sparse(src, dst, gptfloats);
	}
};
#endif
