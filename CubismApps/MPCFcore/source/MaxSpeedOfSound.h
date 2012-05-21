/*
 *  MaxSpeedOfSound.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cstdio>		
#include <algorithm>

#include "common.h"
#include "SOA2D.h"

class MaxSpeedOfSound_CPP
{
protected:
	
	const Real gamma1, gamma2, smoothlength;
    const Real pc1, pc2;
	
	inline Real _getgamma(const Real phi) const { return reconstruct(gamma1, gamma2, phi, 1/smoothlength); } 
	inline Real _getPC(const Real phi) const { return reconstruct(pc1, pc2, phi, 1/smoothlength); }
	
public:
	
	MaxSpeedOfSound_CPP(const Real gamma1, const Real gamma2, const Real smoothlength, const Real pc1, const Real pc2):
		gamma1(gamma1), gamma2(gamma2), smoothlength(smoothlength), pc1(pc1), pc2(pc2) { }
	
	Real compute(const Real * const src, const int gptfloats);
	
	static void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NT, const int NBLOCKS, float MEASUREDTIME, const bool bAwk=false)
	{
		const float PEAKPERF = PEAKPERF_CORE*NCORES;
		
		//FLOP estimation
		const double GFLOPSOS = NBLOCKS*_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*122.*1.e-9; 
		
		//FLOP/s estimation
		const float OISOS = 122./(6*sizeof(float));
		const float AISOS = 122./(6*sizeof(float)*3);
		const double EPERFSOS = min(OISOS*PEAKBAND, PEAKPERF);
		const double EPERFSOSAI = min(AISOS*PEAKBAND, PEAKPERF);

		const double TSOS = 1.e9*GFLOPSOS/EPERFSOS;
		
		printPerformanceTitle();
		printf("\tSOS TOTAL GFLOPS: %.2f\n", GFLOPSOS);
		printf("\tSOS ASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tSOS RIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tSOS THIS ONE IS %.2f GFLOP/s,\t\"PER SOS\" %.2f FLOP/B\n", GFLOPSOS/MEASUREDTIME, OISOS);
		printf("\tSOS TIME PER BLOCK: %.5f ms (expected %.5f ms)\n",  1e3*MEASUREDTIME/NBLOCKS, 1e3*TSOS/NBLOCKS);
		printf("\tSOS Expected Performance is: %.2f GFLOP/s [AI], %.2f GFLOP/s [OI]\n", GFLOPSOS/(1.e9*GFLOPSOS/EPERFSOSAI), GFLOPSOS/TSOS);
		printf("\tSOS EFFICIENCY: %.2f%% [AI] - %.2f%% [OI], HW-UTILIZATION: %.2f%%\n", 100.*(1.e9*GFLOPSOS/EPERFSOSAI)/MEASUREDTIME, 100.*TSOS/MEASUREDTIME, 100*(GFLOPSOS/MEASUREDTIME*1e9)/PEAKPERF);		
		printEndLine();
	}
};

#if defined(_SSE_) 
class MaxSpeedOfSound_SSE: public MaxSpeedOfSound_CPP
{		
	SOA2D<0, _BLOCKSIZE_, 0, _BLOCKSIZE_> r, u, v, w, p, l;
	
	inline __m128 _heaviside(const __m128 phi, const __m128 invh, const __m128 one, const __m128 phalf, const __m128 mhalf) const;

	
	inline __m128 _getgamma(const __m128 phi, const __m128 inv_smoothlength, 
							const __m128 gamma1, const __m128 gamma2,
							const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const;
	
    inline __m128 _getPC(const __m128 phi, const __m128 inv_smoothlength, 
							const __m128 pc1, const __m128 pc2,
							const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const;
    
	void _sse_convert(const float * const gptfirst, const int gptfloats, float * const r, 
					  float * const u, float * const v, float * const w, float * const p, float * const l);
	
	void _sse_maxsos(const float * const r, const float * const  u, const float * const  v, 
					 const float * const w, const float * const  p, const float * const l, float * const sos);
	
public:
	
	MaxSpeedOfSound_SSE(const Real gamma1, const Real gamma2, const Real smoothlength, const Real pc1, const Real pc2):
		MaxSpeedOfSound_CPP(gamma1, gamma2, smoothlength, pc1, pc2) { }
	
	Real compute(const Real * const src, const int gptfloats)
	{
#if defined(_SP_COMP_)
		float __attribute__((aligned(_ALIGNBYTES_))) sos[4] = {0,0,0,0};
		
		const int slicesrcs = _BLOCKSIZE_*_BLOCKSIZE_;
		
		for(int islice=0; islice<_BLOCKSIZE_; islice++)
		{
			_sse_convert(src + islice*gptfloats*slicesrcs, gptfloats, &r.ref(0,0), &u.ref(0,0), &v.ref(0,0), &w.ref(0,0), &p.ref(0,0), &l.ref(0,0));
			
			_sse_maxsos(r.ptr(0,0), u.ptr(0,0), v.ptr(0,0), w.ptr(0,0), p.ptr(0,0), l.ptr(0,0), sos);
		}
		
		return  std::max( std::max(sos[0], sos[1]), std::max(sos[2], sos[3]) );
#else
		return MaxSpeedOfSound_CPP::compute(src, gptfloats); //C++ fallback for double precision
#endif
	}
};
#endif
