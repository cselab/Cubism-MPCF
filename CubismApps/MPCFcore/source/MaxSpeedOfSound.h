/*
 *  MaxSpeedOfSound.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include <xmmintrin.h>
#include <algorithm>

#include "FlowStep_CPP.h"

class MaxSpeedOfSound_CPP
{
protected:
	
	const Real m_gamma1, m_gamma2, m_smoothlength;
    Real m_pc1, m_pc2;
	
	inline Real _getgamma(const Real phi);
    inline Real _getPC(const Real phi);
	
public:
	
	MaxSpeedOfSound_CPP(const Real gamma1=2.5, const Real gamma2=2.1, const Real smoothlength=1, const Real pc1=0, const Real pc2=0):
		m_gamma1(gamma1), m_gamma2(gamma2), m_smoothlength(smoothlength), m_pc1(pc1), m_pc2(pc2) { }
	
	Real compute(const Real * const src, const int gptfloats);
	
	Real gamma1() const { return m_gamma1; }
	Real gamma2() const { return m_gamma2; }
	Real smoothlength() const { return m_smoothlength; }
	
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
		cout << "\tSOS TOTAL GFLOPS " << GFLOPSOS << endl;
		printf("\tSOS ASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tSOS RIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tSOS THIS ONE IS %.2f GFLOP/s,\t\"PER SOS\" %.2f FLOP/B\n", GFLOPSOS/MEASUREDTIME, OISOS);
		printf("\tSOS TIME PER BLOCK: %.5f ms (expected %.5f ms)\n",  1e3*MEASUREDTIME/NBLOCKS, 1e3*TSOS/NBLOCKS);
		cout << "\tSOS Expected Performance is: " << GFLOPSOS/(1.e9*GFLOPSOS/EPERFSOSAI) << " GFLOP/s [AI], " << GFLOPSOS/TSOS << " GFLOP/s [OI]"  << endl;
		printf("\tSOS EFFICIENCY: %.2f%% [AI] - %.2f%% [OI], HW-UTILIZATION: %.2f%%\n", 100.*(1.e9*GFLOPSOS/EPERFSOSAI)/MEASUREDTIME, 100.*TSOS/MEASUREDTIME, 100*(GFLOPSOS/MEASUREDTIME*1e9)/PEAKPERF);
		if (bAwk)
		{
			string kernelname = "SOS";
			
			awkMCorePredictions(kernelname, OISOS, EPERFSOS);
			awkMCore(kernelname, GFLOPSOS, PEAKPERF_CORE, PEAKPERF, PEAKBAND, MEASUREDTIME, TSOS, NBLOCKS, NT);
		}
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
	
	MaxSpeedOfSound_SSE(const Real gamma1=2.5, const Real gamma2=2.1, const Real smoothlength=1, const Real pc1=0, const Real pc2=0):
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
		return MaxSpeedOfSound_CPP::compute(src, gptfloats);//C++ fallback for double precision
#endif
	}
};
#endif

