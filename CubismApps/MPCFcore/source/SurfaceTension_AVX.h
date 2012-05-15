/*
 *  SurfaceTension_AVX.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 3/12/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "common.h"
#include "SurfaceTension_SSE.h"
#include "DivTensor_AVX.h"

class SurfaceTension_AVX: public virtual SurfaceTension_SSE, public virtual DivTensor_AVX
{	
	void _tensor_xface(const TempSOAf_ST& nx0, const TempSOAf_ST& nx1, 
					   const TempSOAf_ST& ny0, const TempSOAf_ST& ny1, 
					   const TempSOAf_ST& nz0, const TempSOAf_ST& nz1)
	{
#ifndef _SP_COMP_
		printf("Diffusion_SSE::_tensor_xface: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const __m256 F_1_3 = _mm256_set1_ps((Real)(1./3));
		
		_average_xface<false>(1, nx0, nx1, txx);
		_average_xface<false>(1, ny0, ny1, txy);
		_average_xface<false>(1, nz0, nz1, txz);
		
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiXSOA_ST::PITCH;
		
		float * const txxbase = &txx.ref(0,0); 
		float * const txybase = &txy.ref(0,0); 
		float * const txzbase = &txz.ref(0,0); 
		
		for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
		{
			float * const txxptr = txxbase + iy*PITCHOUT;
			float * const txyptr = txybase + iy*PITCHOUT;
			float * const txzptr = txzbase + iy*PITCHOUT;
			
			for(int ix=0; ix<TempPiXSOA_ST::NX; ix+=8)
			{
				const __m256 nx = _mm256_load_ps(txxptr + ix);
				const __m256 ny = _mm256_load_ps(txyptr + ix);
				const __m256 nz = _mm256_load_ps(txzptr + ix);
				
				const __m256 mag2 = nx*nx + ny*ny + nz*nz;
				const __m256 tmp = myrsqrt(mag2);
				const __m256 inv_mag = _mm256_and_ps(tmp, _mm256_cmp_ps(mag2, _mm256_setzero_ps(), _CMP_GT_OS));
				
				_mm256_store_ps(txxptr + ix, (F_1_3*mag2 - nx*nx)*inv_mag);
				_mm256_store_ps(txyptr + ix, _mm256_setzero_ps() - nx*ny*inv_mag);
				_mm256_store_ps(txzptr + ix, _mm256_setzero_ps() - nx*nz*inv_mag);				
			}
		}
#endif
	}
	
	void _tensor_yface(const TempSOAf_ST& nx0, const TempSOAf_ST& nx1, 
					   const TempSOAf_ST& ny0, const TempSOAf_ST& ny1, 
					   const TempSOAf_ST& nz0, const TempSOAf_ST& nz1)
	{
#ifndef _SP_COMP_
		printf("Diffusion_SSE::_tensor_yface: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const __m256 F_1_3 = _mm256_set1_ps((Real)(1./3));
		
		_average_yface<false>(1, nx0, nx1, tyx);
		_average_yface<false>(1, ny0, ny1, tyy);
		_average_yface<false>(1, nz0, nz1, tyz);
		
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiYSOA_ST::PITCH;
		
		float * const tyxbase = &tyx.ref(0,0); 
		float * const tyybase = &tyy.ref(0,0); 
		float * const tyzbase = &tyz.ref(0,0); 
		
		for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		{
			float * const tyxptr = tyxbase + iy*PITCHOUT;
			float * const tyyptr = tyybase + iy*PITCHOUT;
			float * const tyzptr = tyzbase + iy*PITCHOUT;
			
			for(int ix=0; ix<TempPiYSOA_ST::NX; ix+=8)
			{
				const __m256 nx = _mm256_load_ps(tyxptr + ix);
				const __m256 ny = _mm256_load_ps(tyyptr + ix);
				const __m256 nz = _mm256_load_ps(tyzptr + ix);
				
				const __m256 mag2 = nx*nx + ny*ny + nz*nz;
				const __m256 tmp = myrsqrt(mag2);
				const __m256 inv_mag = _mm256_and_ps(tmp, _mm256_cmp_ps(mag2, _mm256_setzero_ps(), _CMP_GT_OS));
				
				_mm256_store_ps(tyxptr + ix, _mm256_setzero_ps() - ny*nx*inv_mag);
				_mm256_store_ps(tyyptr + ix, (F_1_3*mag2 - ny*ny)*inv_mag);
				_mm256_store_ps(tyzptr + ix, _mm256_setzero_ps() - ny*nz*inv_mag);				
			}
		}
#endif
	}
	
	void _tensor_zface(const TempSOAf_ST& nx, const TempSOAf_ST& ny, const TempSOAf_ST& nz,
					   TempPiZSOAf_ST& tzx, TempPiZSOAf_ST& tzy, TempPiZSOAf_ST& tzz)
	{
#ifndef _SP_COMP_
		printf("Diffusion_SSE::_tensor_zface: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const __m256 F_1_3 = _mm256_set1_ps((Real)(1./3));
		
		_average_zface<false>(1, nx, tzx);
		_average_zface<false>(1, ny, tzy);
		_average_zface<false>(1, nz, tzz);
		
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHOUT = TempPiZSOA_ST::PITCH;
		
		float * const tzxbase = &tzx.ref(0,0); 
		float * const tzybase = &tzy.ref(0,0); 
		float * const tzzbase = &tzz.ref(0,0); 
		
		for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		{
			float * const tzxptr = tzxbase + iy*PITCHOUT;
			float * const tzyptr = tzybase + iy*PITCHOUT;
			float * const tzzptr = tzzbase + iy*PITCHOUT;
			
			for(int ix=0; ix<TempPiZSOA_ST::NX; ix+=8)
			{
				const __m256 nx = _mm256_load_ps(tzxptr + ix);
				const __m256 ny = _mm256_load_ps(tzyptr + ix);
				const __m256 nz = _mm256_load_ps(tzzptr + ix);
				
				const __m256 mag2 = nx*nx + ny*ny + nz*nz;
				const __m256 tmp = myrsqrt(mag2);
				const __m256 inv_mag = _mm256_and_ps(tmp, _mm256_cmp_ps(mag2, _mm256_setzero_ps(), _CMP_GT_OS));
				
				_mm256_store_ps(tzxptr + ix, _mm256_setzero_ps() - nz*nx*inv_mag);
				_mm256_store_ps(tzyptr + ix, _mm256_setzero_ps() - nz*ny*inv_mag);
				_mm256_store_ps(tzzptr + ix, (F_1_3*mag2 - nz*nz)*inv_mag);				
			}
		}
#endif
	}
	
public:
	
	SurfaceTension_AVX(const Real a = 1, const Real dtinvh = 1, const Real G1 = 1/(2.5-1), const Real G2 = 1/(2.1-1), const Real h = 1, const Real smoothing_length = 1, const Real sigma=1):
	SurfaceTension_SSE(a, dtinvh, G1, G2, h, smoothing_length, sigma),
	     SurfaceTension_CPP(a, dtinvh, G1, G2, h, smoothing_length, sigma),
	     DivTensor_AVX(a, dtinvh, h, sigma),
	     DivTensor_SSE(a, dtinvh, h, sigma),
	     DivTensor_CPP(a, dtinvh, h, sigma)
	{ 
	}
};
