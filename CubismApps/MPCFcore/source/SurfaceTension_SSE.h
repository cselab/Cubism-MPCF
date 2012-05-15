/*
 *  SurfaceTension_SSE.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 3/9/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "common.h"
#include "SurfaceTension_CPP.h"
#include "DivTensor_SSE.h"

class SurfaceTension_SSE: public virtual DivTensor_SSE, public virtual SurfaceTension_CPP
{
protected:
		
	//used only for biphase
	inline __m128 _compute_ls(const __m128 G,  const __m128 G2, 
							  const __m128 ls_factor, const __m128 one)
	{
		const __m128 x = (G-G2)*ls_factor;
		const __m128 lambda = _mm_min_ps(_mm_max_ps(x ,_mm_setzero_ps()), one);
		
		return one-lambda;
	}
	
	template<bool biphase>
	inline void _convert(const float * const pt, const float ls_factor, float&u, float&v, float&w, float& ls)
	{
		if (biphase)
		{
			const float inv_rho = (float)1/pt[0];
			
			u = pt[1]*inv_rho;
			v = pt[2]*inv_rho;
			w = pt[3]*inv_rho;
			
			const float x = (pt[5]-G2)*ls_factor;
			const float lambda = min(max(x ,(float)0), (float)1);
			ls = 1 - lambda;
		}
		else
		{
			const float inv_rho = (float)1/pt[0];
			
			u = pt[1]*inv_rho;
			v = pt[2]*inv_rho;
			w = pt[3]*inv_rho;
			ls = 1; 
		}
	}
	
	template<bool biphase>
	void _convert_sse(const float * const gptfirst, const int gptfloats, const int rowgpts, 
					  InputSOAf_ST& _u, InputSOAf_ST& _v, InputSOAf_ST& _w, InputSOAf_ST& _mu)
	{
		const float ls_factor = ((float)1/(G1 - G2));
		
		const __m128 F_1 = _mm_set_ps1(1);
		const __m128 F_G2 = _mm_set_ps1(G2);		
		const __m128 F_LS = _mm_set_ps1(ls_factor);		
		
		const int stride = gptfloats*rowgpts;
		static const int PITCHOUT = _u.PITCH;
		
		float * const ubase = &_u.ref(-1,-1);
		float * const vbase = &_v.ref(-1,-1);
		float * const wbase = &_w.ref(-1,-1);
		float * const mubase = &_mu.ref(-1,-1);
		
		for(int iy=0; iy<_BLOCKSIZE_+2; iy++)
		{
			float * const u = ubase + iy*PITCHOUT;
			float * const v = vbase + iy*PITCHOUT;
			float * const w = wbase + iy*PITCHOUT;
			float * const mu = mubase + iy*PITCHOUT;
			
			const float * const leftghost = gptfirst + stride*iy;
			
			_convert<biphase>(leftghost, ls_factor, u[0], v[0], w[0], mu[0]);
			
			for(int ix=1; ix<_BLOCKSIZE_+1; ix+=4)
			{
				const float * const gpt0 = leftghost + gptfloats*(ix+0);
				const float * const gpt1 = leftghost + gptfloats*(ix+1);
				const float * const gpt2 = leftghost + gptfloats*(ix+2);
				const float * const gpt3 = leftghost + gptfloats*(ix+3);
				
				__m128 data0 = _mm_loadu_ps(gpt0);
				__m128 data1 = _mm_loadu_ps(gpt1);
				__m128 data2 = _mm_loadu_ps(gpt2);
				__m128 data3 = _mm_loadu_ps(gpt3);
				
				_MM_TRANSPOSE4_PS(data0, data1, data2, data3);
				
#ifdef _PREC_DIV_
				const __m128 inv_rho = F_1/data0;			
#else
				const __m128 inv_rho = better_rcp(data0);
#endif
				_mm_store_ps(u + ix, data1*inv_rho);
				_mm_store_ps(v + ix, data2*inv_rho);
				_mm_store_ps(w + ix, data3*inv_rho);
				
				if (biphase)
					_mm_store_ps(mu + ix, _compute_ls(_mm_set_ps(gpt3[5], gpt2[5], gpt1[5], gpt0[5]), F_G2, F_LS, F_1));
				else
					_mm_store_ps(mu + ix, _mm_setzero_ps());
			}
			
			_convert<biphase>(leftghost + gptfloats*(_BLOCKSIZE_+1), ls_factor, u[_BLOCKSIZE_+1], v[_BLOCKSIZE_+1], w[_BLOCKSIZE_+1], mu[_BLOCKSIZE_+1]);
		}
	}
	
	void _convert(const float * const gptfirst, const int gptfloats, const int rowgpts)
	{	
#ifndef _SP_COMP_
		printf("Diffusion_SSE::_convert: you should not be here in double precision. Aborting.\n");
		abort();
#else
		InputSOA_ST& u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &mu = ringls.ref();

		if (G1 == G2)
			_convert_sse<false>(gptfirst, gptfloats, rowgpts, u, v, w, mu);
		else
			_convert_sse<true>(gptfirst, gptfloats, rowgpts, u, v, w, mu);
#endif
	}
	
	void _tensor_xface(const TempSOAf_ST& nx0, const TempSOAf_ST& nx1, 
				  const TempSOAf_ST& ny0, const TempSOAf_ST& ny1, 
				  const TempSOAf_ST& nz0, const TempSOAf_ST& nz1)
	{
#ifndef _SP_COMP_
		printf("Diffusion_SSE::_tensor_xface: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const __m128 F_1_3 = _mm_set1_ps((Real)(1./3));
				
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
						
			for(int ix=0; ix<TempPiXSOA_ST::NX; ix+=4)
			{
				const __m128 nx = _mm_load_ps(txxptr + ix);
				const __m128 ny = _mm_load_ps(txyptr + ix);
				const __m128 nz = _mm_load_ps(txzptr + ix);
				
				const __m128 mag2 = nx*nx + ny*ny + nz*nz;
				const __m128 tmp = myrsqrt(mag2);
				const __m128 inv_mag = _mm_and_ps(tmp, _mm_cmpgt_ps(mag2, _mm_setzero_ps()));
				
				_mm_store_ps(txxptr + ix, (F_1_3*mag2 - nx*nx)*inv_mag);
				_mm_store_ps(txyptr + ix, _mm_setzero_ps() - nx*ny*inv_mag);
				_mm_store_ps(txzptr + ix, _mm_setzero_ps() - nx*nz*inv_mag);				
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
		const __m128 F_1_3 = _mm_set1_ps((Real)(1./3));
		
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
			
			for(int ix=0; ix<TempPiYSOA_ST::NX; ix+=4)
			{
				const __m128 nx = _mm_load_ps(tyxptr + ix);
				const __m128 ny = _mm_load_ps(tyyptr + ix);
				const __m128 nz = _mm_load_ps(tyzptr + ix);
				
				const __m128 mag2 = nx*nx + ny*ny + nz*nz;
				const __m128 tmp = myrsqrt(mag2);
				const __m128 inv_mag = _mm_and_ps(tmp, _mm_cmpgt_ps(mag2, _mm_setzero_ps()));
				
				_mm_store_ps(tyxptr + ix, _mm_setzero_ps() - ny*nx*inv_mag);
				_mm_store_ps(tyyptr + ix, (F_1_3*mag2 - ny*ny)*inv_mag);
				_mm_store_ps(tyzptr + ix, _mm_setzero_ps() - ny*nz*inv_mag);				
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
		const __m128 F_1_3 = _mm_set1_ps((Real)(1./3));
		
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
			
			for(int ix=0; ix<TempPiZSOA_ST::NX; ix+=4)
			{
				const __m128 nx = _mm_load_ps(tzxptr + ix);
				const __m128 ny = _mm_load_ps(tzyptr + ix);
				const __m128 nz = _mm_load_ps(tzzptr + ix);
				
				const __m128 mag2 = nx*nx + ny*ny + nz*nz;
				const __m128 tmp = myrsqrt(mag2);
				const __m128 inv_mag = _mm_and_ps(tmp, _mm_cmpgt_ps(mag2, _mm_setzero_ps()));
				
				_mm_store_ps(tzxptr + ix, _mm_setzero_ps() - nz*nx*inv_mag);
				_mm_store_ps(tzyptr + ix, _mm_setzero_ps() - nz*ny*inv_mag);
				_mm_store_ps(tzzptr + ix, (F_1_3*mag2 - nz*nz)*inv_mag);				
			}
		}
#endif
	}
	
public:
	
	SurfaceTension_SSE(const Real a = 1, const Real dtinvh = 1,
					   const Real G1 = 1/(2.5-1), const Real G2 = 1/(2.1-1), const Real h = 1, const Real sigma=1):
	SurfaceTension_CPP(a, dtinvh, G1, G2, h, sigma),
	DivTensor_SSE(a, dtinvh, h, sigma), DivTensor_CPP(a, dtinvh, h, sigma)
	{ 
	}
};
