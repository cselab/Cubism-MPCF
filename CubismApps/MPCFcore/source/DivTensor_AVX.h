/*
 *  DivTensor_AVX.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 3/12/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "common.h"

class DivTensor_AVX: public virtual DivTensor_SSE
{	
public:

 DivTensor_AVX(const Real a = 1, const Real dtinvh = 1, const Real h = 1, const Real sigma=1):
  DivTensor_SSE(a, dtinvh, h, sigma), DivTensor_CPP(a, dtinvh, h, sigma)
	{
	}
	
protected:
	
	inline __m256 _LEFT(const __m128 W, const __m256 C) const
	{
		const __m256 Cpermuted = _mm256_permute2f128_ps(C, C, 1);
		const __m256 WC0 = _mm256_insertf128_ps(Cpermuted, W, 0);
		return _mm256_shuffle_ps(_mm256_shuffle_ps(WC0, C, _MM_SHUFFLE(0,0,3,3)), C, _MM_SHUFFLE(2,1,3,0));
	}
	
	inline __m256 _LEFT(const __m256 W, const __m256 C) const
	{
		__m256 WEST = _mm256_insertf128_ps(_mm256_setzero_ps(), _mm256_extractf128_ps(W, 1), 0);
		WEST = _mm256_insertf128_ps(WEST, _mm256_extractf128_ps(C, 0), 1);
		return  _mm256_shuffle_ps(_mm256_shuffle_ps(WEST, C, _MM_SHUFFLE(0,0,3,3)), C, _MM_SHUFFLE(2,1,3,0));
	}
	
	inline __m256 _RIGHT(const __m256 C, const __m128 E) const
	{
		const __m256 EAST = _mm256_insertf128_ps(_mm256_permute2f128_ps(C, C, 1), E, 1);
		return _mm256_shuffle_ps(C, _mm256_shuffle_ps(C, EAST, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1));
	}
	
	inline __m256 _RIGHT(const __m256 C, const __m256 E) const
	{
		__m256 EAST = _mm256_insertf128_ps(_mm256_setzero_ps(), _mm256_extractf128_ps(C, 1), 0);
		EAST = _mm256_insertf128_ps(EAST, _mm256_extractf128_ps(E, 0), 1);
		return _mm256_shuffle_ps(C, _mm256_shuffle_ps(C, EAST, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1));
	}
	
	void _transpose4x2 (__m256& ax, __m256& ay, __m256& az, __m256& aw) const
	{
		const __m256 t0 = _mm256_unpacklo_ps(ax, ay);
		const __m256 t2 = _mm256_unpacklo_ps(az, aw);
		const __m256 t1 = _mm256_unpackhi_ps(ax, ay);
		const __m256 t3 = _mm256_unpackhi_ps(az, aw);
		
		ax = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1,0,1,0));
		ay = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3,2,3,2));
		az = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1,0,1,0));
		aw = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3,2,3,2));			 
	}
	
	//virtual method of DivTensor_CPP overridden here
    void _corners(const InputSOAf_ST& _ls0, const InputSOAf_ST& _ls1, TempSOAf_ST& _nx, TempSOAf_ST& _ny, TempSOAf_ST& _nz)
	{	
		static const int PITCHIN = _ls0.PITCH;
		const float * const ls0base = _ls0.ptr(0,0);
		const float * const ls1base = _ls1.ptr(0,0);
		
		static const int PITCHOUT = _nx.PITCH;
		float * const nxbase = &_nx.ref(0,0);
		float * const nybase = &_ny.ref(0,0);
		float * const nzbase = &_nz.ref(0,0);
		
		for(int iy=0; iy<TempSOA_ST::NY; iy++)
		{
			const float * const ls0 = ls0base + iy*PITCHIN;
			const float * const ls1 = ls1base + iy*PITCHIN;
			
			float * const nx = nxbase + iy*PITCHOUT;
			float * const ny = nybase + iy*PITCHOUT;
			float * const nz = nzbase + iy*PITCHOUT;
			
			__m128 WA0 = _mm_load_ps(ls0 - 4);
			__m128 WB0 = _mm_load_ps(ls0 - PITCHIN -4);
			__m128 WA1 = _mm_load_ps(ls1 - 4);
			__m128 WB1 = _mm_load_ps(ls1 - PITCHIN -4);
			
			for(int ix=0; ix<TempSOA_ST::NX; ix+=8)
			{	
				const __m256 c00m1 = _mm256_load_ps(ls0 + ix);
				const __m256 cm10m1 = _LEFT(WA0, c00m1);
				const __m256 c0m1m1 = _mm256_load_ps(ls0 + ix - PITCHIN);
				const __m256 cm1m1m1 = _LEFT(WB0, c0m1m1);
				
				const __m256 c000 = _mm256_load_ps(ls1 + ix);
				const __m256 cm100 = _LEFT(WA1, c000);
				const __m256 c0m10 = _mm256_load_ps(ls1 + ix - PITCHIN);
				const __m256 cm1m10 = _LEFT(WB1, c0m10);
				
				_mm256_store_ps(nx + ix, c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1);
				_mm256_store_ps(ny + ix, cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1  + c00m1-c0m1m1);
				_mm256_store_ps(nz + ix, cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1  + c000-c00m1);
				
				WA0 = _mm256_extractf128_ps(c00m1, 1);
				WB0 = _mm256_extractf128_ps(c0m1m1, 1);
				WA1 = _mm256_extractf128_ps(c000, 1);
				WB1 = _mm256_extractf128_ps(c0m10, 1);
			}
		}
	}
	
	void _copyback(float * const gptfirst, const int gptfloats, const int rowgpts)
	{	
#ifndef _SP_COMP_
		printf("DivTensor_AVX::_copyback: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const __m256 F = _mm256_set1_ps(sigma*dtinvh/(16*h));
		const __m128 A = _mm_set1_ps(a);
		const __m256 F2= F*_mm256_set1_ps((float)0.5);
		
		static const int PITCHIN = OutputSOA::PITCH;
		const int stride = gptfloats*rowgpts;
		
		const float * const baserhsu = rhsu.ptr(0,0);
		const float * const baserhsv = rhsv.ptr(0,0);
		const float * const baserhsw = rhsw.ptr(0,0);
		const float * const baserhss = rhss.ptr(0,0);
		
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{
			const float * const ptrrhsu = baserhsu + iy*PITCHIN;
			const float * const ptrrhsv = baserhsv + iy*PITCHIN;
			const float * const ptrrhsw = baserhsw + iy*PITCHIN;
			const float * const ptrrhss = baserhss + iy*PITCHIN;
			
			float * const basegp = 1 + gptfirst + iy*stride;
			
			for(int ix=0; ix<OutputSOA::NX; ix+=8)
			{
				float * const gp0 = basegp + gptfloats*(ix+0);
				float * const gp1 = basegp + gptfloats*(ix+1);
				float * const gp2 = basegp + gptfloats*(ix+2);
				float * const gp3 = basegp + gptfloats*(ix+3);
				float * const gp4 = basegp + gptfloats*(ix+4);
				float * const gp5 = basegp + gptfloats*(ix+5);
				float * const gp6 = basegp + gptfloats*(ix+6);
				float * const gp7 = basegp + gptfloats*(ix+7);
				
				__m256 rhs0 = F*_mm256_load_ps(ptrrhsu + ix);
				__m256 rhs1 = F*_mm256_load_ps(ptrrhsv + ix);
				__m256 rhs2 = F*_mm256_load_ps(ptrrhsw + ix);
				__m256 rhs3 = F2*_mm256_load_ps(ptrrhss + ix);
				
				_transpose4x2(rhs0, rhs1, rhs2, rhs3);
				
				_mm_storeu_ps(gp0, A*_mm_loadu_ps(gp0) + _mm256_extractf128_ps(rhs0, 0));
				_mm_storeu_ps(gp1, A*_mm_loadu_ps(gp1) + _mm256_extractf128_ps(rhs1, 0));
				_mm_storeu_ps(gp2, A*_mm_loadu_ps(gp2) + _mm256_extractf128_ps(rhs2, 0));
				_mm_storeu_ps(gp3, A*_mm_loadu_ps(gp3) + _mm256_extractf128_ps(rhs3, 0));
				_mm_storeu_ps(gp4, A*_mm_loadu_ps(gp4) + _mm256_extractf128_ps(rhs0, 1));
				_mm_storeu_ps(gp5, A*_mm_loadu_ps(gp5) + _mm256_extractf128_ps(rhs1, 1));
				_mm_storeu_ps(gp6, A*_mm_loadu_ps(gp6) + _mm256_extractf128_ps(rhs2, 1));
				_mm_storeu_ps(gp7, A*_mm_loadu_ps(gp7) + _mm256_extractf128_ps(rhs3, 1));
			}
		}
#endif
	}
	
	template<int PITCHIN0, int PITCHIN1, int PITCHOUT>
	void _div_dxy_avx(const float * const basesrc0, const float * const basesrc1, float * const basedest)
	{
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{
			const float * const src0 = basesrc0 + iy*PITCHIN0;
			const float * const src1 = basesrc1 + iy*PITCHIN1;
			float * const dest = basedest + iy*PITCHOUT;
			
			__m256 C = _mm256_load_ps(src0);
			
			for(int ix=0; ix<OutputSOA::NX; ix+=8)
			{
				const __m256 E = _mm256_load_ps(src0 + ix + 8);
				
				_mm256_store_ps(dest + ix,
								_RIGHT(C, E) - C + 
								_mm256_load_ps(src1 + ix + PITCHIN1) - _mm256_load_ps(src1 + ix));
				
				C = E;
			}
		}
	}
	
	template<int PITCHIN, int PITCHOUT>
	void _div_dz_avx(const float * const basesrc0, const float * const basesrc1, float * const basedest)
	{
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{	
			const float * const src0 = basesrc0 + iy*PITCHIN;
			const float * const src1 = basesrc1 + iy*PITCHIN;
			float * const dest = basedest + iy*PITCHOUT;
			
			//The following code is re-written to avoid unintended write-to-load forwarding stalls
			//for(int ix=0; ix<OutputSOA::NX; ix+=8)
			//	_mm256_store_ps(dest + ix, _mm256_load_ps(dest + ix) + _mm256_load_ps(src0 + ix) - _mm256_load_ps(src1 + ix)); 
			//VTUNE told me that.
			
			float __attribute__((aligned(_ALIGNBYTES_))) tmp[PITCHOUT];
			
			for(int ix=0; ix<OutputSOA::NX; ix+=8)
				_mm256_store_ps(tmp + ix, _mm256_load_ps(dest + ix) + _mm256_load_ps(src0 + ix) - _mm256_load_ps(src1 + ix)); 
			
			for(int ix=0; ix<OutputSOA::NX; ix+=8)
				_mm256_store_ps(dest + ix, _mm256_load_ps(tmp + ix)); 
		}
	}
	
#ifdef _SP_COMP_
	void _div_dxy()
	{	
		static const int PIN0 = TempPiXSOA_ST::PITCH;
		static const int PIN1 = TempPiYSOA_ST::PITCH;
		static const int POUT = OutputSOA::PITCH;
		
		_div_dxy_avx<PIN0, PIN1, POUT>(txx.ptr(0,0), tyx.ptr(0,0), &rhsu.ref(0,0));
		_div_dxy_avx<PIN0, PIN1, POUT>(txy.ptr(0,0), tyy.ptr(0,0), &rhsv.ref(0,0));
		_div_dxy_avx<PIN0, PIN1, POUT>(txz.ptr(0,0), tyz.ptr(0,0), &rhsw.ref(0,0));
		_div_dxy_avx<PIN0, PIN1, POUT>(utx.ptr(0,0), uty.ptr(0,0), &rhss.ref(0,0));
	}
	
	void _div_dz(const TempPiZSOAf_ST& tzx0, const TempPiZSOAf_ST& tzy0, const TempPiZSOAf_ST& tzz0, const TempPiZSOAf_ST& utz0,
				 const TempPiZSOAf_ST& tzx1, const TempPiZSOAf_ST& tzy1, const TempPiZSOAf_ST& tzz1, const TempPiZSOAf_ST& utz1)
	{
		static const int PIN0 = TempPiZSOA_ST::PITCH;
		static const int POUT = OutputSOA::PITCH;
		
		_div_dz_avx<PIN0, POUT>(tzx1.ptr(0,0), tzx0.ptr(0,0), &rhsu.ref(0,0));
		_div_dz_avx<PIN0, POUT>(tzy1.ptr(0,0), tzy0.ptr(0,0), &rhsv.ref(0,0));
		_div_dz_avx<PIN0, POUT>(tzz1.ptr(0,0), tzz0.ptr(0,0), &rhsw.ref(0,0));
		_div_dz_avx<PIN0, POUT>(utz1.ptr(0,0), utz0.ptr(0,0), &rhss.ref(0,0));
	}
#endif
	
	void _udot_tx(const InputSOAf_ST& u, const InputSOAf_ST& v, const InputSOAf_ST& w)
	{
#ifndef _SP_COMP_
		printf("DivTensor_AVX::_udot_tx: you should not be here in double precision. Aborting.\n");
		abort();
#else
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHTENSOR = TempPiXSOA_ST::PITCH;
		
		const float * const ubase = u.ptr(0,0);
		const float * const vbase = v.ptr(0,0);
		const float * const wbase = w.ptr(0,0);
		
		const float * const txxbase = txx.ptr(0,0);
		const float * const txybase = txy.ptr(0,0);
		const float * const txzbase = txz.ptr(0,0);
		
		float * const utxbase = &utx.ref(0,0);
		
		for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
		{
			const float * const uptr = ubase + iy*PITCHIN;
			const float * const vptr = vbase + iy*PITCHIN;
			const float * const wptr = wbase + iy*PITCHIN;
			
			const float * const txxptr = txxbase + iy*PITCHTENSOR;
			const float * const txyptr = txybase + iy*PITCHTENSOR;
			const float * const txzptr = txzbase + iy*PITCHTENSOR;
			
			float * const utxptr = utxbase + iy*PITCHTENSOR;
			
			__m256 WU = _mm256_load_ps(uptr - 8);
			__m256 WV = _mm256_load_ps(vptr - 8);
			__m256 WW = _mm256_load_ps(wptr - 8);
			
			for(int ix=0; ix<TempPiXSOA_ST::NX; ix += 8)
			{
				const __m256 CU = _mm256_load_ps(uptr + ix);
				const __m256 CV = _mm256_load_ps(vptr + ix);
				const __m256 CW = _mm256_load_ps(wptr + ix);
				
				_mm256_store_ps(utxptr + ix , 
								_mm256_load_ps(txxptr + ix)*(CU + _LEFT(WU, CU)) +
								_mm256_load_ps(txyptr + ix)*(CV + _LEFT(WV, CV)) +
								_mm256_load_ps(txzptr + ix)*(CW + _LEFT(WW, CW)));
				
				WU = CU;
				WV = CV;
				WW = CW;
			}
		}
#endif
	}
	
	void _udot_ty(const InputSOAf_ST& u, const InputSOAf_ST& v, const InputSOAf_ST& w)
	{
#ifndef _SP_COMP_
		printf("DivTensor_SSE::_udot_ty: you should not be here in double precision. Aborting.\n");
		abort();
#else
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHTENSOR = TempPiYSOA_ST::PITCH;
		
		const float * const ubase = u.ptr(0,0);
		const float * const vbase = v.ptr(0,0);
		const float * const wbase = w.ptr(0,0);
		
		const float * const tyxbase = tyx.ptr(0,0);
		const float * const tyybase = tyy.ptr(0,0);
		const float * const tyzbase = tyz.ptr(0,0);
		
		float * const utybase = &uty.ref(0,0);
		
		for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		{
			const float * const uptr = ubase + iy*PITCHIN;
			const float * const vptr = vbase + iy*PITCHIN;
			const float * const wptr = wbase + iy*PITCHIN;
			
			const float * const tyxptr = tyxbase + iy*PITCHTENSOR;
			const float * const tyyptr = tyybase + iy*PITCHTENSOR;
			const float * const tyzptr = tyzbase + iy*PITCHTENSOR;
			
			float * const utyptr = utybase + iy*PITCHTENSOR;
			
			for(int ix=0; ix<TempPiYSOA_ST::NX; ix +=  8)
			{
				_mm256_store_ps(utyptr + ix , 
								_mm256_load_ps(tyxptr + ix)*(_mm256_load_ps(uptr + ix) + _mm256_load_ps(uptr + ix - PITCHIN)) +
								_mm256_load_ps(tyyptr + ix)*(_mm256_load_ps(vptr + ix) + _mm256_load_ps(vptr + ix - PITCHIN)) +
								_mm256_load_ps(tyzptr + ix)*(_mm256_load_ps(wptr + ix) + _mm256_load_ps(wptr + ix - PITCHIN)));
			}
		}
#endif
	}
	
	void _udot_tz(const InputSOAf_ST& u0, const InputSOAf_ST& v0, const InputSOAf_ST& w0,
				  const InputSOAf_ST& u1, const InputSOAf_ST& v1, const InputSOAf_ST& w1,
				  const TempPiZSOAf_ST& tzx, const TempPiZSOAf_ST& tzy, const TempPiZSOAf_ST& tzz,
				  TempPiZSOAf_ST& utz)
	{
#ifndef _SP_COMP_
		printf("DivTensor_SSE::_udot_tz: you should not be here in double precision. Aborting.\n");
		abort();
#else
		static const int PITCHIN = InputSOA_ST::PITCH;
		static const int PITCHTENSOR = TempPiZSOA_ST::PITCH;
		
		const float * const ubase0 = u0.ptr(0,0);
		const float * const vbase0 = v0.ptr(0,0);
		const float * const wbase0 = w0.ptr(0,0);
		const float * const ubase1 = u1.ptr(0,0);
		const float * const vbase1 = v1.ptr(0,0);
		const float * const wbase1 = w1.ptr(0,0);
		
		const float * const tzxbase = tzx.ptr(0,0);
		const float * const tzybase = tzy.ptr(0,0);
		const float * const tzzbase = tzz.ptr(0,0);
		
		float * const utzbase = &utz.ref(0,0);
		
		for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		{
			const float * const uptr0 = ubase0 + iy*PITCHIN;
			const float * const vptr0 = vbase0 + iy*PITCHIN;
			const float * const wptr0 = wbase0 + iy*PITCHIN;
			const float * const uptr1 = ubase1 + iy*PITCHIN;
			const float * const vptr1 = vbase1 + iy*PITCHIN;
			const float * const wptr1 = wbase1 + iy*PITCHIN;
			
			const float * const tzxptr = tzxbase + iy*PITCHTENSOR;
			const float * const tzyptr = tzybase + iy*PITCHTENSOR;
			const float * const tzzptr = tzzbase + iy*PITCHTENSOR;
			
			float * const utzptr = utzbase + iy*PITCHTENSOR;
			
			for(int ix=0; ix<TempPiZSOA_ST::NX; ix += 8)
			{
				_mm256_store_ps(utzptr + ix , 
								_mm256_load_ps(tzxptr + ix)*(_mm256_load_ps(uptr0 + ix) + _mm256_load_ps(uptr1 + ix)) +
								_mm256_load_ps(tzyptr + ix)*(_mm256_load_ps(vptr0 + ix) + _mm256_load_ps(vptr1 + ix)) +
								_mm256_load_ps(tzzptr + ix)*(_mm256_load_ps(wptr0 + ix) + _mm256_load_ps(wptr1 + ix)));
			}
		}
#endif
	}
	
	template <bool accum> inline void write256(float * const dest, const __m256 src)
	{
		if (accum)
			_mm256_store_ps(dest, _mm256_load_ps(dest) + src);
		else
			_mm256_store_ps(dest, src);
	}
	
	/* The method _average_xface(.)
	 is left unoveridden as it was observed to be slower than the base implementation (ILP-issues).
	 */
	
	template<bool accum>
	void _average_yface(const float factor, const TempSOA_ST& src0, const TempSOA_ST& src1, TempPiYSOA_ST& dest)
	{
		static const int PITCHIN = TempSOA_ST::PITCH;
		static const int PITCHOUT = TempPiYSOA_ST::PITCH;
		
		const float * const src0base = src0.ptr(0,0); 
		const float * const src1base = src1.ptr(0,0); 
		float * const destbase = &dest.ref(0,0); 
		
		const __m256 F = _mm256_set1_ps(factor);
		
		for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		{
			const float * const src0ptr = src0base + iy*PITCHIN; 
			const float * const src1ptr = src1base + iy*PITCHIN; 
			float * const destptr = destbase + iy*PITCHOUT;
			
			for(int ix=0; ix<TempPiYSOA_ST::NX; ix+=8)
			{
				const __m256 C0 = _mm256_load_ps(src0ptr + ix);
				const __m256 C1 = _mm256_load_ps(src1ptr + ix);
				const __m128 E0 = _mm_load_ps(src0ptr + ix + 8);
				const __m128 E1 = _mm_load_ps(src1ptr + ix + 8);
				
				write256<accum>(destptr + ix, F*(C0 + _RIGHT(C0, E0) + C1 + _RIGHT(C1, E1)));
			}
		}
	}
	
	template<bool accum>
	void _average_zface(const float factor, const TempSOA_ST& src, TempPiZSOA_ST& dest)
	{
		static const int PITCHIN = TempSOA_ST::PITCH;
		static const int PITCHOUT = TempPiZSOA_ST::PITCH;
		
		const float * const srcbase = src.ptr(0,0); 
		float * const destbase = &dest.ref(0,0); 
		
		const __m256 F = _mm256_set1_ps(factor);
		
		for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		{
			const float * const srcptr = srcbase + iy*PITCHIN; 
			float * const destptr = destbase + iy*PITCHOUT;
			
			for(int ix=0; ix<TempPiZSOA_ST::NX; ix+=8)
			{
				const __m256 C0 = _mm256_load_ps(srcptr + ix);
				const __m256 C1 = _mm256_load_ps(srcptr + PITCHIN + ix);
				const __m128 E0 = _mm_load_ps(srcptr + ix + 8);
				const __m128 E1 = _mm_load_ps(srcptr + ix + PITCHIN + 8);
				
				write256<accum>(destptr + ix, F*(C0 +  _RIGHT(C0, E0) + C1 + _RIGHT(C1, E1)));
			}
		}
	}
};
