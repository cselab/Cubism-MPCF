/*
 *  DivTensor_SSE.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 3/2/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <xmmintrin.h>
#include "DivTensor_CPP.h"

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

class DivTensor_SSE: public virtual DivTensor_CPP
{	
#define	LEFT(W,C) \
_mm_shuffle_ps(_mm_shuffle_ps(W,C, _MM_SHUFFLE(0,0,3,3)), C, _MM_SHUFFLE(2,1,3,0))
	
#define RIGHT(C, E)	\
_mm_shuffle_ps(C, _mm_shuffle_ps(C,E, _MM_SHUFFLE(0,0,3,3)), _MM_SHUFFLE(3,0,2,1))
	
protected:
	
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
			
			for(int ix=0; ix<TempSOA_ST::NX; ix+=4)
			{	
				const __m128 c00m1 = _mm_load_ps(ls0 + ix);
				const __m128 cm10m1 = LEFT(WA0, c00m1);
				const __m128 c0m1m1 = _mm_load_ps(ls0 + ix - PITCHIN);
				const __m128 cm1m1m1 = LEFT(WB0, c0m1m1);
				
				const __m128 c000 = _mm_load_ps(ls1 + ix);
				const __m128 cm100 = LEFT(WA1, c000);
				const __m128 c0m10 = _mm_load_ps(ls1 + ix - PITCHIN);
				const __m128 cm1m10 = LEFT(WB1, c0m10);
				
				_mm_store_ps(nx + ix, c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1);
				_mm_store_ps(ny + ix, cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1  + c00m1-c0m1m1);
				_mm_store_ps(nz + ix, cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1  + c000-c00m1);
				
				WA0 = c00m1;
				WB0 = c0m1m1;
				WA1 = c000;
				WB1 = c0m10;
			}
		}
	}
	
	void _copyback(float * const gptfirst, const int gptfloats, const int rowgpts)
	{	
#ifndef _SP_COMP_
		printf("DivTensor_SSE::_copyback: you should not be here in double precision. Aborting.\n");
		abort();
#else
		const __m128 F = _mm_set1_ps(sigma*dtinvh/(16*h));
		const __m128 A = _mm_set1_ps(a);
		const __m128 F2= F*_mm_set1_ps((float)0.5);
		
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
			
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				float * const gp0 = basegp + gptfloats*(ix+0);
				float * const gp1 = basegp + gptfloats*(ix+1);
				float * const gp2 = basegp + gptfloats*(ix+2);
				float * const gp3 = basegp + gptfloats*(ix+3);
				
				__m128 rhs0 = F*_mm_load_ps(ptrrhsu + ix);
				__m128 rhs1 = F*_mm_load_ps(ptrrhsv + ix);
				__m128 rhs2 = F*_mm_load_ps(ptrrhsw + ix);
				__m128 rhs3 = F2*_mm_load_ps(ptrrhss + ix);
				
				_MM_TRANSPOSE4_PS(rhs0, rhs1, rhs2, rhs3);
				
				_mm_storeu_ps(gp0, A*_mm_loadu_ps(gp0) + rhs0);
				_mm_storeu_ps(gp1, A*_mm_loadu_ps(gp1) + rhs1);
				_mm_storeu_ps(gp2, A*_mm_loadu_ps(gp2) + rhs2);
				_mm_storeu_ps(gp3, A*_mm_loadu_ps(gp3) + rhs3);
			}
		}
#endif
	}
	
	template<int PITCHIN0, int PITCHIN1, int PITCHOUT>
	void _div_dxy_sse(const float * const basesrc0, const float * const basesrc1, float * const basedest)
	{
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{
			const float * const src0 = basesrc0 + iy*PITCHIN0;
			const float * const src1 = basesrc1 + iy*PITCHIN1;
			float * const dest = basedest + iy*PITCHOUT;
			
			__m128 C = _mm_load_ps(src0);
			
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
			{
				const __m128 E = _mm_load_ps(src0 + ix + 4);
				
				_mm_store_ps(dest + ix,
							 RIGHT(C, E) - C + 
							 _mm_load_ps(src1 + ix + PITCHIN1) - _mm_load_ps(src1 + ix));
				
				C = E;
			}
		}
	}
	
	template<int PITCHIN, int PITCHOUT>
	void _div_dz_sse(const float * const basesrc0, const float * const basesrc1, float * const basedest)
	{
		for(int iy=0; iy<OutputSOA::NY; iy++)
		{
			const float * const src0 = basesrc0 + iy*PITCHIN;
			const float * const src1 = basesrc1 + iy*PITCHIN;
			float * const dest = basedest + iy*PITCHOUT;
			
			for(int ix=0; ix<OutputSOA::NX; ix+=4)
				_mm_store_ps(dest + ix, _mm_load_ps(dest + ix) + _mm_load_ps(src0 + ix) - _mm_load_ps(src1 + ix)); 
		}
	}
	
#ifdef _SP_COMP_
	void _div_dxy()
	{	
		static const int PIN0 = TempPiXSOA_ST::PITCH;
		static const int PIN1 = TempPiYSOA_ST::PITCH;
		static const int POUT = OutputSOA::PITCH;
		
		_div_dxy_sse<PIN0, PIN1, POUT>(txx.ptr(0,0), tyx.ptr(0,0), &rhsu.ref(0,0));
		_div_dxy_sse<PIN0, PIN1, POUT>(txy.ptr(0,0), tyy.ptr(0,0), &rhsv.ref(0,0));
		_div_dxy_sse<PIN0, PIN1, POUT>(txz.ptr(0,0), tyz.ptr(0,0), &rhsw.ref(0,0));
		_div_dxy_sse<PIN0, PIN1, POUT>(utx.ptr(0,0), uty.ptr(0,0), &rhss.ref(0,0));
	}
	
	void _div_dz(const TempPiZSOAf_ST& tzx0, const TempPiZSOAf_ST& tzy0, const TempPiZSOAf_ST& tzz0, const TempPiZSOAf_ST& utz0,
				 const TempPiZSOAf_ST& tzx1, const TempPiZSOAf_ST& tzy1, const TempPiZSOAf_ST& tzz1, const TempPiZSOAf_ST& utz1)
	{
		static const int PIN0 = TempPiZSOA_ST::PITCH;
		static const int POUT = OutputSOA::PITCH;
		
		_div_dz_sse<PIN0, POUT>(tzx1.ptr(0,0), tzx0.ptr(0,0), &rhsu.ref(0,0));
		_div_dz_sse<PIN0, POUT>(tzy1.ptr(0,0), tzy0.ptr(0,0), &rhsv.ref(0,0));
		_div_dz_sse<PIN0, POUT>(tzz1.ptr(0,0), tzz0.ptr(0,0), &rhsw.ref(0,0));
		_div_dz_sse<PIN0, POUT>(utz1.ptr(0,0), utz0.ptr(0,0), &rhss.ref(0,0));
	}
#endif
	
	void _udot_tx(const InputSOAf_ST& u, const InputSOAf_ST& v, const InputSOAf_ST& w)
	{
#ifndef _SP_COMP_
		printf("DivTensor_SSE::_udot_tx: you should not be here in double precision. Aborting.\n");
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
			
			for(int ix=0; ix<TempPiXSOA_ST::NX; ix += 4)
			{
				const __m128 C0 = _mm_load_ps(uptr + ix);
				const __m128 C1 = _mm_load_ps(vptr + ix);
				const __m128 C2 = _mm_load_ps(wptr + ix);
				
				const __m128 W0 = _mm_load_ps(uptr + ix - 4);
				const __m128 W1 = _mm_load_ps(vptr + ix - 4);
				const __m128 W2 = _mm_load_ps(wptr + ix - 4);
				
				_mm_store_ps(utxptr + ix , 
							 _mm_load_ps(txxptr + ix)*(C0 + LEFT(W0, C0)) +
							 _mm_load_ps(txyptr + ix)*(C1 + LEFT(W1, C1)) +
							 _mm_load_ps(txzptr + ix)*(C2 + LEFT(W2, C2)));
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
			
			for(int ix=0; ix<TempPiYSOA_ST::NX; ix += 4)
			{
				_mm_store_ps(utyptr + ix , 
							 _mm_load_ps(tyxptr + ix)*(_mm_load_ps(uptr + ix) + _mm_load_ps(uptr + ix - PITCHIN)) +
							 _mm_load_ps(tyyptr + ix)*(_mm_load_ps(vptr + ix) + _mm_load_ps(vptr + ix - PITCHIN)) +
							 _mm_load_ps(tyzptr + ix)*(_mm_load_ps(wptr + ix) + _mm_load_ps(wptr + ix - PITCHIN)));
			}
		}
#endif
	}
	
	void _udot_tz(const InputSOAf_ST& u0, const InputSOAf_ST& v0, const InputSOAf_ST& w0,
				  const InputSOAf_ST& u1, const InputSOAf_ST& v1, const InputSOAf_ST& w1,
				  const TempPiZSOAf_ST& tzx, const TempPiZSOAf_ST& tzy, const TempPiZSOAf_ST& tzz,
				  TempPiZSOAf_ST& utz)
	{
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
			
			for(int ix=0; ix<TempPiZSOA_ST::NX; ix += 4)
			{
				_mm_store_ps(utzptr + ix , 
							 _mm_load_ps(tzxptr + ix)*(_mm_load_ps(uptr0 + ix) + _mm_load_ps(uptr1 + ix)) +
							 _mm_load_ps(tzyptr + ix)*(_mm_load_ps(vptr0 + ix) + _mm_load_ps(vptr1 + ix)) +
							 _mm_load_ps(tzzptr + ix)*(_mm_load_ps(wptr0 + ix) + _mm_load_ps(wptr1 + ix)));
			}
		}
	}
	
protected:
	
	template <bool accum> inline void write128(float * const dest, const __m128 src)
	{
		if (accum)
			_mm_store_ps(dest, _mm_load_ps(dest) + src);
		else
			_mm_store_ps(dest, src);
	}
	
	template<bool accum>
	void _average_xface(const float factor, const TempSOA_ST& src0, const TempSOA_ST& src1, TempPiXSOA_ST& dest)
	{
		static const int PITCHIN = TempSOA_ST::PITCH;
		static const int PITCHOUT = TempPiXSOA_ST::PITCH;
		
		const float * const src0base = src0.ptr(0,0); 
		const float * const src1base = src1.ptr(0,0); 
		float * const destbase = &dest.ref(0,0); 
		
		const __m128 F = _mm_set1_ps(factor);
		
		for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
		{
			const float * const src0ptr = src0base + iy*PITCHIN; 
			const float * const src1ptr = src1base + iy*PITCHIN; 
			float * const destptr = destbase + iy*PITCHOUT;
			
			for(int ix=0; ix<TempPiXSOA_ST::NX; ix+=4)
				write128<accum>(destptr + ix, F*(_mm_load_ps(src0ptr + ix) + _mm_load_ps(src0ptr + ix + PITCHIN) + 
												 _mm_load_ps(src1ptr + ix) + _mm_load_ps(src1ptr + ix + PITCHIN)));
		}
	}
	
	template<bool accum>
	void _average_yface(const float factor, const TempSOA_ST& src0, const TempSOA_ST& src1, TempPiYSOA_ST& dest)
	{
		static const int PITCHIN = TempSOA_ST::PITCH;
		static const int PITCHOUT = TempPiYSOA_ST::PITCH;
		
		const float * const src0base = src0.ptr(0,0); 
		const float * const src1base = src1.ptr(0,0); 
		float * const destbase = &dest.ref(0,0); 
		
		const __m128 F = _mm_set1_ps(factor);
		
		for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
		{
			const float * const src0ptr = src0base + iy*PITCHIN; 
			const float * const src1ptr = src1base + iy*PITCHIN; 
			float * const destptr = destbase + iy*PITCHOUT;
			
			__m128 C0 = _mm_load_ps(src0ptr);
			__m128 C1 = _mm_load_ps(src1ptr);
			
			for(int ix=0; ix<TempPiYSOA_ST::NX; ix+=4)
			{
				const __m128 E0 = _mm_load_ps(src0ptr + ix + 4);
				const __m128 E1 = _mm_load_ps(src1ptr + ix + 4);
				
				write128<accum>(destptr + ix, F*(C0 + RIGHT(C0, E0) + 
												 C1 + RIGHT(C1, E1)));
				C0 = E0;
				C1 = E1;
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
		
		const __m128 F = _mm_set1_ps(factor);
		
		for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
		{
			const float * const srcptr = srcbase + iy*PITCHIN; 
			float * const destptr = destbase + iy*PITCHOUT;
			
			__m128 C0 = _mm_load_ps(srcptr);
			__m128 C1 = _mm_load_ps(srcptr + PITCHIN);
			
			for(int ix=0; ix<TempPiZSOA_ST::NX; ix+=4)
			{
				const __m128 E0 = _mm_load_ps(srcptr + ix + 4);
				const __m128 E1 = _mm_load_ps(srcptr + ix + PITCHIN + 4);
				
				write128<accum>(destptr + ix, F*(C0 + RIGHT(C0, E0) + 
												 C1 + RIGHT(C1, E1)));
				C0 = E0;
				C1 = E1;
			}
		}
	}
	
public:
	
	DivTensor_SSE(const Real a = 1, const Real dtinvh = 1, const Real h = 1, const Real sigma=1):
	DivTensor_CPP(a, dtinvh, h, sigma)
	{
	}
	
#undef LEFT
#undef RIGHT	
};
