/*
 *  Convection_AVX.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/3/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <cstdio>
#include <immintrin.h>
#include <xmmintrin.h>

#include "Convection_SSE.h"

class Convection_AVX : public Convection_SSE
{

public:
	
	Convection_AVX(Real a, Real dtinvh, Real gamma1, Real gamma2, Real smoothlength, Real pc1, Real pc2):
	Convection_SSE(a, dtinvh, gamma1, gamma2, smoothlength, pc1, pc2) { }

protected:

	//these methods override the polymorphic ones of Convection_CPP
	void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts, const int slicegpts)
	{
		const bool bAligned = ((size_t)gptfirst & 0x1f) == 0;
		const bool b8Multiple = gptfloats % 8 == 0;
		
		if (bAligned && b8Multiple)
			_avx_convert_aligned(gptfirst, gptfloats, rowgpts, slicegpts, 
								 & ringrho.ref().ref(-8,-3),
								 & ringu.ref().ref(-8,-3),
								 & ringv.ref().ref(-8,-3),
								 & ringw.ref().ref(-8,-3),
								 & ringp.ref().ref(-8,-3),
								 & ringls.ref().ref(-8,-3));
		else
			_avx_convert(gptfirst, gptfloats, rowgpts, slicegpts, 
						 & ringrho.ref().ref(-8,-3),
						 & ringu.ref().ref(-8,-3),
						 & ringv.ref().ref(-8,-3),
						 & ringw.ref().ref(-8,-3),
						 & ringp.ref().ref(-8,-3),
						 & ringls.ref().ref(-8,-3));
	}
	
	void _xweno_minus(const InputSOA& in, TempSOA& out) { _avx_xweno_minus(in.ptr(-8,0), &out.ref(0,0)); }
	void _xweno_pluss(const InputSOA& in, TempSOA& out){ _avx_xweno_pluss(in.ptr(-8,0), &out.ref(0,0)); }
	
	void _yweno_minus(const InputSOA& in, TempSOA& out){ _avx_yweno_minus(in.ptr(0,-3), &out.ref(0,0)); }
	void _yweno_pluss(const InputSOA& in, TempSOA& out){ _avx_yweno_pluss(in.ptr(0,-3), &out.ref(0,0)); }
	
	void _zweno_minus(const int r, const RingInputSOA& in, TempSOA& out)
	{
		_avx_zweno_minus(in(r-3).ptr(0,0), in(r-2).ptr(0,0), in(r-1).ptr(0,0), in(r).ptr(0,0), in(r+1).ptr(0,0), &out.ref(0,0)); 
	}
	
	void _zweno_pluss(const int r, const RingInputSOA& in, TempSOA& out)
	{
		_avx_zweno_pluss(in(r-2).ptr(0,0), in(r-1).ptr(0,0), in(r).ptr(0,0), in(r+1).ptr(0,0), in(r+2).ptr(0,0), &out.ref(0,0)); 
	}
	
	void _hlle_rho(const TempSOA& rm, const TempSOA& rp,
				   const TempSOA& vm, const TempSOA& vp,
				   const TempSOA& am, const TempSOA& ap,
				   TempSOA& out)
	{
		_avx_hlle_rho(rm.ptr(0,0), rp.ptr(0,0), vm.ptr(0,0), vp.ptr(0,0), am.ptr(0,0), ap.ptr(0,0), &out.ref(0,0));
	}
	
	void _hlle_vel(const TempSOA& rminus, const TempSOA& rplus,
				   const TempSOA& vminus, const TempSOA& vplus,
				   const TempSOA& vdminus, const TempSOA& vdplus,
				   const TempSOA& aminus, const TempSOA& aplus,
				   TempSOA& out)
	{
		_avx_hlle_vel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0),
					  vdminus.ptr(0,0), vdplus.ptr(0,0), aminus.ptr(0,0), aplus.ptr(0,0), &out.ref(0,0));
	}
	
	
	void _hlle_pvel(const TempSOA& rminus, const TempSOA& rplus,
					const TempSOA& vminus, const TempSOA& vplus,
					const TempSOA& pminus, const TempSOA& pplus,
					const TempSOA& aminus, const TempSOA& aplus,
					TempSOA& out)
	{
		_avx_hlle_pvel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0), 
					   pminus.ptr(0,0), pplus.ptr(0,0), aminus.ptr(0,0), aplus.ptr(0,0), &out.ref(0,0));
	}
	
	void _hlle_e(const TempSOA& rminus, const TempSOA& rplus,
				 const TempSOA& vdminus, const TempSOA& vdplus,
				 const TempSOA& v1minus, const TempSOA& v1plus,
				 const TempSOA& v2minus, const TempSOA& v2plus,
				 const TempSOA& pminus, const TempSOA& pplus,
				 const TempSOA& lsminus, const TempSOA& lsplus, 
				 const TempSOA& aminus, const TempSOA& aplus,
				 TempSOA& out)
	{
		_avx_hlle_e(rminus.ptr(0,0), rplus.ptr(0,0), vdminus.ptr(0,0), vdplus.ptr(0,0), v1minus.ptr(0,0), v1plus.ptr(0,0),
					v2minus.ptr(0,0), v2plus.ptr(0,0), pminus.ptr(0,0), pplus.ptr(0,0), lsminus.ptr(0,0), lsplus.ptr(0,0), 
					aminus.ptr(0,0), aplus.ptr(0,0), &out.ref(0,0));
	}
	
	void _char_vel(const TempSOA& rminus, const TempSOA& rplus, 
				   const TempSOA& vminus, const TempSOA& vplus,
				   const TempSOA& pminus, const TempSOA& pplus,
				   const TempSOA& lminus, const TempSOA& lplus, 
				   TempSOA& out_minus, TempSOA& out_plus)
	{
		_avx_char_vel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0), 
					  pminus.ptr(0,0), pplus.ptr(0,0), lminus.ptr(0,0), lplus.ptr(0,0), 
					  &out_minus.ref(0,0), &out_plus.ref(0,0));
	}
	
	void _xrhs()
	{	
		_avx_xrhsadd(fluxrho().ptr(0,0), &rhsrho.ref(0,0));
		_avx_xrhsadd(fluxu().ptr(0,0), &rhsu.ref(0,0));
		_avx_xrhsadd(fluxv().ptr(0,0), &rhsv.ref(0,0));
		_avx_xrhsadd(fluxw().ptr(0,0), &rhsw.ref(0,0));
		_avx_xrhsadd(fluxp().ptr(0,0), &rhss.ref(0,0));
		_avx_xrhsadd(fluxls().ptr(0,0), &rhsls.ref(0,0));
	}
	
	void _yrhs()
	{	
		_avx_yrhsadd(fluxrho().ptr(0,0), &rhsrho.ref(0,0));
		_avx_yrhsadd(fluxu().ptr(0,0), &rhsu.ref(0,0));
		_avx_yrhsadd(fluxv().ptr(0,0), &rhsv.ref(0,0));
		_avx_yrhsadd(fluxw().ptr(0,0), &rhsw.ref(0,0));
		_avx_yrhsadd(fluxp().ptr(0,0), &rhss.ref(0,0));
		_avx_yrhsadd(fluxls().ptr(0,0), &rhsls.ref(0,0));
	}
	
	void _zrhs()
	{	
		_avx_zrhsadd(fluxrho(-1).ptr(0,0), fluxrho(0).ptr(0,0), &rhsrho.ref(0,0));
		_avx_zrhsadd(fluxu(-1).ptr(0,0), fluxu(0).ptr(0,0), &rhsu.ref(0,0));
		_avx_zrhsadd(fluxv(-1).ptr(0,0), fluxv(0).ptr(0,0), &rhsv.ref(0,0));
		_avx_zrhsadd(fluxw(-1).ptr(0,0), fluxw(0).ptr(0,0), &rhsw.ref(0,0));
		_avx_zrhsadd(fluxp(-1).ptr(0,0), fluxp(0).ptr(0,0), &rhss.ref(0,0));
		_avx_zrhsadd(fluxls(-1).ptr(0,0), fluxls(0).ptr(0,0), &rhsls.ref(0,0));
	}
	
	//scratchpad data used to transpose the computation in y-direction
	static const int _CPER32BYTES = 32/sizeof(Real);
	float __attribute__((aligned(32))) tmp[TempSOA::NX][_CPER32BYTES];
	
	//these methods are called by the methods declared/defined above
	void _avx_convert_aligned(const float * const gptfirst, const int gptfloats, const int rowgpts, const int slicegpts,
							  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l);
	void _avx_convert(const float * const gptfirst, const int gptfloats, const int rowgpts, const int slicegpts,
					  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l);
	
	void _avx_xweno_minus(const float * const in, float * const out) const;	
	void _avx_xweno_pluss(const float * const in, float * const out) const;
	
	void _avx_yweno_minus(const float * const in, float * const out);	
	void _avx_yweno_pluss(const float * const in, float * const out);	
	
	void _avx_zweno_minus(const float * const a, const float * const b, const float * const c, const float * const d, const float * const e , float * const out) const;
	void _avx_zweno_pluss(const float * const a, const float * const b, const float * const c, const float * const d, const float * const e , float * const out) const;
	
	void _avx_hlle_rho(const float * const rm, const float * const rp,
					   const float * const vm, const float * const vp,
					   const float * const am, const float * const ap,
					   float * const out);
	
	void _avx_hlle_vel(const float * const rminus, const float * const rplus,
					   const float * const vminus, const float * const vplus,
					   const float * const vdminus, const float * const vdplus,
					   const float * const aminus, const float * const aplus,
					   float * const out);
	
	void _avx_hlle_pvel(const float * const rminus, const float * const rplus,
						const float * const vminus, const float * const vplus,
						const float * const pminus, const float * const pplus,
						const float * const aminus, const float * const aplus,
						float * const out);
	
	void _avx_hlle_e(const float * const rminus, const float * const rplus,
					 const float * const vdminus, const float * const vdplus,
					 const float * const v1minus, const float * const v1plus,
					 const float * const v2minus, const float * const v2plus,
					 const float * const pminus, const float * const pplus,
					 const float * const lsminus, const float * const lsplus, 
					 const float * const aminus, const float * const aplus,
					 float * const out);
	
	void _avx_char_vel(const float * const rminus, const float * const rplus, 
					   const float * const vminus, const float * const vplus,
					   const float * const pminus, const float * const pplus,
					   const float * const lminus, const float * const lplus, 
					   float * const out_minus, float * const out_plus);
	
	void _avx_xrhsadd(const float * const flux, float * const rhs);	
	void _avx_yrhsadd(const float * const flux, float * const rhs);
	void _avx_zrhsadd(const float * const fback, const float * const fforward, float * const rhs);	

	inline __m256 _getgamma(const __m256 phi, const __m256 inv_smoothlength, 
							const __m256 gamma1, const __m256 gamma2,
							const __m256 F_1, const __m256 F_1_2, const __m256 M_1_2) const
	{
		return reconstruct(gamma1, gamma2, phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	}
	
    inline __m256 _getPC(const __m256 phi, const __m256 inv_smoothlength, 
						 const __m256 pc1, const __m256 pc2,
						 const __m256 F_1, const __m256 F_1_2, const __m256 M_1_2) const
	{
		return reconstruct(pc1, pc2, phi, inv_smoothlength, F_1, F_1_2, M_1_2);
	}
};
