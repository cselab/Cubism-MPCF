/*
 *  Convection_SSE.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/19/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cstdio>
#include <xmmintrin.h>

#include "common.h"
#include "SOA2D.h"
#include "Convection_CPP.h"

class Convection_SSE : public Convection_CPP
{
	static const int _CPER16BYTES = 16/sizeof(Real);
	float __attribute__((aligned(16))) tmp[TempSOA::NX][_CPER16BYTES];
	
public:
	
	Convection_SSE(Real a, Real dtinvh, 
				   Real gamma1, Real gamma2, Real smoothlength, 
				   Real pc1, Real pc2):
	Convection_CPP(a, dtinvh, gamma1, gamma2, smoothlength, pc1, pc2) { }
	
protected:
	
	//these methods override the polymorphic ones of Convection_CPP
	void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
	{
		const bool bAligned = ((size_t)gptfirst & 0xf) == 0;
		const bool b4Multiple = gptfloats % 4 == 0;
		
		if (bAligned && b4Multiple)
			_sse_convert_aligned(gptfirst, gptfloats, rowgpts, 
								 & ringrho.ref().ref(-4,-3),
								 & ringu.ref().ref(-4,-3),
								 & ringv.ref().ref(-4,-3),
								 & ringw.ref().ref(-4,-3),
								 & ringp.ref().ref(-4,-3),
								 & ringls.ref().ref(-4,-3));
		else
			_sse_convert(gptfirst, gptfloats, rowgpts, 
						 & ringrho.ref().ref(-4,-3),
						 & ringu.ref().ref(-4,-3),
						 & ringv.ref().ref(-4,-3),
						 & ringw.ref().ref(-4,-3),
						 & ringp.ref().ref(-4,-3),
						 & ringls.ref().ref(-4,-3));
	}
	
	void _xweno_minus(const InputSOA& in, TempSOA& out) { _sse_xweno_minus(in.ptr(-4,0), &out.ref(0,0)); }
	void _xweno_pluss(const InputSOA& in, TempSOA& out){ _sse_xweno_pluss(in.ptr(-4,0), &out.ref(0,0)); }
	void _yweno_minus(const InputSOA& in, TempSOA& out){ _sse_yweno_minus(in.ptr(0,-3), &out.ref(0,0)); }
	void _yweno_pluss(const InputSOA& in, TempSOA& out){ _sse_yweno_pluss(in.ptr(0,-3), &out.ref(0,0)); }
	
	void _zweno_minus(const int r, const RingInputSOA& in, TempSOA& out)
	{
		_sse_zweno_minus(in(r-3).ptr(0,0), in(r-2).ptr(0,0), in(r-1).ptr(0,0), in(r).ptr(0,0), in(r+1).ptr(0,0), &out.ref(0,0)); 
	}
	
	void _zweno_pluss(const int r, const RingInputSOA& in, TempSOA& out)
	{
		_sse_zweno_pluss(in(r-2).ptr(0,0), in(r-1).ptr(0,0), in(r).ptr(0,0), in(r+1).ptr(0,0), in(r+2).ptr(0,0), &out.ref(0,0)); 
	}
	
	void _hlle_rho(const TempSOA& rm, const TempSOA& rp,
				   const TempSOA& vm, const TempSOA& vp,
				   const TempSOA& am, const TempSOA& ap,
				   TempSOA& out)
	{
		_sse_hlle_rho(rm.ptr(0,0), rp.ptr(0,0), vm.ptr(0,0), vp.ptr(0,0), am.ptr(0,0), ap.ptr(0,0), &out.ref(0,0));
	}
	
	void _hlle_vel(const TempSOA& rminus, const TempSOA& rplus,
				   const TempSOA& vminus, const TempSOA& vplus,
				   const TempSOA& vdminus, const TempSOA& vdplus,
				   const TempSOA& aminus, const TempSOA& aplus,
				   TempSOA& out)
	{
		_sse_hlle_vel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0),
					  vdminus.ptr(0,0), vdplus.ptr(0,0), aminus.ptr(0,0), aplus.ptr(0,0), &out.ref(0,0));
	}
	
	
	void _hlle_pvel(const TempSOA& rminus, const TempSOA& rplus,
					const TempSOA& vminus, const TempSOA& vplus,
					const TempSOA& pminus, const TempSOA& pplus,
					const TempSOA& aminus, const TempSOA& aplus,
					TempSOA& out)
	{
		_sse_hlle_pvel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0), 
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
		_sse_hlle_e(rminus.ptr(0,0), rplus.ptr(0,0), vdminus.ptr(0,0), vdplus.ptr(0,0), v1minus.ptr(0,0), v1plus.ptr(0,0),
					v2minus.ptr(0,0), v2plus.ptr(0,0), pminus.ptr(0,0), pplus.ptr(0,0), lsminus.ptr(0,0), lsplus.ptr(0,0), 
					aminus.ptr(0,0), aplus.ptr(0,0), &out.ref(0,0));
	}
	
	void _char_vel(const TempSOA& rminus, const TempSOA& rplus, 
				   const TempSOA& vminus, const TempSOA& vplus,
				   const TempSOA& pminus, const TempSOA& pplus,
				   const TempSOA& lminus, const TempSOA& lplus, 
				   TempSOA& out_minus, TempSOA& out_plus)
	{
		_sse_char_vel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0), 
					  pminus.ptr(0,0), pplus.ptr(0,0), lminus.ptr(0,0), lplus.ptr(0,0), 
					  &out_minus.ref(0,0), &out_plus.ref(0,0));
	}
	
	void _xrhs()
	{	
		_sse_xrhsadd(fluxrho().ptr(0,0), &rhsrho.ref(0,0));
		_sse_xrhsadd(fluxu().ptr(0,0), &rhsu.ref(0,0));
		_sse_xrhsadd(fluxv().ptr(0,0), &rhsv.ref(0,0));
		_sse_xrhsadd(fluxw().ptr(0,0), &rhsw.ref(0,0));
		_sse_xrhsadd(fluxp().ptr(0,0), &rhss.ref(0,0));
		_sse_xrhsadd(fluxls().ptr(0,0), &rhsls.ref(0,0));
	}
	
	void _yrhs()
	{	
		_sse_yrhsadd(fluxrho().ptr(0,0), &rhsrho.ref(0,0));
		_sse_yrhsadd(fluxu().ptr(0,0), &rhsu.ref(0,0));
		_sse_yrhsadd(fluxv().ptr(0,0), &rhsv.ref(0,0));
		_sse_yrhsadd(fluxw().ptr(0,0), &rhsw.ref(0,0));
		_sse_yrhsadd(fluxp().ptr(0,0), &rhss.ref(0,0));
		_sse_yrhsadd(fluxls().ptr(0,0), &rhsls.ref(0,0));
	}
	
	void _zrhs()
	{	
		_sse_zrhsadd(fluxrho(-1).ptr(0,0), fluxrho(0).ptr(0,0), &rhsrho.ref(0,0));
		_sse_zrhsadd(fluxu(-1).ptr(0,0), fluxu(0).ptr(0,0), &rhsu.ref(0,0));
		_sse_zrhsadd(fluxv(-1).ptr(0,0), fluxv(0).ptr(0,0), &rhsv.ref(0,0));
		_sse_zrhsadd(fluxw(-1).ptr(0,0), fluxw(0).ptr(0,0), &rhsw.ref(0,0));
		_sse_zrhsadd(fluxp(-1).ptr(0,0), fluxp(0).ptr(0,0), &rhss.ref(0,0));
		_sse_zrhsadd(fluxls(-1).ptr(0,0), fluxls(0).ptr(0,0), &rhsls.ref(0,0));
	}

	//these methods are called by the methods declared/defined above
	void _sse_convert_aligned(const Real * const gptfirst, const int gptfloats, const int rowgpts,
							  Real * const rho, Real * const u, Real * const v, Real * const w, Real * const p, Real * const l);
	
	void _sse_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts,
					  Real * const rho, Real * const u, Real * const v, Real * const w, Real * const p, Real * const l);
	
	void _sse_xweno_minus(const Real * const in, Real * const out) const;	
	void _sse_xweno_pluss(const Real * const in, Real * const out) const;
	void _sse_yweno_minus(const Real * const in, Real * const out);	
	void _sse_yweno_pluss(const Real * const in, Real * const out);	
	void _sse_zweno_minus(const Real * const a, const Real * const b, const Real * const c, const Real * const d, const Real * const e , Real * const out) const;
	void _sse_zweno_pluss(const Real * const a, const Real * const b, const Real * const c, const Real * const d, const Real * const e , Real * const out) const;
	
	void _sse_hlle_rho(const Real * const rm, const Real * const rp,
					   const Real * const vm, const Real * const vp,
					   const Real * const am, const Real * const ap,
					   Real * const out);
	
	void _sse_hlle_vel(const Real * const rminus, const Real * const rplus,
					   const Real * const vminus, const Real * const vplus,
					   const Real * const vdminus, const Real * const vdplus,
					   const Real * const aminus, const Real * const aplus,
					   Real * const out);
	
	void _sse_hlle_pvel(const Real * const rminus, const Real * const rplus,
						const Real * const vminus, const Real * const vplus,
						const Real * const pminus, const Real * const pplus,
						const Real * const aminus, const Real * const aplus,
						Real * const out);
	
	void _sse_hlle_e(const Real * const rminus, const Real * const rplus,
					 const Real * const vdminus, const Real * const vdplus,
					 const Real * const v1minus, const Real * const v1plus,
					 const Real * const v2minus, const Real * const v2plus,
					 const Real * const pminus, const Real * const pplus,
					 const Real * const lsminus, const Real * const lsplus, 
					 const Real * const aminus, const Real * const aplus,
					 Real * const out);
	
	void _sse_char_vel(const Real * const rminus, const Real * const rplus, 
					   const Real * const vminus, const Real * const vplus,
					   const Real * const pminus, const Real * const pplus,
					   const Real * const lminus, const Real * const lplus, 
					   Real * const out_minus, Real * const out_plus);
	
	void _sse_xrhsadd(const Real * const flux, Real * const rhs);	
	void _sse_yrhsadd(const Real * const flux, Real * const rhs);
	void _sse_zrhsadd(const Real * const fback, const Real * const fforward, Real * const rhs);
};
