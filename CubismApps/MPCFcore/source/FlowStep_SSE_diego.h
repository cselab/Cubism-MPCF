/*
 *  FlowStep_SSE_diego.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/19/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <cstdio>
#include <xmmintrin.h>

#include "FlowStep_CPP.h"

class FlowStep_SSE_diego : public FlowStep_CPP
{
	static const int _CPER16BYTES = 16/sizeof(Real);
	float __attribute__((aligned(16))) tmp[TempSOA::NX][_CPER16BYTES];

protected:
	
	inline __m128 _heaviside(const __m128 phi, const __m128 invh, const __m128 one, const __m128 phalf, const __m128 mhalf) const;
	inline __m128d _heaviside(const __m128d phi, const __m128d invh, const __m128d one, const __m128d phalf, const __m128d mhalf) const;
	
	inline __m128 _getgamma(const __m128 phi, const __m128 inv_smoothlength, 
							const __m128 gamma1, const __m128 gamma2,
							const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const;
	
	inline __m128d _getgamma(const __m128d phi, const __m128d inv_smoothlength, 
							const __m128d gamma1, const __m128d gamma2,
							const __m128d F_1, const __m128d F_1_2, const __m128d M_1_2) const;

    inline __m128 _getPC(const __m128 phi, const __m128 inv_smoothlength, 
							const __m128 pc1, const __m128 pc2,
							const __m128 F_1, const __m128 F_1_2, const __m128 M_1_2) const;
	
	inline __m128d _getPC(const __m128d phi, const __m128d inv_smoothlength, 
                             const __m128d pc1, const __m128d pc2,
                             const __m128d F_1, const __m128d F_1_2, const __m128d M_1_2) const;
    
	void _sse_convert_aligned(const float * const gptfirst, const int gptfloats, const int rowgpts,
					  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l);
	void _sse_convert_aligned(const double * const gptfirst, const int gptfloats, const int rowgpts,
					  double * const rho, double * const u, double * const v, double * const w, double * const p, double * const l);
	
	void _sse_convert(const float * const gptfirst, const int gptfloats, const int rowgpts,
					  float * const rho, float * const u, float * const v, float * const w, float * const p, float * const l);
	void _sse_convert(const double * const gptfirst, const int gptfloats, const int rowgpts,
					  double * const rho, double * const u, double * const v, double * const w, double * const p, double * const l);

	void _sse_xweno_minus(const float * const in, float * const out) const;	
	void _sse_xweno_minus(const double * const in, double * const out) const;	

	void _sse_xweno_pluss(const float * const in, float * const out) const;
	void _sse_xweno_pluss(const double * const in, double * const out) const;

	void _sse_yweno_minus(const float * const in, float * const out);	
	void _sse_yweno_minus(const double * const in, double * const out);	

	void _sse_yweno_pluss(const float * const in, float * const out);	
	void _sse_yweno_pluss(const double * const in, double * const out);	

	void _sse_zweno_minus(const float * const a, const float * const b, const float * const c, const float * const d, const float * const e , float * const out) const;
	void _sse_zweno_minus(const double * const a, const double * const b, const double * const c, const double * const d, const double * const e , double * const out) const;

	void _sse_zweno_pluss(const float * const a, const float * const b, const float * const c, const float * const d, const float * const e , float * const out) const;
	void _sse_zweno_pluss(const double * const a, const double * const b, const double * const c, const double * const d, const double * const e , double * const out) const;
	
       	void _sse_hlle_rho(const float * const rm, const float * const rp,
					   const float * const vm, const float * const vp,
					   const float * const am, const float * const ap,
					   float * const out);
	void _sse_hlle_rho(const double * const rm, const double * const rp,
					   const double * const vm, const double * const vp,
					   const double * const am, const double * const ap,
					   double * const out);
	
	
	void _sse_hlle_vel(const float * const rminus, const float * const rplus,
					   const float * const vminus, const float * const vplus,
					   const float * const vdminus, const float * const vdplus,
					   const float * const aminus, const float * const aplus,
					   float * const out);
	void _sse_hlle_vel(const double * const rminus, const double * const rplus,
					   const double * const vminus, const double * const vplus,
					   const double * const vdminus, const double * const vdplus,
					   const double * const aminus, const double * const aplus,
					   double * const out);
	
	void _sse_hlle_pvel(const float * const rminus, const float * const rplus,
						const float * const vminus, const float * const vplus,
						const float * const pminus, const float * const pplus,
						const float * const aminus, const float * const aplus,
						float * const out);
	void _sse_hlle_pvel(const double * const rminus, const double * const rplus,
						const double * const vminus, const double * const vplus,
						const double * const pminus, const double * const pplus,
						const double * const aminus, const double * const aplus,
						double * const out);
	
	void _sse_hlle_e(const float * const rminus, const float * const rplus,
					 const float * const vdminus, const float * const vdplus,
					 const float * const v1minus, const float * const v1plus,
					 const float * const v2minus, const float * const v2plus,
					 const float * const pminus, const float * const pplus,
					 const float * const lsminus, const float * const lsplus, 
					 const float * const aminus, const float * const aplus,
					 float * const out);
	void _sse_hlle_e(const double * const rminus, const double * const rplus,
					 const double * const vdminus, const double * const vdplus,
					 const double * const v1minus, const double * const v1plus,
					 const double * const v2minus, const double * const v2plus,
					 const double * const pminus, const double * const pplus,
					 const double * const lsminus, const double * const lsplus, 
					 const double * const aminus, const double * const aplus,
					 double * const out);
	
	void _sse_char_vel(const float * const rminus, const float * const rplus, 
					   const float * const vminus, const float * const vplus,
					   const float * const pminus, const float * const pplus,
					   const float * const lminus, const float * const lplus, 
					   float * const out_minus, float * const out_plus);
	void _sse_char_vel(const double * const rminus, const double * const rplus, 
					   const double * const vminus, const double * const vplus,
					   const double * const pminus, const double * const pplus,
					   const double * const lminus, const double * const lplus, 
					   double * const out_minus, double * const out_plus);
	
	void _sse_xrhsadd(const float * const flux, float * const rhs);	
	void _sse_xrhsadd(const double * const flux, double * const rhs);	

	void _sse_yrhsadd(const float * const flux, float * const rhs);
	void _sse_yrhsadd(const double * const flux, double * const rhs);

	void _sse_zrhsadd(const float * const fback, const float * const fforward, float * const rhs);
	void _sse_zrhsadd(const double * const fback, const double * const fforward, double * const rhs);
	
	void _xflux(const int relsliceid);
	void _yflux(const int relsliceid);
	void _zflux(const int relsliceid);
	
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
	
public:
	
	void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
				 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts);
	
	FlowStep_SSE_diego(Real a=0, Real dtinvh=1, Real gamma1=2.5, Real gamma2=2.1, Real smoothlength=1, Real pc1=0, Real pc2=0):
		FlowStep_CPP(a, dtinvh, gamma1, gamma2, smoothlength, pc1, pc2) { }
};
