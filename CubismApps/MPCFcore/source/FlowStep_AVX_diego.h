/*
 *  FlowStep_AVX_diego.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/3/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <cstdio>
#include <xmmintrin.h>
#include "FlowStep_SSE_diego.h"

#ifdef _AVX_
#include <immintrin.h>
#endif

class FlowStep_AVX_diego : public FlowStep_SSE_diego
{
	static const int _CPER32BYTES = 32/sizeof(Real);
	float __attribute__((aligned(32))) tmp[TempSOA::NX][_CPER32BYTES];
	
	inline __m256 _heaviside(const __m256 phi, const __m256 invh, const __m256 one, const __m256 phalf, const __m256 mhalf) const;
	
	inline __m256 _getgamma(const __m256 phi, const __m256 inv_smoothlength, 
							const __m256 gamma1, const __m256 gamma2,
							const __m256 F_1, const __m256 F_1_2, const __m256 M_1_2) const;
	
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
	
	void _xflux(const int relid)
	{	
		_xweno_minus(ringrho(relid), wenorho.ref(0));
		_xweno_pluss(ringrho(relid), wenorho.ref(1));
		_xweno_minus(ringu(relid), wenou.ref(0));
		_xweno_pluss(ringu(relid), wenou.ref(1));
		_xweno_minus(ringv(relid), wenov.ref(0));
		_xweno_pluss(ringv(relid), wenov.ref(1));
		_xweno_minus(ringw(relid), wenow.ref(0));
		_xweno_pluss(ringw(relid), wenow.ref(1));
		_xweno_minus(ringp(relid), wenop.ref(0));
		_xweno_pluss(ringp(relid), wenop.ref(1));
		_xweno_minus(ringls(relid), wenols.ref(0));
		_xweno_pluss(ringls(relid), wenols.ref(1));
		
		_xextraterm(wenou(0), wenou(1), wenols(0), wenols(1));	
		_char_vel(wenorho(0), wenorho(1), wenou(0), wenou(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel.ref(0), charvel.ref(1));
		
		_hlle_rho(wenorho(0), wenorho(1), wenou(0), wenou(1), charvel(0), charvel(1), fluxrho.ref());
		_hlle_pvel(wenorho(0), wenorho(1), wenou(0), wenou(1), wenop(0), wenop(1), charvel(0), charvel(1), fluxu.ref());
		_hlle_vel(wenorho(0), wenorho(1), wenov(0), wenov(1), wenou(0), wenou(1), charvel(0), charvel(1), fluxv.ref());
		_hlle_vel(wenorho(0), wenorho(1), wenow(0), wenow(1), wenou(0), wenou(1), charvel(0), charvel(1), fluxw.ref());
		_hlle_e(wenorho(0), wenorho(1), wenou(0), wenou(1), wenov(0), wenov(1), wenow(0), wenow(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel(0), charvel(1), fluxp.ref());
		_hlle_rho(wenols(0), wenols(1), wenou(0), wenou(1), charvel(0), charvel(1), fluxls.ref());
	}
	
	void _yflux(const int relid)
	{	
		_yweno_minus(ringrho(relid), wenorho.ref(0));
		_yweno_pluss(ringrho(relid), wenorho.ref(1));
		_yweno_minus(ringu(relid), wenou.ref(0));
		_yweno_pluss(ringu(relid), wenou.ref(1));
		_yweno_minus(ringv(relid), wenov.ref(0));
		_yweno_pluss(ringv(relid), wenov.ref(1));
		_yweno_minus(ringw(relid), wenow.ref(0));
		_yweno_pluss(ringw(relid), wenow.ref(1));
		_yweno_minus(ringp(relid), wenop.ref(0));
		_yweno_pluss(ringp(relid), wenop.ref(1));
		_yweno_minus(ringls(relid), wenols.ref(0));
		_yweno_pluss(ringls(relid), wenols.ref(1));
		
		_yextraterm(wenov(0), wenov(1), wenols(0), wenols(1));
		_char_vel(wenorho(0), wenorho(1), wenov(0), wenov(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel.ref(0), charvel.ref(1));
		
		_hlle_rho(wenorho(0), wenorho(1), wenov(0), wenov(1), charvel(0), charvel(1), fluxrho.ref());
		_hlle_vel(wenorho(0), wenorho(1), wenou(0), wenou(1), wenov(0), wenov(1), charvel(0), charvel(1), fluxu.ref());
		_hlle_pvel(wenorho(0), wenorho(1), wenov(0), wenov(1), wenop(0), wenop(1), charvel(0), charvel(1), fluxv.ref());
		_hlle_vel(wenorho(0), wenorho(1), wenow(0), wenow(1), wenov(0), wenov(1), charvel(0), charvel(1), fluxw.ref());
		_hlle_e(wenorho(0), wenorho(1), wenov(0), wenov(1), wenou(0), wenou(1), wenow(0), wenow(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel(0), charvel(1), fluxp.ref());
		_hlle_rho(wenols(0), wenols(1), wenov(0), wenov(1), charvel(0), charvel(1), fluxls.ref());
	}
	
	void _zflux(const int relid)
	{	
		_zweno_minus(relid, ringrho, wenorho.ref(0));
		_zweno_pluss(relid, ringrho, wenorho.ref(1));
		_zweno_minus(relid, ringu, wenou.ref(0));
		_zweno_pluss(relid, ringu, wenou.ref(1));
		_zweno_minus(relid, ringv, wenov.ref(0));
		_zweno_pluss(relid, ringv, wenov.ref(1));
		_zweno_minus(relid, ringw, wenow.ref(0));
		_zweno_pluss(relid, ringw, wenow.ref(1));
		_zweno_minus(relid, ringp, wenop.ref(0));
		_zweno_pluss(relid, ringp, wenop.ref(1));
		_zweno_minus(relid, ringls, wenols.ref(0));
		_zweno_pluss(relid, ringls, wenols.ref(1));
		
		_zextraterm(wenow(0), wenow(-1), wenols(0), wenols(-1));
		_char_vel(wenorho(0), wenorho(1), wenow(0), wenow(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel.ref(0), charvel.ref(1));
		
		_hlle_rho(wenorho(0), wenorho(1), wenow(0), wenow(1), charvel(0), charvel(1), fluxrho.ref());
		_hlle_vel(wenorho(0), wenorho(1), wenou(0), wenou(1), wenow(0), wenow(1), charvel(0), charvel(1), fluxu.ref());
		_hlle_vel(wenorho(0), wenorho(1), wenov(0), wenov(1), wenow(0), wenow(1), charvel(0), charvel(1), fluxv.ref());
		_hlle_pvel(wenorho(0), wenorho(1), wenow(0), wenow(1), wenop(0), wenop(1), charvel(0), charvel(1), fluxw.ref());
		_hlle_e(wenorho(0), wenorho(1), wenow(0), wenow(1), wenou(0), wenou(1), wenov(0), wenov(1), wenop(0), wenop(1), wenols(0), wenols(1), charvel(0), charvel(1), fluxp.ref());
		_hlle_rho(wenols(0), wenols(1), wenow(0), wenow(1), charvel(0), charvel(1), fluxls.ref());
	}
	
public:
	
	void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
				 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
	{
		for(int islice=0; islice<5; islice++)
		{
			_convert(srcfirst+islice*srcfloats*slicesrcs, srcfloats, rowsrcs, slicesrcs);
			_next();
		}
		
		_convert(srcfirst + 5*srcfloats*slicesrcs, srcfloats, rowsrcs, slicesrcs);
		_zflux(-2);
		_flux_next();
		
		for(int islice=0; islice<_BLOCKSIZE_; islice++)
		{
			_xflux(-2);
			_xrhs();
			
			_yflux(-2);
			_yrhs();
			
			_next();
			_convert(srcfirst + (islice+6)*srcfloats*slicesrcs, srcfloats, rowsrcs, slicesrcs);
			
			_zflux(-2);
			_zrhs();
			
			_copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
			_flux_next();
		}
	}
	
 FlowStep_AVX_diego(Real a=0, Real dtinvh=1, Real gamma1=2.5, Real gamma2=2.1, Real smoothlength=1, Real pc1=0, Real pc2=0):
	FlowStep_SSE_diego(a, dtinvh, gamma1, gamma2, smoothlength, pc1, pc2) { }
};
