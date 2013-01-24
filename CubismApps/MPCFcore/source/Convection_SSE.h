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
	
	Convection_SSE(Real a, Real dtinvh) :
	Convection_CPP(a, dtinvh) { }
	
protected:
	
	//these methods override the polymorphic ones of Convection_CPP
	void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
	{
		const bool bAligned = ((size_t)gptfirst & 0xf) == 0;
		const bool b4Multiple = gptfloats % 4 == 0;
		
		if (bAligned && b4Multiple)
			_sse_convert_aligned(gptfirst, gptfloats, rowgpts,
								 & rho.ring.ref().ref(-4,-3),
								 & u.ring.ref().ref(-4,-3),
								 & v.ring.ref().ref(-4,-3),
								 & w.ring.ref().ref(-4,-3),
								 & p.ring.ref().ref(-4,-3),
								 & G.ring.ref().ref(-4,-3)
                                 , & P.ring.ref().ref(-4,-3)
                                 );
		else
			_sse_convert(gptfirst, gptfloats, rowgpts,
						 & rho.ring.ref().ref(-4,-3),
						 & u.ring.ref().ref(-4,-3),
						 & v.ring.ref().ref(-4,-3),
						 & w.ring.ref().ref(-4,-3),
						 & p.ring.ref().ref(-4,-3),
						 & G.ring.ref().ref(-4,-3)
                         , & P.ring.ref().ref(-4,-3)
                         );
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
	
    void _hllc_rho(const TempSOA& rm, const TempSOA& rp,
				   const TempSOA& vm, const TempSOA& vp,
				   const TempSOA& sm, const TempSOA& sp,
                   const TempSOA& star,
				   TempSOA& out)
	{
		_sse_hllc_rho(rm.ptr(0,0), rp.ptr(0,0), vm.ptr(0,0), vp.ptr(0,0), sm.ptr(0,0), sp.ptr(0,0), star.ptr(0,0), &out.ref(0,0));
	}
    
    void _hllc_phi(const TempSOA& phim, const TempSOA& phip,
				   const TempSOA& vm, const TempSOA& vp,
				   const TempSOA& sm, const TempSOA& sp,
                   const TempSOA& star,
				   TempSOA& out)
	{
		_sse_hllc_phi(phim.ptr(0,0), phip.ptr(0,0), vm.ptr(0,0), vp.ptr(0,0), sm.ptr(0,0), sp.ptr(0,0), star.ptr(0,0), &out.ref(0,0));
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
    
    void _hllc_vel(const TempSOA& rminus, const TempSOA& rplus,
				   const TempSOA& vminus, const TempSOA& vplus,
				   const TempSOA& vdminus, const TempSOA& vdplus,
				   const TempSOA& sminus, const TempSOA& splus,
                   const TempSOA& star,
				   TempSOA& out)
	{
		_sse_hllc_vel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0),
					  vdminus.ptr(0,0), vdplus.ptr(0,0), sminus.ptr(0,0), splus.ptr(0,0), star.ptr(0,0), &out.ref(0,0));
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
    
    void _hllc_pvel(const TempSOA& rminus, const TempSOA& rplus,
					const TempSOA& vminus, const TempSOA& vplus,
					const TempSOA& pminus, const TempSOA& pplus,
					const TempSOA& sminus, const TempSOA& splus,
                    const TempSOA& star,
					TempSOA& out)
	{
		_sse_hllc_pvel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0),
					   pminus.ptr(0,0), pplus.ptr(0,0), sminus.ptr(0,0), splus.ptr(0,0), star.ptr(0,0), &out.ref(0,0));
	}
    
	void _hlle_e(const TempSOA& rminus, const TempSOA& rplus,
				 const TempSOA& vdminus, const TempSOA& vdplus,
				 const TempSOA& v1minus, const TempSOA& v1plus,
				 const TempSOA& v2minus, const TempSOA& v2plus,
				 const TempSOA& pminus, const TempSOA& pplus,
				 const TempSOA& Gminus, const TempSOA& Gplus,
                 const TempSOA& Pminus, const TempSOA& Pplus,
				 const TempSOA& aminus, const TempSOA& aplus,
				 TempSOA& out)
	{
		_sse_hlle_e(rminus.ptr(0,0), rplus.ptr(0,0), vdminus.ptr(0,0), vdplus.ptr(0,0), v1minus.ptr(0,0), v1plus.ptr(0,0),
					v2minus.ptr(0,0), v2plus.ptr(0,0), pminus.ptr(0,0), pplus.ptr(0,0), Gminus.ptr(0,0), Gplus.ptr(0,0),
                    Pminus.ptr(0,0), Pplus.ptr(0,0),
					aminus.ptr(0,0), aplus.ptr(0,0), &out.ref(0,0));
	}
    
    void _hllc_e(const TempSOA& rminus, const TempSOA& rplus,
				 const TempSOA& vdminus, const TempSOA& vdplus,
				 const TempSOA& v1minus, const TempSOA& v1plus,
				 const TempSOA& v2minus, const TempSOA& v2plus,
				 const TempSOA& pminus, const TempSOA& pplus,
				 const TempSOA& Gminus, const TempSOA& Gplus,
                 const TempSOA& Pminus, const TempSOA& Pplus,
				 const TempSOA& sminus, const TempSOA& splus,
                 const TempSOA& star,
				 TempSOA& out)
	{
		_sse_hllc_e(rminus.ptr(0,0), rplus.ptr(0,0), vdminus.ptr(0,0), vdplus.ptr(0,0), v1minus.ptr(0,0), v1plus.ptr(0,0),
					v2minus.ptr(0,0), v2plus.ptr(0,0), pminus.ptr(0,0), pplus.ptr(0,0), Gminus.ptr(0,0), Gplus.ptr(0,0),
                    Pminus.ptr(0,0), Pplus.ptr(0,0),
					sminus.ptr(0,0), splus.ptr(0,0), star.ptr(0,0), &out.ref(0,0));
	}
    
	void _char_vel(const TempSOA& rminus, const TempSOA& rplus,
				   const TempSOA& vminus, const TempSOA& vplus,
				   const TempSOA& pminus, const TempSOA& pplus,
				   const TempSOA& Gminus, const TempSOA& Gplus,
                   const TempSOA& Pminus, const TempSOA& Pplus,
				   TempSOA& out_minus, TempSOA& out_plus)
	{
		_sse_char_vel(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0),
					  pminus.ptr(0,0), pplus.ptr(0,0), Gminus.ptr(0,0), Gplus.ptr(0,0),
                      Pminus.ptr(0,0), Pplus.ptr(0,0),
					  &out_minus.ref(0,0), &out_plus.ref(0,0));
	}
	
    void _char_vel_hllc(const TempSOA& rminus, const TempSOA& rplus,
                        const TempSOA& vminus, const TempSOA& vplus,
                        const TempSOA& pminus, const TempSOA& pplus,
                        const TempSOA& Gminus, const TempSOA& Gplus,
                        const TempSOA& Pminus, const TempSOA& Pplus,
                        TempSOA& out_minus, TempSOA& out_plus, TempSOA& out_star)
	{
        _sse_char_vel_hllc(rminus.ptr(0,0), rplus.ptr(0,0), vminus.ptr(0,0), vplus.ptr(0,0),
                           pminus.ptr(0,0), pplus.ptr(0,0), Gminus.ptr(0,0), Gplus.ptr(0,0),
                           Pminus.ptr(0,0), Pplus.ptr(0,0),
                           &out_minus.ref(0,0), &out_plus.ref(0,0), &out_star.ref(0,0));
    }
    
    void _sse_u_hllc(const float * const v_minus, const float * const v_plus, const float * const s_minus, const float * const s_plus, const float * const u_star, float * const out_u_hllc);
    
	void _xrhs()
	{
		_sse_xrhsadd(rho.flux().ptr(0,0), &rho.rhs.ref(0,0));
		_sse_xrhsadd(u.flux().ptr(0,0), &u.rhs.ref(0,0));
		_sse_xrhsadd(v.flux().ptr(0,0), &v.rhs.ref(0,0));
		_sse_xrhsadd(w.flux().ptr(0,0), &w.rhs.ref(0,0));
		_sse_xrhsadd(p.flux().ptr(0,0), &p.rhs.ref(0,0));
		_sse_xrhsadd(G.flux().ptr(0,0), &G.rhs.ref(0,0));
		_sse_xrhsadd(P.flux().ptr(0,0), &P.rhs.ref(0,0));
	}
	
	void _yrhs()
	{
		_sse_yrhsadd(rho.flux().ptr(0,0), &rho.rhs.ref(0,0));
		_sse_yrhsadd(u.flux().ptr(0,0), &u.rhs.ref(0,0));
		_sse_yrhsadd(v.flux().ptr(0,0), &v.rhs.ref(0,0));
		_sse_yrhsadd(w.flux().ptr(0,0), &w.rhs.ref(0,0));
		_sse_yrhsadd(p.flux().ptr(0,0), &p.rhs.ref(0,0));
		_sse_yrhsadd(G.flux().ptr(0,0), &G.rhs.ref(0,0));
		_sse_yrhsadd(P.flux().ptr(0,0), &P.rhs.ref(0,0));
	}
	
	void _zrhs()
	{
		_sse_zrhsadd(rho.flux(-1).ptr(0,0), rho.flux(0).ptr(0,0), &rho.rhs.ref(0,0));
		_sse_zrhsadd(u.flux(-1).ptr(0,0), u.flux(0).ptr(0,0), &u.rhs.ref(0,0));
		_sse_zrhsadd(v.flux(-1).ptr(0,0), v.flux(0).ptr(0,0), &v.rhs.ref(0,0));
		_sse_zrhsadd(w.flux(-1).ptr(0,0), w.flux(0).ptr(0,0), &w.rhs.ref(0,0));
		_sse_zrhsadd(p.flux(-1).ptr(0,0), p.flux(0).ptr(0,0), &p.rhs.ref(0,0));
		_sse_zrhsadd(G.flux(-1).ptr(0,0), G.flux(0).ptr(0,0), &G.rhs.ref(0,0));
		_sse_zrhsadd(P.flux(-1).ptr(0,0), P.flux(0).ptr(0,0), &P.rhs.ref(0,0));
	}
    
	//these methods are called by the methods declared/defined above
	void _sse_convert_aligned(const Real * const gptfirst, const int gptfloats, const int rowgpts,
                              Real * const rho, Real * const u, Real * const v, Real * const w, Real * const p, Real * const G, Real * const P);
    
    void _sse_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts,
                      Real * const rho, Real * const u, Real * const v, Real * const w, Real * const p, Real * const G, Real * const P);
    
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
	
    void _sse_hllc_rho(const Real * const rm, const Real * const rp,
					   const Real * const vm, const Real * const vp,
					   const Real * const sm, const Real * const sp,
                       const Real * const star,
                       Real * const out);
    
    void _sse_hllc_phi(const float * const phim, const float * const phip,
                       const float * const vm, const float * const vp,
                       const float * const sm, const float * const sp,
                       const float * const star,
                       float * const out);
    
	void _sse_hlle_vel(const Real * const rminus, const Real * const rplus,
					   const Real * const vminus, const Real * const vplus,
					   const Real * const vdminus, const Real * const vdplus,
					   const Real * const aminus, const Real * const aplus,
					   Real * const out);
	
    void _sse_hllc_vel(const float * const rm, const float * const rp,
                       const float * const vm, const float * const vp,
                       const float * const vdm, const float * const vdp,
                       const float * const sm, const float * const sp,
                       const float * const star,
                       float * const out);
    
	void _sse_hlle_pvel(const Real * const rminus, const Real * const rplus,
						const Real * const vminus, const Real * const vplus,
						const Real * const pminus, const Real * const pplus,
						const Real * const aminus, const Real * const aplus,
						Real * const out);
	
    void _sse_hllc_pvel(const float * const rm, const float * const rp,
                        const float * const vm, const float * const vp,
                        const float * const pm, const float * const pp,
                        const float * const sm, const float * const sp,
                        const float * const star,
                        float * const out);
    
	void _sse_hlle_e(const Real * const rminus, const Real * const rplus,
					 const Real * const vdminus, const Real * const vdplus,
					 const Real * const v1minus, const Real * const v1plus,
					 const Real * const v2minus, const Real * const v2plus,
					 const Real * const pminus, const Real * const pplus,
					 const Real * const Gminus, const Real * const Gplus,
                     const Real * const Pminus, const Real * const Pplus,
					 const Real * const aminus, const Real * const aplus,
					 Real * const out);
	
    void _sse_hllc_e(const float * const rm, const float * const rp,
                     const float * const vdm, const float * const vdp,
                     const float * const v1m, const float * const v1p,
                     const float * const v2m, const float * const v2p,
                     const float * const pm, const float * const pp,
                     const float * const Gm, const float * const Gp,
                     const float * const Pm, const float * const Pp,
                     const float * const sm, const float * const sp,
                     const float * const star,
                     float * const out);
    
	void _sse_char_vel(const Real * const rminus, const Real * const rplus,
					   const Real * const vminus, const Real * const vplus,
					   const Real * const pminus, const Real * const pplus,
					   const Real * const Gminus, const Real * const Gplus,
                       const Real * const Pminus, const Real * const Pplus,
					   Real * const out_minus, Real * const out_plus);
	
    void _sse_char_vel_hllc(const float * const rm, const float * const rp,
                            const float * const vm, const float * const vp,
                            const float * const pm, const float * const pp,
                            const float * const Gm, const float * const Gp,
                            const float * const Pm, const float * const Pp,
                            float * const outm, float * const outp, float * const out_star);
    
	void _sse_xrhsadd(const Real * const flux, Real * const rhs);
	void _sse_yrhsadd(const Real * const flux, Real * const rhs);
	void _sse_zrhsadd(const Real * const fback, const Real * const fforward, Real * const rhs);
};
