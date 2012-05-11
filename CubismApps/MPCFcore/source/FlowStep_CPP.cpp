/*
 *  FlowStep_CPP.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/6/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include <cmath>
#include <cassert>
#include <xmmintrin.h>
#include <algorithm>
#include <cstdio>
#include <iostream>

using namespace std;

#include "FlowStep_CPP.h"

template<typename X> inline X mysqrt(X x){ abort(); return sqrt(x);}
template<>  inline float mysqrt<float>(float x){ return sqrtf(x);}
template<>  inline double mysqrt<double>(double x){ return sqrt(x);}

float FlowStep_CPP::kB(const bool bVerbose) const
{
	const float size = sizeof(*this)/1024.;
	
	if (bVerbose)
		printf("FlowStep_CPP::kB(): %.2f kiloBytes\n", size);
	
	return size;
}

Real FlowStep_CPP::_getgamma(const Real phi)
{
	const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/m_smoothlength)));
	const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
	const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
	const Real hs = x<0 ? val_xneg : val_xpos;
	
	return m_gamma1*hs + m_gamma2*(((Real)1)-hs);
}

Real FlowStep_CPP::_getPC(const Real phi)
{
	const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/m_smoothlength)));
	const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
	const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
	const Real hs = x<0 ? val_xneg : val_xpos;
	
	return m_pc1*hs + m_pc2*(((Real)1)-hs);
}

struct AIGP { Real r, u, v, w, s, l; }; //Assumed Interleaved Grid Point type

void FlowStep_CPP::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{	
	InputSOA& rho = ringrho.ref(), &u=ringu.ref(), &v=ringv.ref(), &w=ringw.ref(), &p=ringp.ref(), &l = ringls.ref();
	
	for(int sy=0; sy<_BLOCKSIZE_+6; sy++)
		for(int sx=0; sx<_BLOCKSIZE_+6; sx++)
		{
			AIGP pt = *(AIGP*)(gptfirst + gptfloats*(sx + sy*rowgpts));
			
			const int dx = sx-3;
			const int dy = sy-3;
			
			rho.ref(dx, dy) = pt.r;
			u.ref(dx, dy) = pt.u/pt.r;
			v.ref(dx, dy) = pt.v/pt.r;
			w.ref(dx, dy) = pt.w/pt.r;
			p.ref(dx, dy) = (pt.s - (pt.u*pt.u + pt.v*pt.v + pt.w*pt.w)*(((Real)0.5)/pt.r))*(_getgamma(pt.l)-(Real)1) - _getgamma(pt.l)*_getPC(pt.l);
			l.ref(dx, dy) = pt.l;
		}
}

inline Real weno_minus(const Real a, const Real b, const Real c, const Real d, const Real e)
{
	const Real is0 = a*(a*(Real)(4./3.)  - b*(Real)(19./3.)  + c*(Real)(11./3.)) + b*(b*(Real)(25./3.)  - c*(Real)(31./3.)) + c*c*(Real)(10./3.);
	const Real is1 = b*(b*(Real)(4./3.)  - c*(Real)(13./3.)  + d*(Real)(5./3.)) + c*(c*(Real)(13./3.)  - d*(Real)(13./3.)) + d*d*(Real)(4./3.);
	const Real is2 = c*(c*(Real)(10./3.)  - d*(Real)(31./3.)  + e*(Real)(11./3.)) + d*(d*(Real)(25./3.)  - e*(Real)(19./3.)) + e*e*(Real)(4./3.);
	
	const Real is0plus = is0+ (Real)WENOEPS;
	const Real is1plus = is1+ (Real)WENOEPS;
	const Real is2plus = is2+ (Real)WENOEPS;
	
	const Real alpha0 = (Real)(0.1)*((Real)1/(is0plus*is0plus));
	const Real alpha1 = (Real)(0.6)*((Real)1/(is1plus*is1plus));
	const Real alpha2 = (Real)(0.3)*((Real)1/(is2plus*is2plus));
	const Real alphasum = alpha0+alpha1+alpha2;
	
	const Real omega0=alpha0 * (((Real)1)/alphasum);
	const Real omega1=alpha1 * (((Real)1)/alphasum);
	const Real omega2= 1-omega0-omega1;
	
	return omega0*((Real)(1.0/3.)*a-(Real)(7./6.)*b+(Real)(11./6.)*c) + omega1*(-(Real)(1./6.)*b+(Real)(5./6.)*c+(Real)(1./3.)*d) + omega2*((Real)(1./3.)*c+(Real)(5./6.)*d-(Real)(1./6.)*e);
}

inline Real weno_plus(const Real b, const Real c, const Real d, const Real e, const Real f)
{	
	const Real is0 = d*(d*(Real)(10./3.) - e*(Real)(31./3.)  + f*(Real)(11./3.)) + e*(e*(Real)(25./3.)  - f*(Real)(19./3.)) +	f*f*(Real)(4./3.);
	const Real is1 = c*(c*(Real)(4./3.) - d*(Real)(13./3.)  + e*(Real)(5./3.)) + d*(d*(Real)(13./3.)  - e*(Real)(13./3.)) +	e*e*(Real)(4./3.);
	const Real is2 = b*(b*(Real)(4./3.) - c*(Real)(19./3.)  + d*(Real)(11./3.)) + c*(c*(Real)(25./3.)  - d*(Real)(31./3.)) +	d*d*(Real)(10./3.);
	
	const Real is0plus = is0+ (Real)WENOEPS;
	const Real is1plus = is1+ (Real)WENOEPS;
	const Real is2plus = is2+ (Real)WENOEPS;
	
	const Real alpha0 = (Real)(0.1)*(((Real)1)/(is0plus*is0plus));
	const Real alpha1 = (Real)(0.6)*(((Real)1)/(is1plus*is1plus));
	const Real alpha2 = (Real)(0.3)*(((Real)1)/(is2plus*is2plus));
	const Real alphasum = alpha0+alpha1+alpha2;
	
	const Real omega0=alpha0 * (((Real)1)/alphasum);
	const Real omega1=alpha1 * (((Real)1)/alphasum);
	const Real omega2= 1-omega0-omega1;
	
	return omega0*((Real)(1./3.)*f-(Real)(7./6.)*e+(Real)(11./6.)*d) + omega1*(-(Real)(1./6.)*e+(Real)(5./6.)*d+(Real)(1./3.)*c) + omega2*((Real)(1./3.)*d+(Real)(5./6.)*c-(Real)(1./6.)*b);
}

void FlowStep_CPP::_xweno_minus(const InputSOA& _in, TempSOA& _out)
{	
	for(int iy=0; iy<TempSOA::NY; iy++)
	{
		const Real * const in = _in.ptr(-3, iy);
		Real * const out = & _out.ref(0, iy);
		
		for(int ix=0; ix<TempSOA::NX; ix++)
			out[ix] = weno_minus(in[ix+0], in[ix+1], in[ix+2], in[ix+3], in[ix+4]);
	}
}

void FlowStep_CPP::_xweno_pluss(const InputSOA& _in, TempSOA& _out)
{
	for(int iy=0; iy<TempSOA::NY; iy++)
	{
		const Real * const in = _in.ptr(-2, iy);
		Real * const out = & _out.ref(0, iy);
		
		for(int ix=0; ix<TempSOA::NX; ix++)
			out[ix] = weno_plus(in[ix+0], in[ix+1], in[ix+2], in[ix+3], in[ix+4]);
	}
}

void FlowStep_CPP::_yweno_minus(const InputSOA& _in, TempSOA& _out)
{
	static const int L = InputSOA::PITCH;
	const Real * const in = _in.ptr(0,-3);
	Real * out = &_out.ref(0,0);
	
	for(int iy=0; iy<TempSOA::NY; iy++)
	{
		const Real * ptr = &in[iy];

		for(int ix=0; ix<TempSOA::NX; ix++)
			out[ix + iy*TempSOA::PITCH] = weno_minus(ptr[ix*L], ptr[ix*L+L], ptr[ix*L+2*L], ptr[ix*L+3*L], ptr[ix*L+4*L]);
	}
}

void FlowStep_CPP::_yweno_pluss(const InputSOA& _in, TempSOA& _out)
{
	static const int L = InputSOA::PITCH;

	const Real * const in = _in.ptr(0,-2);
	Real * out = &_out.ref(0,0);

	for(int iy=0; iy<TempSOA::NY; iy++)
	{
		const Real * ptr = &in[iy];
		
		for(int ix=0; ix<TempSOA::NX; ix++)
			out[ix + iy*TempSOA::PITCH] = weno_plus(ptr[ix*L], ptr[ix*L+L], ptr[ix*L+2*L], ptr[ix*L+3*L], ptr[ix*L+4*L]);
	}
}

void FlowStep_CPP::_zweno_minus(const int r, const RingInputSOA& in, TempSOA& out)
{
	static const int L = InputSOA::PITCH;
	
	const Real * const a = in(r-3).ptr(0,0);
	const Real * const b = in(r-2).ptr(0,0);
	const Real * const c = in(r-1).ptr(0,0);
	const Real * const d = in(r).ptr(0,0);
	const Real * const e = in(r+1).ptr(0,0);
	
	Real * const o = &out.ref(0,0);
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX-1; ix++)
			o[ix + TempSOA::PITCH*iy] = weno_minus(a[ix+L*iy], b[ix+L*iy], c[ix+L*iy], d[ix+L*iy], e[ix+L*iy]);
}

void FlowStep_CPP::_zweno_pluss(const int r, const RingInputSOA& in, TempSOA& out)
{
	static const int L = InputSOA::PITCH;
	
	const Real * const a = in(r-2).ptr(0,0);
	const Real * const b = in(r-1).ptr(0,0);
	const Real * const c = in(r).ptr(0,0);
	const Real * const d = in(r+1).ptr(0,0);
	const Real * const e = in(r+2).ptr(0,0);
	
	Real * const o = &out.ref(0,0);
	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX-1; ix++)
			o[ix + TempSOA::PITCH*iy] = weno_plus(a[ix+L*iy], b[ix+L*iy], c[ix+L*iy], d[ix+L*iy], e[ix+L*iy]);
}

void FlowStep_CPP::_xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			divu.ref(ix, iy) = um(ix+1, iy) - up(ix, iy);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumls.ref(ix, iy) = lp(ix, iy) + lm(ix+1, iy);
}

void FlowStep_CPP::_yextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			divu.ref(ix, iy) += um(iy+1, ix) - up(iy, ix);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumls.ref(ix, iy) += lp(iy, ix) + lm(iy+1, ix);
}

void FlowStep_CPP::_zextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			divu.ref(ix, iy) += um(ix, iy) - up(ix, iy);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumls.ref(ix, iy) += lp(ix, iy) + lm(ix, iy);
}

void FlowStep_CPP::_char_vel(const TempSOA& rm, const TempSOA& rp, 
							 const TempSOA& vm, const TempSOA& vp,
							 const TempSOA& pm, const TempSOA& pp,
							 const TempSOA& lm, const TempSOA& lp, 
							 TempSOA& outm, TempSOA& outp)
{	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix++)
		{
			const Real cminus = mysqrt(_getgamma(lm(ix, iy))* max((pm(ix, iy)+_getPC(lm(ix, iy)))*((Real)1/rm(ix, iy)), (Real)0));
			const Real cplus  = mysqrt(_getgamma(lp(ix, iy))* max((pp(ix, iy)+_getPC(lp(ix, iy)))*((Real)1/rp(ix, iy)), (Real)0));
			
			outm.ref(ix, iy) = min(vm(ix, iy) - cminus, vm(ix, iy) - cplus);
			outp.ref(ix, iy) = max(vp(ix, iy) + cminus, vp(ix, iy) + cplus);
		}
	
}

void FlowStep_CPP::_hlle_rho(const TempSOA& rm, const TempSOA& rp,
							 const TempSOA& vm, const TempSOA& vp,
							 const TempSOA& am, const TempSOA& ap,
							 TempSOA& out)
{
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix++)
		{
			const bool flagminus = am(ix, iy) > 0;
			const bool flagplus  = ap(ix, iy) < 0; 
			const bool flagother = !(flagminus || flagplus);
			
			const Real fminus = vm(ix, iy)*rm(ix, iy);
			const Real fplus  = vp(ix, iy)*rp(ix, iy);
			
			const Real aminus = am(ix, iy);
			const Real aplus  = ap(ix, iy);
			const Real fother = (aplus*fminus-aminus*fplus+aminus*aplus*(rp(ix, iy)-rm(ix, iy)))*((Real)1/(aplus-aminus));
			
			out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagplus)*fplus + ((Real)flagother)*fother;
		}
}

void FlowStep_CPP::_hlle_vel(const TempSOA& rm, const TempSOA& rp,
							 const TempSOA& vm, const TempSOA& vp,
							 const TempSOA& vdm, const TempSOA& vdp,
							 const TempSOA& am, const TempSOA& ap,
							 TempSOA& out)
{
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix++)
		{
			const bool flagminus = am(ix, iy) > 0;
			const bool flagpluss = ap(ix, iy) < 0; 
			const bool flagother = !(flagminus || flagpluss);
			
			const Real uminus = vm(ix, iy)*rm(ix, iy);
			const Real upluss = vp(ix, iy)*rp(ix, iy);
			
			const Real fminus = vdm(ix, iy)*uminus;
			const Real fpluss = vdp(ix, iy)*upluss;
			
			const Real aminus = am(ix, iy);
			const Real apluss  = ap(ix, iy);
			
			const Real fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))*((Real)1/(apluss-aminus));
			
			out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagpluss)*fpluss + ((Real)flagother)*fother;
		}
}

void FlowStep_CPP::_hlle_pvel(const TempSOA& rm, const TempSOA& rp,
							  const TempSOA& vm, const TempSOA& vp,
							  const TempSOA& pm, const TempSOA& pp,
							  const TempSOA& am, const TempSOA& ap,
							  TempSOA& out)
{
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix++)
		{
			const bool flagminus = am(ix, iy) > 0;
			const bool flagpluss = ap(ix, iy) < 0; 
			const bool flagother = !(flagminus || flagpluss);
			
			const Real myvminus = vm(ix, iy);
			const Real myvpluss = vp(ix, iy);
			
			const Real uminus = myvminus*rm(ix, iy);
			const Real upluss = myvpluss*rp(ix, iy);
			
			const Real fminus = myvminus*uminus + pm(ix, iy);
			const Real fpluss = myvpluss*upluss + pp(ix, iy);
			
			const Real aminus = am(ix, iy);
			const Real apluss  = ap(ix, iy);
			
			const Real fother = (apluss*fminus-aminus*fpluss+aminus*apluss*(upluss-uminus))*((Real)1/(apluss-aminus));
			
			out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagpluss)*fpluss + ((Real)flagother)*fother;
		}
}

void FlowStep_CPP::_hlle_e(const TempSOA& rm, const TempSOA& rp,
						   const TempSOA& vdm, const TempSOA& vdp,
						   const TempSOA& v1m, const TempSOA& v1p,
						   const TempSOA& v2m, const TempSOA& v2p,
						   const TempSOA& pm, const TempSOA& pp,
						   const TempSOA& lm, const TempSOA& lp, 
						   const TempSOA& am, const TempSOA& ap,
						   TempSOA& out)
{	
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix++)
		{
			const bool flagminus = am(ix, iy) > 0;
			const bool flagplus  = ap(ix, iy) < 0; 
			const bool flagother = !(flagminus || flagplus);			
			
			const Real vdminus = vdm(ix, iy);
			const Real v1minus = v1m(ix, iy);
			const Real v2minus = v2m(ix, iy);
			const Real pminus = pm(ix, iy);
			const Real eminus = pminus*(((Real)1)/(_getgamma(lm(ix, iy))-(Real)1)) +
			((Real)0.5)*rm(ix, iy)*(vdminus*vdminus + v1minus*v1minus + v2minus*v2minus);
			
			const Real vdplus = vdp(ix, iy);
			const Real v1plus = v1p(ix, iy);
			const Real v2plus = v2p(ix, iy);
			const Real pplus = pp(ix, iy);
			const Real eplus = pplus*(((Real)1)/(_getgamma(lp(ix, iy))-((Real)1))) +
			((Real)0.5)*rp(ix, iy)*(vdplus*vdplus + v1plus*v1plus + v2plus*v2plus);
			
			const Real fminus = vdminus*(pminus + eminus);
			const Real fpluss = vdplus *(pplus + eplus);
			
			const Real aminus = am(ix, iy);
			const Real aplus  = ap(ix, iy);
			
			const Real fother = (aplus*fminus-aminus*fpluss+aminus*aplus*(eplus-eminus))*((Real)1/(aplus-aminus));
			
			out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagplus)*fpluss + ((Real)flagother)*fother;			
		}
}

void FlowStep_CPP::_xrhsadd(const TempSOA& flux, OutputSOA& rhs)
{
	const Real * const f = flux.ptr(0,0);
	Real * const r = &rhs.ref(0,0);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			r[ix + OutputSOA::PITCH*iy] = f[ix + 1 + TempSOA::PITCH*iy] - f[ix + TempSOA::PITCH*iy];
}

void FlowStep_CPP::_yrhsadd(const TempSOA& flux, OutputSOA& rhs)
{
	const Real * const f = flux.ptr(0,0);
	Real * const r = &rhs.ref(0,0);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			r[ix + OutputSOA::PITCH*iy] += f[iy +  1 + TempSOA::PITCH*ix] - f[iy + TempSOA::PITCH*ix];
}

void FlowStep_CPP::_zrhsadd(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs)
{
	const Real * const ff = fforward.ptr(0,0);
	const Real * const fb = fback.ptr(0,0);

	Real * const r = &rhs.ref(0,0);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			r[ix + OutputSOA::PITCH*iy] += ff[ix +  TempSOA::PITCH*iy] - fb[ix + TempSOA::PITCH*iy];
}

void FlowStep_CPP::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
{	
	const Real factor2 = ((Real)1.)/6;
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
		{
			AIGP& rhs = *(AIGP*)(gptfirst + gptfloats*(ix + iy*rowgpts));
			
			rhs.r = m_a*rhs.r - m_dtinvh*rhsrho(ix, iy);
			rhs.u = m_a*rhs.u - m_dtinvh*rhsu(ix, iy);
			rhs.v = m_a*rhs.v - m_dtinvh*rhsv(ix, iy);
			rhs.w = m_a*rhs.w - m_dtinvh*rhsw(ix, iy);
			rhs.s = m_a*rhs.s - m_dtinvh*rhss(ix, iy);
			rhs.l = m_a*rhs.l - m_dtinvh*(rhsls(ix, iy) - divu(ix,iy)*sumls(ix,iy)*factor2);
		}
}

void FlowStep_CPP::_xflux(const int relid)
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

void FlowStep_CPP::_yflux(const int relid)
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

void FlowStep_CPP::_zflux(const int relid)
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

void FlowStep_CPP::_xrhs()
{	
	_xrhsadd(fluxrho(), rhsrho);
	_xrhsadd(fluxu(), rhsu);
	_xrhsadd(fluxv(), rhsv);
	_xrhsadd(fluxw(), rhsw);
	_xrhsadd(fluxp(), rhss);
	_xrhsadd(fluxls(), rhsls);
}

void FlowStep_CPP::_yrhs()
{	
	_yrhsadd(fluxrho(), rhsrho);
	_yrhsadd(fluxu(), rhsu);
	_yrhsadd(fluxv(), rhsv);
	_yrhsadd(fluxw(), rhsw);
	_yrhsadd(fluxp(), rhss);
	_yrhsadd(fluxls(), rhsls);
}

void FlowStep_CPP::_zrhs()
{	
	_zrhsadd(fluxrho(-1), fluxrho(0), rhsrho);
	_zrhsadd(fluxu(-1), fluxu(0), rhsu);
	_zrhsadd(fluxv(-1), fluxv(0), rhsv);
	_zrhsadd(fluxw(-1), fluxw(0), rhsw);
	_zrhsadd(fluxp(-1), fluxp(0), rhss);
	_zrhsadd(fluxls(-1), fluxls(0), rhsls);
}

void FlowStep_CPP::compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
						   Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
{
	for(int islice=0; islice<5; islice++)
	{
		_convert(srcfirst+islice*srcfloats*slicesrcs, srcfloats, rowsrcs);
		_next();
	}
	
	_convert(srcfirst + 5*srcfloats*slicesrcs, srcfloats, rowsrcs);
	_zflux(-2);
	_flux_next();
	
	for(int islice=0; islice<_BLOCKSIZE_; islice++)
	{
		_xflux(-2);
		_xrhs();
		
		_yflux(-2);
		_yrhs();
		
		_next();
		_convert(srcfirst + (islice+6)*srcfloats*slicesrcs, srcfloats, rowsrcs);
		
		_zflux(-2);
		_zrhs();
		
		_copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
		_flux_next();
	}
}
