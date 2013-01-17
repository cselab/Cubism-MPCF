/*
 *  Convection_CPP.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/6/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cmath>
#include <cassert>

#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>

using namespace std;

#include "Convection_CPP.h"

Convection_CPP::Convection_CPP(const Real a, const Real dtinvh): a(a), dtinvh(dtinvh) { }

void Convection_CPP::compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
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

void Convection_CPP::hpc_info(float& flop_convert, int& traffic_convert,
							  float& flop_weno, int& traffic_weno,
							  float& flop_extraterm, int& traffic_extraterm,
							  float& flop_charvel, int& traffic_charvel,
							  float& flop_hlle, int& traffic_hlle,
							  float& flop_div, int& traffic_div,
							  float& flop_copyback, int& traffic_copyback,
							  size_t& footprint)
{
	const int ninputs = (int)powf(_BLOCKSIZE_ + 6, 3);
	const int nfaces = (int)powf(_BLOCKSIZE_, 2) * (_BLOCKSIZE_ + 1);
	const int ncells = (int)powf(_BLOCKSIZE_, 3);
	const int nquantities = 6;
	const int ndirections = 3;
	
	flop_convert = (3 + 12 + 3 * 15) * ninputs;
	traffic_convert = (6 + 6) * sizeof(Real) * ninputs;
	flop_weno =  82 * 2 * nfaces * nquantities * ndirections;
	traffic_weno = (5 + 1) * sizeof(Real) * (2 * nfaces * nquantities * ndirections);
	flop_extraterm = 2 * ndirections * ncells;
	traffic_extraterm = (2 + 4) * sizeof(Real) * ndirections * ncells;
	flop_charvel = (18 + 15 * 4) * ndirections * nfaces;
	traffic_charvel = (8 + 2) * sizeof(Real) * ndirections * nfaces;
	flop_hlle = (17 * 2 + 19 * 2 + 21 + 73) * ndirections * nfaces;
	traffic_hlle = ((6 + 1) * 2 + (8 + 1) * 2 + (8 + 1) + (14 + 1)) * sizeof(Real) * ndirections * nfaces;
	flop_div = (1 + 2 * 2) * ncells;
	traffic_div = (2 + 1) * sizeof(Real) * ndirections * ncells;
	flop_copyback = (6 * 3 + 3) * ncells;
	traffic_copyback = (8 + 6) * sizeof(Real) * ncells;
	footprint = sizeof(Convection_CPP);
}

void Convection_CPP::printflops(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES,
								const int NT, const int NBLOCKS, const float MEASUREDTIME)
{
	const float PEAKPERF = PEAKPERF_CORE*NCORES;
	
	float flop_convert, flop_weno, flop_extraterm, flop_charvel, flop_hlle, flop_div, flop_copyback;
	int traffic_convert, traffic_weno, traffic_extraterm, traffic_charvel, traffic_hlle, traffic_div, traffic_copyback;
	size_t footprint;
	
	hpc_info(flop_convert, traffic_convert,
			 flop_weno, traffic_weno,
			 flop_extraterm, traffic_extraterm,
			 flop_charvel, traffic_charvel,
			 flop_hlle, traffic_hlle,
			 flop_div, traffic_div,
			 flop_copyback, traffic_copyback,
			 footprint);
	
	double texpected_ai = 0;
	double totflop = 0;
	
	//compute texpected_ai
	{
		std::vector<float> ai(7), flop(7);
		
		ai[0] = flop_convert/traffic_convert;
		ai[1] = flop_weno/traffic_weno;
		ai[2] = flop_extraterm/traffic_extraterm;
		ai[3] = flop_charvel/traffic_charvel;
		ai[4] = flop_hlle/traffic_hlle;
		ai[5] = flop_div/traffic_div;
		ai[6] = flop_copyback/traffic_copyback;
		
		flop[0] = flop_convert;
		flop[1] = flop_weno;
		flop[2] = flop_extraterm;
		flop[3] = flop_charvel;
		flop[4] = flop_hlle;
		flop[5] = flop_div;
		flop[6] = flop_copyback;
		
		for(int i=0; i<ai.size(); ++i)
			texpected_ai += NT * NBLOCKS * flop[i] / min(PEAKPERF, PEAKBAND*ai[i]);
		
		for(int i=0; i<ai.size(); ++i)
			totflop += NT * NBLOCKS * flop[i];
	}
	
	const double ai_overall = min((double)PEAKPERF , (totflop / texpected_ai) / PEAKBAND);
    const int nquantities = 7;
    
	const double inout_footprint =  NT * NBLOCKS * nquantities * (size_t)sizeof(Real) *
	(powf(_BLOCKSIZE_ + 6, 3) + 2 * powf(_BLOCKSIZE_, 3));
	
	const double oi_overall = totflop/(inout_footprint + NT * (2 + 1) * footprint);
	const double texpected_oi = totflop/min((double)PEAKPERF, PEAKBAND*oi_overall);
	
	const double perf_measured = 1e-9*totflop/MEASUREDTIME;
	
	printPerformanceTitle();
	printf("\tINTERMEDIATE MEMORY FOOTPRINT: %.4f MB\tTOTAL TRAFFIC: %.4f MB\n", footprint/1024./1024, (NT * NBLOCKS * (2 + 1) * footprint + inout_footprint)/1024./1024);
	printf("\tASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
	printf("\tRIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
	printf("\tCONVECTION THIS ONE IS %.2f GFLOP/s,\t\"per block\" %.2f FLOP/B [AI] - %.2f FLOP/B [OI]\n", perf_measured, ai_overall, oi_overall);
	printf("\tTIME PER BLOCK: %.5f ms (expected %.5f [AI] - %.5f [OI] ms)\n",  1e3*MEASUREDTIME/(NT * NBLOCKS), 1e3*texpected_ai/(NT * NBLOCKS), 1e3*texpected_oi/(NT * NBLOCKS));
	printf("\tExpected Performance is: %.2f  GFLOP/s [AI], %.2f  GFLOP/s [OI]\n", totflop*1e-9/texpected_ai, totflop*1e-9/texpected_oi);
	printf("\tEFFICIENCY: %.2f%% [AI] - %.2f%% [OI], HW-UTILIZATION: %.2f%%\n", 100.*min(1., texpected_ai/MEASUREDTIME), 100.*texpected_oi/MEASUREDTIME, 100*perf_measured*1e9/PEAKPERF);
	printEndLine();
}

void Convection_CPP::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
	InputSOA& rho = this->rho.ring.ref(), &u = this->u.ring.ref(), &v = this->v.ring.ref(),
	&w = this->w.ring.ref(), &p = this->p.ring.ref(), &G = this->G.ring.ref();
    
	InputSOA& P = this->P.ring.ref();
	
	for(int sy=0; sy<_BLOCKSIZE_+6; sy++)
		for(int sx=0; sx<_BLOCKSIZE_+6; sx++)
		{
			AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));
			
			const int dx = sx-3;
			const int dy = sy-3;
			
			rho.ref(dx, dy) = pt.r;
			u.ref(dx, dy) = pt.u/pt.r;
			v.ref(dx, dy) = pt.v/pt.r;
			w.ref(dx, dy) = pt.w/pt.r;
			p.ref(dx, dy) = (pt.s - (pt.u*pt.u + pt.v*pt.v + pt.w*pt.w)*(((Real)0.5)/pt.r)-pt.P)/pt.G;
			G.ref(dx, dy) = pt.G;
			P.ref(dx, dy) = pt.P;
			
			/*assert(!isnan(rho.ref(dx, dy)));
             assert(!isnan(u.ref(dx, dy)));
             assert(!isnan(v.ref(dx, dy)));
             assert(!isnan(w.ref(dx, dy)));
             assert(!isnan(p.ref(dx, dy)));
             assert(!isnan(G.ref(dx, dy)));*/
			//assert(!isnan(P.ref(dx, dy)));
		}
}

inline Real weno_minus(const Real a, const Real b, const Real c, const Real d, const Real e) //82 FLOP
{
	const Real is0 = a*(a*(Real)(4./3.)  - b*(Real)(19./3.)  + c*(Real)(11./3.)) + b*(b*(Real)(25./3.)  - c*(Real)(31./3.)) + c*c*(Real)(10./3.);
	const Real is1 = b*(b*(Real)(4./3.)  - c*(Real)(13./3.)  + d*(Real)(5./3.))  + c*(c*(Real)(13./3.)  - d*(Real)(13./3.)) + d*d*(Real)(4./3.);
	const Real is2 = c*(c*(Real)(10./3.) - d*(Real)(31./3.)  + e*(Real)(11./3.)) + d*(d*(Real)(25./3.)  - e*(Real)(19./3.)) + e*e*(Real)(4./3.);
	
	const Real is0plus = is0 + (Real)WENOEPS;
	const Real is1plus = is1 + (Real)WENOEPS;
	const Real is2plus = is2 + (Real)WENOEPS;
	
	const Real alpha0 = (Real)(0.1)*((Real)1/(is0plus*is0plus));
	const Real alpha1 = (Real)(0.6)*((Real)1/(is1plus*is1plus));
	const Real alpha2 = (Real)(0.3)*((Real)1/(is2plus*is2plus));
	const Real alphasum = alpha0+alpha1+alpha2;
	
	const Real omega0=alpha0 * (((Real)1)/alphasum);
	const Real omega1=alpha1 * (((Real)1)/alphasum);
	const Real omega2= 1-omega0-omega1;
	
	return omega0*((Real)(1.0/3.)*a-(Real)(7./6.)*b+(Real)(11./6.)*c) + omega1*(-(Real)(1./6.)*b+(Real)(5./6.)*c+(Real)(1./3.)*d) + omega2*((Real)(1./3.)*c+(Real)(5./6.)*d-(Real)(1./6.)*e);
}

inline Real weno_plus(const Real b, const Real c, const Real d, const Real e, const Real f) //82 FLOP
{
	const Real is0 = d*(d*(Real)(10./3.)- e*(Real)(31./3.) + f*(Real)(11./3.)) + e*(e*(Real)(25./3.) - f*(Real)(19./3.)) +	f*f*(Real)(4./3.);
	const Real is1 = c*(c*(Real)(4./3.) - d*(Real)(13./3.) + e*(Real)(5./3.)) + d*(d*(Real)(13./3.)  - e*(Real)(13./3.)) +	e*e*(Real)(4./3.);
	const Real is2 = b*(b*(Real)(4./3.) - c*(Real)(19./3.) + d*(Real)(11./3.)) + c*(c*(Real)(25./3.) - d*(Real)(31./3.)) +	d*d*(Real)(10./3.);
	
	const Real is0plus = is0 + (Real)WENOEPS;
	const Real is1plus = is1 + (Real)WENOEPS;
	const Real is2plus = is2 + (Real)WENOEPS;
	
	const Real alpha0 = (Real)(0.1)*(((Real)1)/(is0plus*is0plus));
	const Real alpha1 = (Real)(0.6)*(((Real)1)/(is1plus*is1plus));
	const Real alpha2 = (Real)(0.3)*(((Real)1)/(is2plus*is2plus));
	const Real alphasum = alpha0+alpha1+alpha2;
	
	const Real omega0=alpha0 * (((Real)1)/alphasum);
	const Real omega1=alpha1 * (((Real)1)/alphasum);
	const Real omega2= 1-omega0-omega1;
	
	return omega0*((Real)(1./3.)*f-(Real)(7./6.)*e+(Real)(11./6.)*d) + omega1*(-(Real)(1./6.)*e+(Real)(5./6.)*d+(Real)(1./3.)*c) + omega2*((Real)(1./3.)*d+(Real)(5./6.)*c-(Real)(1./6.)*b);
}

void Convection_CPP::_xweno_minus(const InputSOA& _in, TempSOA& _out)
{
	for(int iy=0; iy<TempSOA::NY; iy++)
	{
		const Real * const in = _in.ptr(-3, iy);
		Real * const out = & _out.ref(0, iy);
		
		for(int ix=0; ix<TempSOA::NX; ix++)
			out[ix] = weno_minus(in[ix+0], in[ix+1], in[ix+2], in[ix+3], in[ix+4]);
	}
}

void Convection_CPP::_xweno_pluss(const InputSOA& _in, TempSOA& _out)
{
	for(int iy=0; iy<TempSOA::NY; iy++)
	{
		const Real * const in = _in.ptr(-2, iy);
		Real * const out = & _out.ref(0, iy);
		
		for(int ix=0; ix<TempSOA::NX; ix++)
			out[ix] = weno_plus(in[ix+0], in[ix+1], in[ix+2], in[ix+3], in[ix+4]);
	}
}

void Convection_CPP::_yweno_minus(const InputSOA& _in, TempSOA& _out)
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

void Convection_CPP::_yweno_pluss(const InputSOA& _in, TempSOA& _out)
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

void Convection_CPP::_zweno_minus(const int r, const RingInputSOA& in, TempSOA& out)
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

void Convection_CPP::_zweno_pluss(const int r, const RingInputSOA& in, TempSOA& out)
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

void Convection_CPP::_xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& Gp
								 , const TempSOA& Pm, const TempSOA& Pp
                                 , const TempSOA& am, const TempSOA& ap)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			divu.ref(ix, iy) = (ap(ix+1, iy)*um(ix+1, iy)-am(ix+1, iy)*up(ix+1, iy))/(ap(ix+1, iy)-am(ix+1, iy))-(ap(ix, iy)*um(ix, iy)-am(ix, iy)*up(ix, iy))/(ap(ix, iy)-am(ix, iy));
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumG.ref(ix, iy) = Gp(ix, iy) + Gm(ix+1, iy);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumP.ref(ix, iy) = Pp(ix, iy) + Pm(ix+1, iy);
}

Real Convection_CPP::_u_hllc(const Real rho_minus, const Real rho_plus, const Real v_minus, const Real v_plus
                             , const Real p_minus, const Real p_plus, const Real G_minus, const Real G_plus)
{
    const Real a_minus = std::sqrt((1/G_minus+1)* max((p_minus)*((Real)1/rho_minus), (Real)0));//NOT CONSIDERING Pc HERE
    const Real a_plus  = std::sqrt((1/G_plus+1)* max((p_plus)*((Real)1/rho_plus), (Real)0));//NOT CONSIDERING Pc HERE
    
    const Real rho_hat = 0.5*(rho_minus+rho_plus);
    const Real a_hat   = 0.5*(a_minus+a_plus);
    const Real rho_hat_a_hat = rho_hat*a_hat;
    
    const Real u_star = 0.5*(v_minus+v_plus+(p_minus-p_plus)/rho_hat_a_hat);
    const Real p_star = 0.5*(p_minus+p_plus)+0.5*(v_minus-v_plus)*rho_hat_a_hat;
    
    const Real q_minus = (p_star <= p_minus)? 1 : sqrt(1+(2*G_minus+1)/(G_minus+1)*(p_star/p_minus-1));//NOT CONSIDERING Pc HERE
    const Real q_plus  = (p_star <= p_plus )? 1 : sqrt(1+(2*G_plus +1)/(G_plus +1)*(p_star/p_plus -1));//NOT CONSIDERING Pc HERE
    
    const Real s_minus = std::min(v_minus - a_minus*q_minus, v_plus - a_plus*q_plus);
    const Real s_plus  = std::max(v_minus + a_minus*q_minus, v_plus + a_plus*q_plus);
    
    const Real signum = (u_star > 0) ? 1 : ((u_star < 0) ? -1 : 0);
    
    const Real xi_minus = (s_minus-v_minus)/(s_minus-u_star);
    const Real xi_plus  = (s_plus -v_plus )/(s_plus -u_star);

    const Real term1 = 0.5*(1+signum)*(v_minus+s_minus*(xi_minus-1));
    const Real term2 = 0.5*(1-signum)*(v_plus +s_plus *(xi_plus -1));
    
    return (term1+term2);
}

void Convection_CPP::_xextraterm_hllc(const TempSOA& rm, const TempSOA& rp
                                      , const TempSOA& vm, const TempSOA& vp
                                      , const TempSOA& pm, const TempSOA& pp
                                      , const TempSOA& Gm, const TempSOA& Gp
                                      , const TempSOA& Pm, const TempSOA& Pp)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
        {
            Real rho_minus = rm(ix, iy);
			Real rho_plus  = rp(ix, iy);
            
            Real v_minus = vm(ix, iy);
			Real v_plus  = vp(ix, iy);
            
            Real p_minus = pm(ix, iy);
			Real p_plus  = pp(ix, iy);
            
            Real G_minus = Gm(ix, iy);
			Real G_plus  = Gp(ix, iy);
            
            const Real u_hllc_minus = _u_hllc(rho_minus, rho_plus, v_minus, v_plus, p_minus, p_plus, G_minus, G_plus);

            rho_minus = rm(ix+1, iy);
			rho_plus  = rp(ix+1, iy);
            
            v_minus = vm(ix+1, iy);
			v_plus  = vp(ix+1, iy);
            
            p_minus = pm(ix+1, iy);
			p_plus  = pp(ix+1, iy);
            
            G_minus = Gm(ix+1, iy);
            G_plus  = Gp(ix+1, iy);
            
            const Real u_hllc_plus = _u_hllc(rho_minus, rho_plus, v_minus, v_plus, p_minus, p_plus, G_minus, G_plus);
            
			divu.ref(ix, iy) = u_hllc_plus-u_hllc_minus;
        }
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumG.ref(ix, iy) = Gp(ix, iy) + Gm(ix+1, iy);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumP.ref(ix, iy) = Pp(ix, iy) + Pm(ix+1, iy);
}

void Convection_CPP::_yextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& Gp
								 , const TempSOA& Pm, const TempSOA& Pp
                                 , const TempSOA& am, const TempSOA& ap)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			divu.ref(ix, iy) += (ap(iy+1, ix)*um(iy+1, ix)-am(iy+1, ix)*up(iy+1, ix))/(ap(iy+1, ix)-am(iy+1, ix))-(ap(iy, ix)*um(iy, ix)-am(iy, ix)*up(iy, ix))/(ap(iy, ix)-am(iy, ix));
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumG.ref(ix, iy) += Gp(iy, ix) + Gm(iy+1, ix);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumP.ref(ix, iy) += Pp(iy, ix) + Pm(iy+1, ix);
}

void Convection_CPP::_zextraterm(const TempSOA& um0, const TempSOA& up0, const TempSOA& um1, const TempSOA& up1, const TempSOA& Gm, const TempSOA& Gp
								 , const TempSOA& Pm, const TempSOA& Pp
                                 , const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1)
{
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			divu.ref(ix, iy) += (ap1(ix, iy)*um1(ix, iy)-am1(ix, iy)*up1(ix, iy))/(ap1(ix, iy)-am1(ix, iy))-(ap0(ix, iy)*um0(ix, iy)-am0(ix, iy)*up0(ix, iy))/(ap0(ix, iy)-am0(ix, iy));
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumG.ref(ix, iy) += Gp(ix, iy) + Gm(ix, iy);
	
	for(int iy=0; iy<OutputSOA::NY; iy++)
		for(int ix=0; ix<OutputSOA::NX; ix++)
			sumP.ref(ix, iy) += Pp(ix, iy) + Pm(ix, iy);
}

void Convection_CPP::_char_vel(const TempSOA& rm, const TempSOA& rp,
							   const TempSOA& vm, const TempSOA& vp,
							   const TempSOA& pm, const TempSOA& pp,
							   const TempSOA& Gm, const TempSOA& Gp,
   							   const TempSOA& Pm, const TempSOA& Pp,
							   TempSOA& outm, TempSOA& outp)
{
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix++)
		{
			const Real cminus = std::sqrt((1/Gm(ix,iy)+1)* max((pm(ix, iy)+Pm(ix,iy)/Gm(ix,iy)/(1/Gm(ix,iy)+1))*((Real)1/rm(ix, iy)), (Real)0));
			const Real cplus  = std::sqrt((1/Gp(ix,iy)+1)* max((pp(ix, iy)+Pp(ix,iy)/Gp(ix,iy)/(1/Gp(ix,iy)+1))*((Real)1/rp(ix, iy)), (Real)0));
            
            outm.ref(ix, iy) = min(vm(ix, iy) - cminus, vp(ix, iy) - cplus);
            outp.ref(ix, iy) = max(vm(ix, iy) + cminus, vp(ix, iy) + cplus);
		}
	
}

void Convection_CPP::_hlle_rho(const TempSOA& rm, const TempSOA& rp,
							   const TempSOA& vm, const TempSOA& vp,
							   const TempSOA& am, const TempSOA& ap,
							   TempSOA& out) //17 FLOP
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

void Convection_CPP::_hllc_rho(const TempSOA& rm, const TempSOA& rp,
							   const TempSOA& vm, const TempSOA& vp,
                               const TempSOA& pm, const TempSOA& pp,
                               const TempSOA& Gm, const TempSOA& Gp,
							   TempSOA& out) //? FLOP
{
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix++)
		{
            const Real rho_minus = rm(ix, iy);
			const Real rho_plus  = rp(ix, iy);
            
            const Real v_minus = vm(ix, iy);
			const Real v_plus  = vp(ix, iy);
            
            const Real p_minus = pm(ix, iy);
			const Real p_plus  = pp(ix, iy);
            
            const Real G_minus = Gm(ix, iy);
			const Real G_plus  = Gp(ix, iy);
            
            const Real a_minus = std::sqrt((1/Gm(ix,iy)+1)* max((pm(ix, iy))*((Real)1/rm(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
			const Real a_plus  = std::sqrt((1/Gp(ix,iy)+1)* max((pp(ix, iy))*((Real)1/rp(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
            
            const Real rho_hat = 0.5*(rho_minus+rho_plus);
            const Real a_hat   = 0.5*(a_minus+a_plus);
            const Real rho_hat_a_hat = rho_hat*a_hat;
            
            const Real u_star = 0.5*(v_minus+v_plus+(p_minus-p_plus)/rho_hat_a_hat);
            const Real p_star = 0.5*(p_minus+p_plus)+0.5*(v_minus-v_plus)*rho_hat_a_hat;
            
            const Real q_minus = (p_star <= p_minus)? 1 : sqrt(1+(2*G_minus+1)/(G_minus+1)*(p_star/p_minus-1));//NOT CONSIDERING Pc HERE
            const Real q_plus  = (p_star <= p_plus )? 1 : sqrt(1+(2*G_plus +1)/(G_plus +1)*(p_star/p_plus -1));//NOT CONSIDERING Pc HERE
            
            const Real s_minus = v_minus - a_minus*q_minus;
            const Real s_plus  = v_plus  + a_plus *q_plus;
            
            const Real q_star_minus = rho_minus*(s_minus-v_minus)/(s_minus-u_star);
            const Real q_star_plus  = rho_plus *(s_plus -v_plus )/(s_plus -u_star);
            
            const Real fminus = vm(ix, iy)*rm(ix, iy);
			const Real fplus  = vp(ix, iy)*rp(ix, iy);
            
            const Real f_star_minus = fminus + s_minus*(q_star_minus-rho_minus);
            const Real f_star_plus  = fplus  + s_plus *(q_star_plus -rho_plus );
            
            const bool flagminus = s_minus > 0;
			const bool flagplus  = s_plus  < 0;
			const bool flagother_minus = !(flagminus || flagplus) && (u_star>=0 && s_minus<=0 );
            const bool flagother_plus  = !(flagminus || flagplus) && (u_star<0  && s_plus >=0);
            
            out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagplus)*fplus + ((Real)flagother_minus)*f_star_minus + ((Real)flagother_plus)*f_star_plus;
        }
}

void Convection_CPP::_hllc_phi(const TempSOA& phim, const TempSOA& phip,
                               const TempSOA& rm, const TempSOA& rp,
							   const TempSOA& vm, const TempSOA& vp,
                               const TempSOA& pm, const TempSOA& pp,
                               const TempSOA& Gm, const TempSOA& Gp,
							   TempSOA& out) //? FLOP
{
	for(int iy=0; iy<TempSOA::NY; iy++)
		for(int ix=0; ix<TempSOA::NX; ix++)
		{
            const Real phi_minus = phim(ix, iy);
			const Real phi_plus  = phip(ix, iy);
            
            const Real rho_minus = rm(ix, iy);
			const Real rho_plus  = rp(ix, iy);
            
            const Real v_minus = vm(ix, iy);
			const Real v_plus  = vp(ix, iy);
            
            const Real p_minus = pm(ix, iy);
			const Real p_plus  = pp(ix, iy);
            
            const Real G_minus = Gm(ix, iy);
			const Real G_plus  = Gp(ix, iy);
            
            const Real a_minus = std::sqrt((1/Gm(ix,iy)+1)* max((pm(ix, iy))*((Real)1/rm(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
			const Real a_plus  = std::sqrt((1/Gp(ix,iy)+1)* max((pp(ix, iy))*((Real)1/rp(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
            
            const Real rho_hat = 0.5*(rho_minus+rho_plus);
            const Real a_hat   = 0.5*(a_minus+a_plus);
            const Real rho_hat_a_hat = rho_hat*a_hat;
            
            const Real u_star = 0.5*(v_minus+v_plus+(p_minus-p_plus)/rho_hat_a_hat);
            const Real p_star = 0.5*(p_minus+p_plus)+0.5*(v_minus-v_plus)*rho_hat_a_hat;
            
            const Real q_minus = (p_star <= p_minus)? 1 : sqrt(1+(2*G_minus+1)/(G_minus+1)*(p_star/p_minus-1));//NOT CONSIDERING Pc HERE
            const Real q_plus  = (p_star <= p_plus )? 1 : sqrt(1+(2*G_plus +1)/(G_plus +1)*(p_star/p_plus -1));//NOT CONSIDERING Pc HERE
            
            const Real s_minus = v_minus - a_minus*q_minus;
            const Real s_plus  = v_plus  + a_plus *q_plus;
            
            const Real q_star_minus = phi_minus*(s_minus-v_minus)/(s_minus-u_star);
            const Real q_star_plus  = phi_plus *(s_plus -v_plus )/(s_plus -u_star);
            
            const Real fminus = v_minus*phi_minus;
			const Real fplus  = v_plus *phi_plus;
            
            const Real f_star_minus = fminus + s_minus*(q_star_minus-phi_minus);
            const Real f_star_plus  = fplus  + s_plus *(q_star_plus -phi_plus );
            
            const bool flagminus = s_minus > 0;
			const bool flagplus  = s_plus  < 0;
			const bool flagother_minus = !(flagminus || flagplus) && (u_star>=0 && s_minus<=0 );
            const bool flagother_plus  = !(flagminus || flagplus) && (u_star<0  && s_plus >=0);
            
            out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagplus)*fplus + ((Real)flagother_minus)*f_star_minus + ((Real)flagother_plus)*f_star_plus;
        }
}

void Convection_CPP::_hlle_vel(const TempSOA& rm, const TempSOA& rp,
                               const TempSOA& vm, const TempSOA& vp,
                               const TempSOA& vdm, const TempSOA& vdp,
                               const TempSOA& am, const TempSOA& ap,
                               TempSOA& out) //19 FLOP
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

void Convection_CPP::_hllc_vel(const TempSOA& rm, const TempSOA& rp,
                               const TempSOA& vm, const TempSOA& vp,
                               const TempSOA& vdm, const TempSOA& vdp,
                               const TempSOA& pm, const TempSOA& pp,
                               const TempSOA& Gm, const TempSOA& Gp,
                               TempSOA& out) //? FLOP
{
    for(int iy=0; iy<TempSOA::NY; iy++)
        for(int ix=0; ix<TempSOA::NX; ix++)
        {
            const Real rho_minus = rm(ix, iy);
			const Real rho_plus  = rp(ix, iy);
            
            const Real v_minus = vm(ix, iy);
			const Real v_plus  = vp(ix, iy);
            
            const Real vd_minus = vdm(ix, iy);
			const Real vd_plus  = vdp(ix, iy);
            
            const Real p_minus = pm(ix, iy);
			const Real p_plus  = pp(ix, iy);
            
            const Real G_minus = Gm(ix, iy);
			const Real G_plus  = Gp(ix, iy);
            
            const Real a_minus = std::sqrt((1/Gm(ix,iy)+1)* max((pm(ix, iy))*((Real)1/rm(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
			const Real a_plus  = std::sqrt((1/Gp(ix,iy)+1)* max((pp(ix, iy))*((Real)1/rp(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
            
            const Real rho_hat = 0.5*(rho_minus+rho_plus);
            const Real a_hat   = 0.5*(a_minus+a_plus);
            const Real rho_hat_a_hat = rho_hat*a_hat;
            
            const Real u_star = 0.5*(v_minus+v_plus+(p_minus-p_plus)/rho_hat_a_hat);
            const Real p_star = 0.5*(p_minus+p_plus)+0.5*(v_minus-v_plus)*rho_hat_a_hat;
            
            const Real q_minus = (p_star <= p_minus)? 1 : sqrt(1+(2*G_minus+1)/(G_minus+1)*(p_star/p_minus-1));//NOT CONSIDERING Pc HERE
            const Real q_plus  = (p_star <= p_plus )? 1 : sqrt(1+(2*G_plus +1)/(G_plus +1)*(p_star/p_plus -1));//NOT CONSIDERING Pc HERE
            
            const Real s_minus = v_minus - a_minus*q_minus;
            const Real s_plus  = v_plus  + a_plus *q_plus;
            
            const Real q_star_minus = rho_minus*vd_minus*(s_minus-v_minus)/(s_minus-u_star);
            const Real q_star_plus  = rho_plus *vd_plus *(s_plus -v_plus )/(s_plus -u_star);
            
            const Real fminus = v_minus*vd_minus*rho_minus;
			const Real fplus  = v_plus *vd_plus *rho_plus ;
            
            const Real f_star_minus = fminus + s_minus*(q_star_minus-rho_minus*vd_minus);
            const Real f_star_plus  = fplus  + s_plus *(q_star_plus -rho_plus *vd_plus );
            
            const bool flagminus = s_minus > 0;
			const bool flagplus  = s_plus  < 0;
			const bool flagother_minus = !(flagminus || flagplus) && (u_star>=0 && s_minus<=0 );
            const bool flagother_plus  = !(flagminus || flagplus) && (u_star<0  && s_plus >=0);
            
            out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagplus)*fplus + ((Real)flagother_minus)*f_star_minus + ((Real)flagother_plus)*f_star_plus;
        }
}

void Convection_CPP::_hlle_pvel(const TempSOA& rm, const TempSOA& rp,
                                const TempSOA& vm, const TempSOA& vp,
                                const TempSOA& pm, const TempSOA& pp,
                                const TempSOA& am, const TempSOA& ap,
                                TempSOA& out) //21 FLOP
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

void Convection_CPP::_hllc_pvel(const TempSOA& rm, const TempSOA& rp,
                                const TempSOA& vm, const TempSOA& vp,
                                const TempSOA& pm, const TempSOA& pp,
                                const TempSOA& Gm, const TempSOA& Gp,
                                TempSOA& out) //? FLOP
{
    for(int iy=0; iy<TempSOA::NY; iy++)
        for(int ix=0; ix<TempSOA::NX; ix++)
        {
            const Real rho_minus = rm(ix, iy);
			const Real rho_plus  = rp(ix, iy);
            
            const Real v_minus = vm(ix, iy);
			const Real v_plus  = vp(ix, iy);
            
            const Real p_minus = pm(ix, iy);
			const Real p_plus  = pp(ix, iy);
            
            const Real G_minus = Gm(ix, iy);
			const Real G_plus  = Gp(ix, iy);
            
            const Real a_minus = std::sqrt((1/Gm(ix,iy)+1)* max((pm(ix, iy))*((Real)1/rm(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
			const Real a_plus  = std::sqrt((1/Gp(ix,iy)+1)* max((pp(ix, iy))*((Real)1/rp(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
            
            const Real rho_hat = 0.5*(rho_minus+rho_plus);
            const Real a_hat   = 0.5*(a_minus+a_plus);
            const Real rho_hat_a_hat = rho_hat*a_hat;
            
            const Real u_star = 0.5*(v_minus+v_plus+(p_minus-p_plus)/rho_hat_a_hat);
            const Real p_star = 0.5*(p_minus+p_plus+(v_minus-v_plus)*rho_hat_a_hat);
            
            const Real q_minus = (p_star <= p_minus)? 1 : sqrt(1+(2*G_minus+1)/(G_minus+1)*(p_star/p_minus-1));//NOT CONSIDERING Pc HERE
            const Real q_plus  = (p_star <= p_plus )? 1 : sqrt(1+(2*G_plus +1)/(G_plus +1)*(p_star/p_plus -1));//NOT CONSIDERING Pc HERE
            
            const Real s_minus = v_minus - a_minus*q_minus;
            const Real s_plus  = v_plus  + a_plus *q_plus;
            
            const Real q_star_minus = rho_minus*u_star*(s_minus-v_minus)/(s_minus-u_star);
            const Real q_star_plus  = rho_plus *u_star*(s_plus -v_plus )/(s_plus -u_star);
            
            const Real fminus = v_minus*v_minus*rho_minus + p_minus;
			const Real fplus  = v_plus *v_plus *rho_plus  + p_plus;
            
            const Real f_star_minus = fminus + s_minus*(q_star_minus-rho_minus*v_minus);
            const Real f_star_plus  = fplus  + s_plus *(q_star_plus -rho_plus *v_plus );
            
            const bool flagminus = s_minus > 0;
			const bool flagplus  = s_plus  < 0;
			const bool flagother_minus = !(flagminus || flagplus) && (u_star>=0 && s_minus<=0 );
            const bool flagother_plus  = !(flagminus || flagplus) && (u_star<0  && s_plus >=0);
            
            out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagplus)*fplus + ((Real)flagother_minus)*f_star_minus + ((Real)flagother_plus)*f_star_plus;
        }
}

void Convection_CPP::_hlle_e(const TempSOA& rm, const TempSOA& rp,
                             const TempSOA& vdm, const TempSOA& vdp,
                             const TempSOA& v1m, const TempSOA& v1p,
                             const TempSOA& v2m, const TempSOA& v2p,
                             const TempSOA& pm, const TempSOA& pp,
                             const TempSOA& Gm, const TempSOA& Gp,
                             const TempSOA& Pm, const TempSOA& Pp,
                             const TempSOA& am, const TempSOA& ap,
                             TempSOA& out) //73 FLOP
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
            const Real eminus = pminus*Gm(ix,iy) +
            ((Real)0.5)*rm(ix, iy)*(vdminus*vdminus + v1minus*v1minus + v2minus*v2minus) + Pm(ix,iy);
            
            const Real vdplus = vdp(ix, iy);
            const Real v1plus = v1p(ix, iy);
            const Real v2plus = v2p(ix, iy);
            const Real pplus = pp(ix, iy);
            
            const Real eplus = pplus*Gp(ix,iy) +
            ((Real)0.5)*rp(ix, iy)*(vdplus*vdplus + v1plus*v1plus + v2plus*v2plus) + Pp(ix,iy);
            
            const Real fminus = vdminus*(pminus + eminus);
            const Real fpluss = vdplus *(pplus + eplus);
            
            const Real aminus = am(ix, iy);
            const Real aplus  = ap(ix, iy);
            
            const Real fother = (aplus*fminus-aminus*fpluss+aminus*aplus*(eplus-eminus))*((Real)1/(aplus-aminus));
            
            out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagplus)*fpluss + ((Real)flagother)*fother;
        }
}

void Convection_CPP::_hllc_e(const TempSOA& rm, const TempSOA& rp,
                             const TempSOA& vdm, const TempSOA& vdp,
                             const TempSOA& v1m, const TempSOA& v1p,
                             const TempSOA& v2m, const TempSOA& v2p,
                             const TempSOA& pm, const TempSOA& pp,
                             const TempSOA& Gm, const TempSOA& Gp,
                             const TempSOA& Pm, const TempSOA& Pp,
                             TempSOA& out) //? FLOP
{
    for(int iy=0; iy<TempSOA::NY; iy++)
        for(int ix=0; ix<TempSOA::NX; ix++)
        {
            const Real rho_minus = rm(ix, iy);
			const Real rho_plus  = rp(ix, iy);
            
            const Real v_minus = vdm(ix, iy);
			const Real v_plus  = vdp(ix, iy);
            
            const Real v1_minus = v1m(ix, iy);
			const Real v1_plus  = v1p(ix, iy);
            
            const Real v2_minus = v2m(ix, iy);
			const Real v2_plus  = v2p(ix, iy);
            
            const Real p_minus = pm(ix, iy);
			const Real p_plus  = pp(ix, iy);
            
            const Real G_minus = Gm(ix, iy);
			const Real G_plus  = Gp(ix, iy);
            
            const Real P_minus = Pm(ix, iy);
			const Real P_plus  = Pp(ix, iy);
            
            const Real a_minus = std::sqrt((1/Gm(ix,iy)+1)* max((pm(ix, iy))*((Real)1/rm(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
			const Real a_plus  = std::sqrt((1/Gp(ix,iy)+1)* max((pp(ix, iy))*((Real)1/rp(ix, iy)), (Real)0));//NOT CONSIDERING Pc HERE
            
            const Real rho_hat = 0.5*(rho_minus+rho_plus);
            const Real a_hat   = 0.5*(a_minus+a_plus);
            const Real rho_hat_a_hat = rho_hat*a_hat;
            
            const Real u_star = 0.5*(v_minus+v_plus+(p_minus-p_plus)/rho_hat_a_hat);
            const Real p_star = 0.5*(p_minus+p_plus)+0.5*(v_minus-v_plus)*rho_hat_a_hat;
            
            const Real q_minus = (p_star <= p_minus)? 1 : sqrt(1+(2*G_minus+1)/(G_minus+1)*(p_star/p_minus-1));//NOT CONSIDERING Pc HERE
            const Real q_plus  = (p_star <= p_plus )? 1 : sqrt(1+(2*G_plus +1)/(G_plus +1)*(p_star/p_plus -1));//NOT CONSIDERING Pc HERE
            
            const Real s_minus = v_minus - a_minus*q_minus;
            const Real s_plus  = v_plus  + a_plus *q_plus;
            
            const Real e_minus = p_minus*G_minus + ((Real)0.5)*rho_minus*(v_minus*v_minus + v1_minus*v1_minus + v2_minus*v2_minus) + P_minus;
            const Real e_plus  = p_plus *G_plus  + ((Real)0.5)*rho_plus *(v_plus *v_plus  + v1_plus *v1_plus  + v2_plus *v2_plus)  + P_plus;
            
            const Real q_star_minus = (s_minus-v_minus)/(s_minus-u_star)*(e_minus + (u_star-v_minus)*(rho_minus*u_star+p_minus/(s_minus-v_minus)));
            const Real q_star_plus  = (s_plus -v_plus )/(s_plus -u_star)*(e_plus  + (u_star-v_plus )*(rho_plus *u_star+p_plus /(s_plus -v_plus )));
            
            const Real fminus = v_minus*(e_minus + p_minus);
			const Real fplus  = v_plus *(e_plus  + p_plus );
            
            const Real f_star_minus = fminus + s_minus*(q_star_minus- e_minus);
            const Real f_star_plus  = fplus  + s_plus *(q_star_plus - e_plus );
            
            const bool flagminus = s_minus > 0;
			const bool flagplus  = s_plus  < 0;
			const bool flagother_minus = !(flagminus || flagplus) && (u_star>=0 && s_minus<=0 );
            const bool flagother_plus  = !(flagminus || flagplus) && (u_star<0  && s_plus >=0);
            
            out.ref(ix, iy) = ((Real)flagminus)*fminus + ((Real)flagplus)*fplus + ((Real)flagother_minus)*f_star_minus + ((Real)flagother_plus)*f_star_plus;
        }
}

void Convection_CPP::_xdivergence(const TempSOA& flux, OutputSOA& rhs)
{
    const Real * const f = flux.ptr(0,0);
    Real * const r = &rhs.ref(0,0);
    
    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
            r[ix + OutputSOA::PITCH*iy] = f[ix + 1 + TempSOA::PITCH*iy] - f[ix + TempSOA::PITCH*iy];
}

void Convection_CPP::_ydivergence(const TempSOA& flux, OutputSOA& rhs)
{
    const Real * const f = flux.ptr(0,0);
    Real * const r = &rhs.ref(0,0);
    
    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
            r[ix + OutputSOA::PITCH*iy] += f[iy +  1 + TempSOA::PITCH*ix] - f[iy + TempSOA::PITCH*ix];
}

void Convection_CPP::_zdivergence(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs)
{
    const Real * const ff = fforward.ptr(0,0);
    const Real * const fb = fback.ptr(0,0);
    
    Real * const r = &rhs.ref(0,0);
    
    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
            r[ix + OutputSOA::PITCH*iy] += ff[ix +  TempSOA::PITCH*iy] - fb[ix + TempSOA::PITCH*iy];
}

void Convection_CPP::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    const Real factor2 = ((Real)1.)/6;
    
    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
        {
            AssumedType& rhs = *(AssumedType*)(gptfirst + gptfloats*(ix + iy*rowgpts));
            
            rhs.r = a*rhs.r - dtinvh*rho.rhs(ix, iy);
            rhs.u = a*rhs.u - dtinvh*u.rhs(ix, iy);
            rhs.v = a*rhs.v - dtinvh*v.rhs(ix, iy);
            rhs.w = a*rhs.w - dtinvh*w.rhs(ix, iy);
            rhs.s = a*rhs.s - dtinvh*p.rhs(ix, iy);
            rhs.G = a*rhs.G - dtinvh*(G.rhs(ix, iy) - divu(ix,iy)*sumG(ix,iy)*factor2);
            rhs.P = 0;//a*rhs.P - dtinvh*(P.rhs(ix, iy) - divu(ix,iy)*sumP(ix,iy)*factor2);
        }
}

void Convection_CPP::_xflux(const int relid)
{
    _xweno_minus(rho.ring(relid), rho.weno.ref(0));
    _xweno_pluss(rho.ring(relid), rho.weno.ref(1));
    _xweno_minus(u.ring(relid), u.weno.ref(0));
    _xweno_pluss(u.ring(relid), u.weno.ref(1));
    _xweno_minus(v.ring(relid), v.weno.ref(0));
    _xweno_pluss(v.ring(relid), v.weno.ref(1));
    _xweno_minus(w.ring(relid), w.weno.ref(0));
    _xweno_pluss(w.ring(relid), w.weno.ref(1));
    _xweno_minus(p.ring(relid), p.weno.ref(0));
    _xweno_pluss(p.ring(relid), p.weno.ref(1));
    _xweno_minus(G.ring(relid), G.weno.ref(0));
    _xweno_pluss(G.ring(relid), G.weno.ref(1));
    _xweno_minus(P.ring(relid), P.weno.ref(0));
    _xweno_pluss(P.ring(relid), P.weno.ref(1));
    _char_vel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1));
    //_xextraterm(u.weno(0), u.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1));
    _xextraterm_hllc(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1));
    
    //_hlle_rho(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), rho.flux.ref());
    _hllc_rho(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), rho.flux.ref());
    //_hlle_pvel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), u.flux.ref());
    _hllc_pvel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), u.flux.ref());
    //_hlle_vel(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), v.flux.ref());
    _hllc_vel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), v.flux.ref());
    //_hlle_vel(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), w.flux.ref());
    _hllc_vel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), w.flux.ref());
    //_hlle_e(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1), p.flux.ref());
    _hllc_e(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), p.flux.ref());
    
    //_hlle_rho(G.weno(0), G.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), G.flux.ref());
    _hllc_phi(G.weno(0), G.weno(1), rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), G.flux.ref());
    //_hlle_rho(P.weno(0), P.weno(1), u.weno(0), u.weno(1), charvel(0), charvel(1), P.flux.ref());
    //_hllc_phi(P.weno(0), P.weno(1), rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.flux.ref());
}

void Convection_CPP::_yflux(const int relid)
{
    _yweno_minus(rho.ring(relid), rho.weno.ref(0));
    _yweno_pluss(rho.ring(relid), rho.weno.ref(1));
    _yweno_minus(u.ring(relid), u.weno.ref(0));
    _yweno_pluss(u.ring(relid), u.weno.ref(1));
    _yweno_minus(v.ring(relid), v.weno.ref(0));
    _yweno_pluss(v.ring(relid), v.weno.ref(1));
    _yweno_minus(w.ring(relid), w.weno.ref(0));
    _yweno_pluss(w.ring(relid), w.weno.ref(1));
    _yweno_minus(p.ring(relid), p.weno.ref(0));
    _yweno_pluss(p.ring(relid), p.weno.ref(1));
    _yweno_minus(G.ring(relid), G.weno.ref(0));
    _yweno_pluss(G.ring(relid), G.weno.ref(1));
    _yweno_minus(P.ring(relid), P.weno.ref(0));
    _yweno_pluss(P.ring(relid), P.weno.ref(1));
    _char_vel(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1));
    _yextraterm(v.weno(0), v.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1));
    
    _hlle_rho(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), rho.flux.ref());
    _hlle_vel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), u.flux.ref());
    _hlle_pvel(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), v.flux.ref());
    _hlle_vel(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), w.flux.ref());
    _hlle_e(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), u.weno(0), u.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1), p.flux.ref());
    
    _hlle_rho(G.weno(0), G.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), G.flux.ref());
    _hlle_rho(P.weno(0), P.weno(1), v.weno(0), v.weno(1), charvel(0), charvel(1), P.flux.ref());
}

void Convection_CPP::_zflux(const int relid)
{
    _zweno_minus(relid, rho.ring, rho.weno.ref(0));
    _zweno_pluss(relid, rho.ring, rho.weno.ref(1));
    _zweno_minus(relid, u.ring, u.weno.ref(0));
    _zweno_pluss(relid, u.ring, u.weno.ref(1));
    _zweno_minus(relid, v.ring, v.weno.ref(0));
    _zweno_pluss(relid, v.ring, v.weno.ref(1));
    _zweno_minus(relid, w.ring, w.weno.ref(0));
    _zweno_pluss(relid, w.ring, w.weno.ref(1));
    _zweno_minus(relid, p.ring, p.weno.ref(0));
    _zweno_pluss(relid, p.ring, p.weno.ref(1));
    _zweno_minus(relid, G.ring, G.weno.ref(0));
    _zweno_pluss(relid, G.ring, G.weno.ref(1));
    
    _zweno_minus(relid, P.ring, P.weno.ref(0));
    _zweno_pluss(relid, P.ring, P.weno.ref(1));
    _char_vel(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel.ref(0), charvel.ref(1));
    _zextraterm(w.weno(-2), w.weno(-1),w.weno(0), w.weno(1), G.weno(0), G.weno(-1), P.weno(0), P.weno(-1), charvel(-2), charvel(-1), charvel(0), charvel(1));
    
    _hlle_rho(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), rho.flux.ref());
    _hlle_vel(rho.weno(0), rho.weno(1), u.weno(0), u.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), u.flux.ref());
    _hlle_vel(rho.weno(0), rho.weno(1), v.weno(0), v.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), v.flux.ref());
    _hlle_pvel(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), p.weno(0), p.weno(1), charvel(0), charvel(1), w.flux.ref());
    _hlle_e(rho.weno(0), rho.weno(1), w.weno(0), w.weno(1), u.weno(0), u.weno(1), v.weno(0), v.weno(1), p.weno(0), p.weno(1), G.weno(0), G.weno(1), P.weno(0), P.weno(1), charvel(0), charvel(1), p.flux.ref());
    
    _hlle_rho(G.weno(0), G.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), G.flux.ref());
    _hlle_rho(P.weno(0), P.weno(1), w.weno(0), w.weno(1), charvel(0), charvel(1), P.flux.ref());
}

void Convection_CPP::_xrhs()
{
    _xdivergence(rho.flux(), rho.rhs);
    _xdivergence(u.flux(), u.rhs);
    _xdivergence(v.flux(), v.rhs);
    _xdivergence(w.flux(), w.rhs);
    _xdivergence(p.flux(), p.rhs);
    _xdivergence(G.flux(), G.rhs);
    _xdivergence(P.flux(), P.rhs);
}

void Convection_CPP::_yrhs()
{
    _ydivergence(rho.flux(), rho.rhs);
    _ydivergence(u.flux(), u.rhs);
    _ydivergence(v.flux(), v.rhs);
    _ydivergence(w.flux(), w.rhs);
    _ydivergence(p.flux(), p.rhs);
    _ydivergence(G.flux(), G.rhs);
    _ydivergence(P.flux(), P.rhs);
}

void Convection_CPP::_zrhs()
{
    _zdivergence(rho.flux(-1), rho.flux(0), rho.rhs);
    _zdivergence(u.flux(-1), u.flux(0), u.rhs);
    _zdivergence(v.flux(-1), v.flux(0), v.rhs);
    _zdivergence(w.flux(-1), w.flux(0), w.rhs);
    _zdivergence(p.flux(-1), p.flux(0), p.rhs);
    _zdivergence(G.flux(-1), G.flux(0), G.rhs);
    _zdivergence(P.flux(-1), P.flux(0), P.rhs);
}
