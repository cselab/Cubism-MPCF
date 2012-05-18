/*
 *  Convection_CPP.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/6/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "common.h"
#include "SOA2D.h"

using namespace std;

class Convection_CPP
{	
public:
	
	const Real a, dtinvh; //LSRK3-related "a" factor, and "lambda"
	const Real gamma1, gamma2, smoothlength; //per-phase specific heat ratios, and smoothing length
    const Real pc1, pc2; //per-phase pressure correction terms, used only for liquids
	
	//the only constructor for this class
	Convection_CPP(const Real a, const Real dtinvh, 
				   const Real gamma1, const Real gamma2, const Real smoothlength, 
				   const Real pc1, const Real pc2);
	
	//main method of the class, it evaluates the convection term of the RHS
	void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
				 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts);
	
	//this provides the amount of flops and memory traffic performed in compute(.)
	void hpc_info(float& flop_convert, int& traffic_convert,
				  float& flop_weno, int& traffic_weno,
				  float& flop_extraterm, int& traffic_extraterm, 
				  float& flop_charvel, int& traffic_charvel,
				  float& flop_hlle, int& traffic_hlle,
				  float& flop_div, int& traffic_div,
				  float& flop_copyback, int& traffic_copyback,
				  size_t& footprint) const;
	
	//this report the performance details of compute(.) given the measured time
	void printflops(const float PEAKPERF_CORE, const float PEAKBAND, 
					const int NCORES, const int NTIMES, const int NBLOCKS, 
					const float MEASUREDTIME) const;
	
protected:
	
	//Assumed input/output grid point type
	struct AssumedType { Real r, u, v, w, s, l; }; 
	
	//working dataset types
	typedef SOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3> InputSOA; //associated with weno5
	typedef RingSOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3, 6> RingInputSOA; //associated with weno5
	typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_> TempSOA; //associated with hlle
	typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_, 2> RingTempSOA; //associated with hlle
	typedef RingSOA2D<0,_BLOCKSIZE_+1, 0, _BLOCKSIZE_, 3> RingTempSOA3;//associated with hlle
	
	//working dataset
	RingInputSOA ringrho, ringu, ringv, ringw, ringp, ringls;
	RingTempSOA wenorho, wenou, wenov, wenop, charvel;
	RingTempSOA3 wenow, wenols;
	RingTempSOA fluxrho, fluxu, fluxv, fluxw, fluxp, fluxls;
	OutputSOA sumls, divu;
	OutputSOA rhsrho, rhsu, rhsv, rhsw, rhss, rhsls;
	
	void _next() 
	{ 
		ringrho.next(); ringu.next(); ringv.next(); ringw.next(); ringp.next(); ringls.next();
	}
	
	void _flux_next()
	{
		fluxrho.next(); fluxu.next(); fluxv.next(); fluxw.next(); fluxp.next(); fluxls.next();
		
		wenow.next(); wenols.next();
	}
	
	inline Real _getgamma(const Real phi) const { return reconstruct(gamma1, gamma2, phi, 1/smoothlength); } 
	inline Real _getPC(const Real phi) const { return reconstruct(pc1, pc2, phi, 1/smoothlength); }
    
	virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
	
	virtual void _xweno_minus(const InputSOA& in, TempSOA& out);
	virtual void _xweno_pluss(const InputSOA& in, TempSOA& out);
	virtual void _yweno_minus(const InputSOA& in, TempSOA& out);
	virtual void _yweno_pluss(const InputSOA& in, TempSOA& out);
	virtual void _zweno_minus(const int relid, const RingInputSOA& in, TempSOA& out);
	virtual void _zweno_pluss(const int relid, const RingInputSOA& in, TempSOA& out);
	
	virtual void _xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp);
	virtual void _yextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp);
	virtual void _zextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp);
	
	virtual void _char_vel(const TempSOA& rminus, const TempSOA& rplus, 
						   const TempSOA& vminus, const TempSOA& vplus,
						   const TempSOA& pminus, const TempSOA& pplus,
						   const TempSOA& lminus, const TempSOA& lplus, 
						   TempSOA& out_minus, TempSOA& out_plus);
	
	virtual void _hlle_rho(const TempSOA& rm, const TempSOA& rp,
						   const TempSOA& vm, const TempSOA& vp,
						   const TempSOA& am, const TempSOA& ap,
						   TempSOA& out);
	
	virtual void _hlle_vel(const TempSOA& rminus, const TempSOA& rplus,
						   const TempSOA& vminus, const TempSOA& vplus,
						   const TempSOA& vdminus, const TempSOA& vdplus,
						   const TempSOA& aminus, const TempSOA& aplus,
						   TempSOA& out);
	
	virtual void _hlle_pvel(const TempSOA& rminus, const TempSOA& rplus,
							const TempSOA& vminus, const TempSOA& vplus,
							const TempSOA& pminus, const TempSOA& pplus,
							const TempSOA& aminus, const TempSOA& aplus,
							TempSOA& out);
	virtual void _hlle_e(const TempSOA& rminus, const TempSOA& rplus,
						 const TempSOA& vdminus, const TempSOA& vdplus,
						 const TempSOA& v1minus, const TempSOA& v1plus,
						 const TempSOA& v2minus, const TempSOA& v2plus,
						 const TempSOA& pminus, const TempSOA& pplus,
						 const TempSOA& lsminus, const TempSOA& lsplus, 
						 const TempSOA& aminus, const TempSOA& aplus,
						 TempSOA& out);
	
	virtual void _xdivergence(const TempSOA& flux, OutputSOA& rhs);
	virtual void _ydivergence(const TempSOA& flux, OutputSOA& rhs);
	virtual void _zdivergence(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs);
	
	virtual void _xflux(const int relsliceid);
	virtual void _yflux(const int relsliceid);
	virtual void _zflux(const int relsliceid);
	
	virtual void _xrhs();	
	virtual void _yrhs();
	virtual void _zrhs();
	
	virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);
};
