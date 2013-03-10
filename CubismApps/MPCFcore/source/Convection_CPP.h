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
	
	//the only constructor for this class
	Convection_CPP(const Real a, const Real dtinvh);
	
	//main method of the class, it evaluates the convection term of the RHS
	void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
				 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts);
	
	//this provides the amount of flops and memory traffic performed in compute(.)
	static void hpc_info(float& flop_convert, int& traffic_convert,
						 float& flop_weno, int& traffic_weno,
						 float& flop_extraterm, int& traffic_extraterm,
						 float& flop_charvel, int& traffic_charvel,
						 float& flop_hlle, int& traffic_hlle,
						 float& flop_div, int& traffic_div,
						 float& flop_copyback, int& traffic_copyback,
						 size_t& footprint);
	
	//this report the performance details of compute(.) given the measured time
	static void printflops(const float PEAKPERF_CORE, const float PEAKBAND,
						   const int NCORES, const int NTIMES, const int NBLOCKS,
						   const float MEASUREDTIME);
	
protected:
	
	//Assumed input/output grid point type
	struct AssumedType { Real r, u, v, w, s, G, P; };
	
	template<int zslices=2>
	struct WorkingSet {
		RingInputSOA ring;
		RingSOA2D<0,_BLOCKSIZE_+1, 0, _BLOCKSIZE_, zslices> weno;
		RingTempSOA flux;
		OutputSOA rhs;
	};
	
	WorkingSet<2> rho, u, v, p;
	WorkingSet<4> w, G;
	WorkingSet<4> P;
	
	//used for HLLE fluxes
	RingSOA2D<0,_BLOCKSIZE_+1, 0, _BLOCKSIZE_, 4> charvel;
	
    //for hllc
    RingSOA2D<0,_BLOCKSIZE_+1, 0, _BLOCKSIZE_, 2> star;
    
	//used for the extra term
	OutputSOA sumG, divu;
	OutputSOA sumP;
	
	void _next()
	{
		rho.ring.next(); u.ring.next(); v.ring.next(); w.ring.next(); p.ring.next(); G.ring.next();
		P.ring.next();
	}
	
	void _flux_next()
	{
		rho.flux.next(); u.flux.next(); v.flux.next(); w.flux.next(); p.flux.next(); G.flux.next();
		
		w.weno.next(); G.weno.next(); w.weno.next(); G.weno.next(); charvel.next(); charvel.next();
		
		P.flux.next();
		P.weno.next();
        P.weno.next();
        
        star.next();
	}
    	
	void _xweno_minus(const InputSOA& in, TempSOA& out);
	void _xweno_pluss(const InputSOA& in, TempSOA& out);
	void _yweno_minus(const InputSOA& in, TempSOA& out);
	void _yweno_pluss(const InputSOA& in, TempSOA& out);
	void _zweno_minus(const int relid, const RingInputSOA& in, TempSOA& out);
	void _zweno_pluss(const int relid, const RingInputSOA& in, TempSOA& out);
	
	void _xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& Gp
							 , const TempSOA& Pm, const TempSOA& Pp
                             , const TempSOA& am, const TempSOA& ap);
    
    virtual void _xextraterm_v2(const TempSOA& um, const TempSOA& up, const InputSOA& G, const InputSOA& P, const TempSOA& am, const TempSOA& ap);
    
	virtual void _yextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& Gp
							 , const TempSOA& Pm, const TempSOA& Pp
                             , const TempSOA& am, const TempSOA& ap);
    
    virtual void _yextraterm_v2(const TempSOA& um, const TempSOA& up, const InputSOA& G
                                , const InputSOA& P
                                , const TempSOA& am, const TempSOA& ap);
    
	virtual void _zextraterm(const TempSOA& um0, const TempSOA& up0, const TempSOA& um1, const TempSOA& up1, const TempSOA& Gm, const TempSOA& Gp
							 , const TempSOA& Pm, const TempSOA& Pp
                             , const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1);
	
    virtual void _zextraterm_v2(const TempSOA& um0, const TempSOA& up0, const TempSOA& um1, const TempSOA& up1,
                                const InputSOA& G, const InputSOA& P,
                                const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1, const bool bFirst=false);
    
    template <int SizeDiMerdaX>
	void _char_vel(const TempSOA& rminus, const TempSOA& rplus,
						   const TempSOA& vminus, const TempSOA& vplus,
						   const TempSOA& pminus, const TempSOA& pplus,
						   const TempSOA& Gminus, const TempSOA& Gplus,
						   const TempSOA& Pminus, const TempSOA& Pplus,
						   TempSOA& out_minus, TempSOA& out_plus, const int relid=0);
	
    template <int SizeDiMerdaX>
	void _hlle_rho(const TempSOA& rm, const TempSOA& rp,
						   const TempSOA& vm, const TempSOA& vp,
						   const TempSOA& am, const TempSOA& ap,
						   TempSOA& out);
	
    template <int SizeDiMerdaX>
	void _hlle_vel(const TempSOA& rminus, const TempSOA& rplus,
						   const TempSOA& vminus, const TempSOA& vplus,
						   const TempSOA& vdminus, const TempSOA& vdplus,
						   const TempSOA& aminus, const TempSOA& aplus,
						   TempSOA& out);
	
    template <int SizeDiMerdaX>
	void _hlle_pvel(const TempSOA& rminus, const TempSOA& rplus,
							const TempSOA& vminus, const TempSOA& vplus,
							const TempSOA& pminus, const TempSOA& pplus,
							const TempSOA& aminus, const TempSOA& aplus,
							TempSOA& out);
    
    template <int SizeDiMerdaX>
	void _hlle_e(const TempSOA& rminus, const TempSOA& rplus,
						 const TempSOA& vdminus, const TempSOA& vdplus,
						 const TempSOA& v1minus, const TempSOA& v1plus,
						 const TempSOA& v2minus, const TempSOA& v2plus,
						 const TempSOA& pminus, const TempSOA& pplus,
						 const TempSOA& Gminus, const TempSOA& Gplus,
						 const TempSOA& Pminus, const TempSOA& Pplus,
						 const TempSOA& aminus, const TempSOA& aplus,
						 TempSOA& out);

    virtual Real _u_hllc(const Real v_minus, const Real v_plus, const Real s_minus, const Real s_plus, const Real star);
    
    virtual void _xextraterm_hllc(const TempSOA& vm, const TempSOA& vp
                                  , const TempSOA& Gm, const TempSOA& Gp
                                  , const TempSOA& Pm, const TempSOA& Pp
                                  , const TempSOA& sm, const TempSOA& sp
                                  , const TempSOA& star);
    
    virtual void _xextraterm_hllc_v2(const TempSOA& vm, const TempSOA& vp,
                                     const InputSOA& G, const InputSOA& P,
                                     const TempSOA& sm, const TempSOA& sp,
                                     const TempSOA& star);
    
    virtual void _yextraterm_hllc(const TempSOA& vm, const TempSOA& vp
                                  , const TempSOA& Gm, const TempSOA& Gp
                                  , const TempSOA& Pm, const TempSOA& Pp
                                  , const TempSOA& sm, const TempSOA& sp
                                  , const TempSOA& star);
    
    virtual void _yextraterm_hllc_v2(const TempSOA& vm, const TempSOA& vp,
                                     const InputSOA& G, const InputSOA& P,
                                     const TempSOA& sm, const TempSOA& sp,
                                     const TempSOA& star);
    
    virtual void _zextraterm_hllc(const TempSOA& v0m, const TempSOA& v0p, const TempSOA& v1m, const TempSOA& v1p
                                  , const TempSOA& Gm, const TempSOA& Gp
                                  , const TempSOA& Pm, const TempSOA& Pp
                                  , const TempSOA& s0m, const TempSOA& s0p , const TempSOA& s1m, const TempSOA& s1p
                                  , const TempSOA& star0, const TempSOA& star1);
    
    virtual void _zextraterm_hllc_v2(const TempSOA& v0m, const TempSOA& v0p, const TempSOA& v1m, const TempSOA& v1p,
                                     const InputSOA& G, const InputSOA& P,
                                     const TempSOA& s0m, const TempSOA& s0p , const TempSOA& s1m, const TempSOA& s1p,
                                     const TempSOA& star0, const TempSOA& star1);
    
    template<int SizeDiMerdaX>
    void _char_vel_hllc(const TempSOA& rm, const TempSOA& rp,
                                const TempSOA& vm, const TempSOA& vp,
                                const TempSOA& pm, const TempSOA& pp,
                                const TempSOA& Gm, const TempSOA& Gp,
                                const TempSOA& Pm, const TempSOA& Pp,
                                TempSOA& outm, TempSOA& outp, TempSOA& out_star);
    
    template<int SizeDiMerdaX>
    void _hllc_rho(const TempSOA& rm, const TempSOA& rp,
                           const TempSOA& vm, const TempSOA& vp,
                           const TempSOA& sm, const TempSOA& sp,
                           const TempSOA& star,
                           TempSOA& out);
    
    template<int SizeDiMerdaX>
    void _hllc_phi(const TempSOA& phim, const TempSOA& phip,
                           const TempSOA& vm, const TempSOA& vp,
                           const TempSOA& sm, const TempSOA& sp,
                           const TempSOA& star,
                           TempSOA& out);
    
    template<int SizeDiMerdaX>
    void _hllc_vel(const TempSOA& rm, const TempSOA& rp,
                           const TempSOA& vm, const TempSOA& vp,
                           const TempSOA& vdm, const TempSOA& vdp,
                           const TempSOA& sm, const TempSOA& sp,
                           const TempSOA& star,
                           TempSOA& out);
    
    template<int SizeDiMerdaX>
    void _hllc_pvel(const TempSOA& rm, const TempSOA& rp,
                            const TempSOA& vm, const TempSOA& vp,
                            const TempSOA& pm, const TempSOA& pp,
                            const TempSOA& sm, const TempSOA& sp,
                            const TempSOA& star,
                            TempSOA& out);
    
    template<int SizeDiMerdaX>
    void _hllc_e(const TempSOA& rm, const TempSOA& rp,
                         const TempSOA& vdm, const TempSOA& vdp,
                         const TempSOA& v1m, const TempSOA& v1p,
                         const TempSOA& v2m, const TempSOA& v2p,
                         const TempSOA& pm, const TempSOA& pp,
                         const TempSOA& Gm, const TempSOA& Gp,
                         const TempSOA& Pm, const TempSOA& Pp,
                         const TempSOA& sm, const TempSOA& sp,
                         const TempSOA& star,
                         TempSOA& out);
    
    /*	virtual void _xdivergence(const TempSOA& flux, OutputSOA& rhs);
	virtual void _ydivergence(const TempSOA& flux, OutputSOA& rhs);
	virtual void _zdivergence(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs);
    */
	void _xdivergence(const TempSOA& flux, OutputSOA& rhs);
	void _ydivergence(const TempSOA& flux, OutputSOA& rhs);
	void _zdivergence(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs);
	
	virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
	
	virtual void _xflux(const int relsliceid);
	virtual void _yflux(const int relsliceid);
	virtual void _zflux(const int relsliceid, const bool bFirst=false);
	
    virtual void _xflux_hllc(const int relsliceid);
    virtual void _yflux_hllc(const int relsliceid);
    virtual void _zflux_hllc(const int relsliceid, const bool bFirst=false);
    
	virtual void _xrhs();
	virtual void _yrhs();
	virtual void _zrhs();
	
	virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);
};

