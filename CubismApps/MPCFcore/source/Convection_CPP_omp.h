/*
 *  Convection_CPP.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/6/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include "Convection_CPP.h"

using namespace std;

class Convection_CPP_omp: public Convection_CPP
{	
public:
	
//	const Real a, dtinvh; //LSRK3-related "a" factor, and "lambda"
	
	//the only constructor for this class
	Convection_CPP_omp(const Real a, const Real dtinvh);
//	void Convection_CPP(const Real a, const Real dtinvh);
	
	//main method of the class, it evaluates the convection term of the RHS
	void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
				 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts);
	
protected:

#if 1
	virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
	
	virtual void _xweno_minus(const InputSOA& in, TempSOA& out);
	virtual void _xweno_pluss(const InputSOA& in, TempSOA& out);
	virtual void _yweno_minus(const InputSOA& in, TempSOA& out);
	virtual void _yweno_pluss(const InputSOA& in, TempSOA& out);
	virtual void _zweno_minus(const int relid, const RingInputSOA& in, TempSOA& out);
	virtual void _zweno_pluss(const int relid, const RingInputSOA& in, TempSOA& out);
	
	virtual void _xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& Gp
							 , const TempSOA& Pm, const TempSOA& Pp
                             , const TempSOA& am, const TempSOA& ap);
    
	virtual void _yextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& Gp
							 , const TempSOA& Pm, const TempSOA& Pp
                             , const TempSOA& am, const TempSOA& ap);
    
	virtual void _zextraterm(const TempSOA& um0, const TempSOA& up0, const TempSOA& um1, const TempSOA& up1, const TempSOA& Gm, const TempSOA& Gp
							 , const TempSOA& Pm, const TempSOA& Pp
                             , const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1);
	
	virtual void _char_vel(const TempSOA& rminus, const TempSOA& rplus, 
						   const TempSOA& vminus, const TempSOA& vplus,
						   const TempSOA& pminus, const TempSOA& pplus,
						   const TempSOA& Gminus, const TempSOA& Gplus,
						   const TempSOA& Pminus, const TempSOA& Pplus,						   
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
						 const TempSOA& Gminus, const TempSOA& Gplus, 
						 const TempSOA& Pminus, const TempSOA& Pplus, 						 
						 const TempSOA& aminus, const TempSOA& aplus,
						 TempSOA& out);
	
	virtual void _xdivergence(const TempSOA& flux, OutputSOA& rhs);
	virtual void _ydivergence(const TempSOA& flux, OutputSOA& rhs);
	virtual void _zdivergence(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs);

#if 0	
	virtual void _xflux(const int relsliceid);
	virtual void _yflux(const int relsliceid);
	virtual void _zflux(const int relsliceid);
	
	virtual void _xrhs();	
	virtual void _yrhs();
	virtual void _zrhs();
#endif
	
	virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);
#endif
};
