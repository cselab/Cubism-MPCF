#pragma once

#include "common.h"

class HLLESOA2D_QPX
{
public:
	void rho(const TempSOA& rm, const TempSOA& rp,
			 const TempSOA& vm, const TempSOA& vp,
			 const TempSOA& am, const TempSOA& ap,
			 TempSOA& out) const;
	
	void vel(const TempSOA& rminus, const TempSOA& rplus,
			 const TempSOA& vminus, const TempSOA& vplus,
			 const TempSOA& vdminus, const TempSOA& vdplus,
			 const TempSOA& aminus, const TempSOA& aplus,
			 TempSOA& out) const;
	
	void pvel(const TempSOA& rminus, const TempSOA& rplus,
			  const TempSOA& vminus, const TempSOA& vplus,
			  const TempSOA& pminus, const TempSOA& pplus,
			  const TempSOA& aminus, const TempSOA& aplus,
			  TempSOA& out) const;
	
	void e(const TempSOA& rminus, const TempSOA& rplus,
		   const TempSOA& vdminus, const TempSOA& vdplus,
		   const TempSOA& v1minus, const TempSOA& v1plus,
		   const TempSOA& v2minus, const TempSOA& v2plus,
		   const TempSOA& pminus, const TempSOA& pplus,
		   const TempSOA& Gminus, const TempSOA& Gplus, 
		   const TempSOA& PIminus, const TempSOA& PIplus,
		   const TempSOA& aminus, const TempSOA& aplus,
		   TempSOA& out) const;
	
	void char_vel(const TempSOA& rminus, const TempSOA& rplus, 
				  const TempSOA& vminus, const TempSOA& vplus,
				  const TempSOA& pminus, const TempSOA& pplus,
				  const TempSOA& Gminus, const TempSOA& Gplus, 
				  const TempSOA& Pminus, const TempSOA& Pplus,
				  TempSOA& out_minus, TempSOA& out_plus) const;
};
