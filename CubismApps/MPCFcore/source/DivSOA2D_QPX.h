#pragma once

#include "common.h"

class DivSOA2D_QPX
{
	public:
		void xrhs(const TempSOA& flux, OutputSOA& rhs) const;
		void yrhs(const TempSOA& flux, OutputSOA& rhs) const;
		void zrhs(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs) const;

		void xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& Gm, const TempSOA& G, 
				const TempSOA& Pm, const TempSOA& Pp,
				const TempSOA& am, const TempSOA& ap,
				OutputSOA& divu, OutputSOA& sumG, OutputSOA& sumP);
	
		void yextraterm(const TempSOA& um, const TempSOA& up, 
				const TempSOA& Gm, const TempSOA& Gp,
				const TempSOA& Pm, const TempSOA& Pp,
				const TempSOA& am, const TempSOA& ap,
				OutputSOA& divu, OutputSOA& sumG, OutputSOA& sumP);
	
		void zextraterm(const TempSOA& um0, const TempSOA& up0, const TempSOA& um1, const TempSOA& up1, 
				const TempSOA& Gm, const TempSOA& Gp, const TempSOA& Pm, const TempSOA& Pp, 
				const TempSOA& am0, const TempSOA& ap0, const TempSOA& am1, const TempSOA& ap1,
				OutputSOA& divu, OutputSOA& sumG, OutputSOA& sumP);
};
