//
//  SurfaceTension.h
//  MPCFcore
//
//  Created by Babak Hejazialhosseini on 8/2/11.
//  Copyright 2011 ETH Zurich. All rights reserved.
//

#pragma once

#include "DivTensor_CPP.h"

class SurfaceTension_CPP: public virtual DivTensor_CPP
{
	struct AssumedType { Real r, u, v, w, s, l; };
	
public:
	const Real G1, G2;
	const Real smoothing_length;

	SurfaceTension_CPP(const Real a = 1,
					   const Real dtinvh = 1, 
			   const Real G1 = 1,
			   const Real G2 = 1,
					   const Real h = 1,
			   const Real smoothing_length = 1,
					   const Real sigma=1):
	DivTensor_CPP(a, dtinvh, h, sigma), G1(G1), G2(G2), smoothing_length(smoothing_length)
	{ 
	}
	
	//here we will "project" the phases
	void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
	
	void hpc_info(float& flop_convert, int& traffic_convert,
				  float& flop_corners, int& traffic_corner,
				  float& flop_tensorface, int& traffic_tensorface, 
				  float& flop_tdotu, int& traffic_tdotu,
				  float& flop_div, int& traffic_div,
				  float& flop_copyback, int& traffic_copyback,
				  size_t& footprint)
	{
		DivTensor_CPP::hpc_info(flop_convert, traffic_convert, flop_corners, traffic_corner, flop_tensorface, traffic_tensorface, 
								flop_tdotu, traffic_tdotu, flop_div, traffic_div, flop_copyback, traffic_copyback, footprint);
		
		const int ninputs = (int)powf(_BLOCKSIZE_ + 2, 3);
		
		flop_convert = 5 * ninputs;
		traffic_convert = (4 + 4) * sizeof(Real) * ninputs;
		
		footprint = sizeof(*this);
	}
};
