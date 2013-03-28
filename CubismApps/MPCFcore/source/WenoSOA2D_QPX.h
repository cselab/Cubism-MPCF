/*
 *  *  WenoSOA2D_QPX.h
 *   *  
 *    *
 *     *  Created by Diego Rossinelli on 2/12/13.
 *      *  Copyright 2013 ETH Zurich. All rights reserved.
 *       *
 *        */

#include "common.h"

#include "../../MPCFthread/source/Weno_CPP.h"
#include "../../MPCFthread/source/Weno_QPX.h"

#ifdef _ACCURATEWENO_
typedef WenoQPX_MinusFunctor SoupMinus;
typedef WenoQPX_PlusFunctor SoupPlus;
#else
typedef WenoQPX_MinusPanos SoupMinus;
typedef WenoQPX_PlusPanos SoupPlus;
#endif

class WenoSOA2D_QPX
{
public:
	
	void xcompute(const InputSOA& in, TempSOA& outm, TempSOA& outp) const;	
	void ycompute(const InputSOA& in, TempSOA& outm, TempSOA& outp) const;	
	void zcompute(const int r, const RingInputSOA& in, TempSOA& outm, TempSOA& outp) const;
};

