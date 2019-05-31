/*
 *  VectorOperator.cpp
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 03/18/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifdef _NONUNIFORM_BLOCK_
#include "VectorOperatorNonuniform.h"
#else
#include "VectorOperator.h"
#endif /* _NONUNIFORM_BLOCK_ */

// NU = NonUniform
#if defined(_SECONDORDER_FD_CORE_) && defined(_NONUNIFORM_BLOCK_)
#define __FDSCHEME(NAME) #NAME"_2nd_NU"
#elif !defined(_SECONDORDER_FD_CORE_) && defined(_NONUNIFORM_BLOCK_)
#define __FDSCHEME(NAME) #NAME"_4th_NU"
#else
#define __FDSCHEME(NAME) #NAME"_4th"
#endif /* _SECONDORDER_FD_CORE_ */

const std::string OVort::NAME = __FDSCHEME(Vort);
const std::string OIgradA2I::NAME = __FDSCHEME(IgradA2I);
const std::string OdivU::NAME = __FDSCHEME(divU);
const std::string OSSd::NAME = __FDSCHEME(SSDeviatoric);
const std::string OQcrit::NAME = __FDSCHEME(Qcrit);
const std::string OgradU::NAME = __FDSCHEME(gradU);
template <> const std::string OgradU_0::NAME = __FDSCHEME(gradU11);
template <> const std::string OgradU_1::NAME = __FDSCHEME(gradU12);
template <> const std::string OgradU_2::NAME = __FDSCHEME(gradU13);
template <> const std::string OgradU_3::NAME = __FDSCHEME(gradU21);
template <> const std::string OgradU_4::NAME = __FDSCHEME(gradU22);
template <> const std::string OgradU_5::NAME = __FDSCHEME(gradU23);
template <> const std::string OgradU_6::NAME = __FDSCHEME(gradU31);
template <> const std::string OgradU_7::NAME = __FDSCHEME(gradU32);
template <> const std::string OgradU_8::NAME = __FDSCHEME(gradU33);
const std::string OUcon::NAME = __FDSCHEME(Ucontraction);
const std::string OPIC::NAME = __FDSCHEME(PIC);

#undef __FDSCHEME
