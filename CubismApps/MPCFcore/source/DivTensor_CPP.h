/*
 *  DivTensor_CPP.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/24/12.
 *  Revision by Fabian Wermelinger on 7/26/2017
 *  Copyright 2012/2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef DIVTENSOR_CPP_H_DENX8LSZ
#define DIVTENSOR_CPP_H_DENX8LSZ

#include <cassert>
#include <cmath>

#include "common.h"
#include "SOA2D.h"

typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_+1, float> TempSOAf_ST;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_+1> TempSOA_ST;
typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_+1, 2> RingTempSOA_ST;

typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_, float> TempPiXSOAf_ST;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_+1, float> TempPiYSOAf_ST;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_, float> TempPiZSOAf_ST;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_> TempPiXSOA_ST;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_+1> TempPiYSOA_ST;
typedef SOA2D<0,_BLOCKSIZE_, 0,_BLOCKSIZE_> TempPiZSOA_ST;
typedef RingSOA2D<0, _BLOCKSIZE_, 0,_BLOCKSIZE_, 2> RingTempPiZSOA_ST;


template <typename TInput>
class DivTensor_CPP
{
public:
    virtual void hpc_info(float& flop_corners, int& traffic_corner,
                          float& flop_tdotu, int& traffic_tdotu,
                          float& flop_div, int& traffic_div,
                          size_t& footprint)
    {
        const int ncorners = (int)std::pow(_BLOCKSIZE_ + 1, 3);
        const int ntensors = 3 * (_BLOCKSIZE_ + 1) * _BLOCKSIZE_ * _BLOCKSIZE_;
        const int ncells = (int)std::pow(_BLOCKSIZE_, 3);

        flop_corners = 21 * ncorners;
        traffic_corner = (8 + 3) * sizeof(Real) * ncorners;
        flop_tdotu = 8 * ntensors;
        traffic_tdotu = (9 + 1) * sizeof(Real) * ntensors;
        flop_div = 5 * 4 * ncells;
        traffic_div = (6 * 4 + 4) * sizeof(Real) * ncells;

        footprint = sizeof(*this);
    }

protected:
    // tensor slices
    TempPiXSOA_ST txx, txy, txz, utx; //slices for the tensors in the x-faces
    TempPiYSOA_ST tyx, tyy, tyz, uty; //slices for the tensors in the y-faces
    RingTempPiZSOA_ST ringtzx, ringtzy, ringtzz, ringutz; //slices for the tensors in the z-faces

    // output slices
    OutputSOA rhsu, rhsv, rhsw, rhse;

    virtual void _corners(const TInput& s0, const TInput& s1, TempSOA_ST& gradx, TempSOA_ST& grady, TempSOA_ST& gradz)
    {
        for(int iy=0; iy<TempSOA_ST::NY; iy++)
            for(int ix=0; ix<TempSOA_ST::NX; ix++)
            {
                const Real c000    = s1(ix,iy);
                const Real cm100   = s1(ix-1,iy);
                const Real c0m10   = s1(ix,iy-1);
                const Real cm1m10  = s1(ix-1,iy-1);
                const Real c00m1   = s0(ix,iy);
                const Real cm10m1  = s0(ix-1,iy);
                const Real c0m1m1  = s0(ix,iy-1);
                const Real cm1m1m1 = s0(ix-1,iy-1);

                gradx.ref(ix, iy) = c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1;
                grady.ref(ix, iy) = cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1  + c00m1-c0m1m1;
                gradz.ref(ix, iy) = cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1  + c000-c00m1;
            }
    }

    virtual void _udot_tx(const TInput& u, const TInput& v, const TInput& w)
    {
        for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
            for(int ix=0; ix<TempPiXSOA_ST::NX; ix++)
            {
                utx.ref(ix, iy) =
                txx(ix, iy)*(u(ix, iy) + u(ix-1, iy)) +
                txy(ix, iy)*(v(ix, iy) + v(ix-1, iy)) +
                txz(ix, iy)*(w(ix, iy) + w(ix-1, iy));
            }
    }

    virtual void _udot_ty(const TInput& u, const TInput& v, const TInput& w)
    {
        for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
            for(int ix=0; ix<TempPiYSOA_ST::NX; ix++)
            {
                uty.ref(ix, iy) =
                tyx(ix, iy)*(u(ix, iy) + u(ix, iy-1)) +
                tyy(ix, iy)*(v(ix, iy) + v(ix, iy-1)) +
                tyz(ix, iy)*(w(ix, iy) + w(ix, iy-1));
            }
    }

    virtual void _udot_tz(const TInput& u0, const TInput& v0, const TInput& w0,
                          const TInput& u1, const TInput& v1, const TInput& w1,
                          const TempPiZSOA_ST& tzx, const TempPiZSOA_ST& tzy, const TempPiZSOA_ST& tzz,
                          TempPiZSOA_ST& utz)
    {
        for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
            for(int ix=0; ix<TempPiZSOA_ST::NX; ix++)
            {
                utz.ref(ix, iy) =
                tzx(ix, iy)*(u0(ix, iy) + u1(ix, iy)) +
                tzy(ix, iy)*(v0(ix, iy) + v1(ix, iy)) +
                tzz(ix, iy)*(w0(ix, iy) + w1(ix, iy));
            }
    }

    virtual void _div_dxy()
    {
        for(int iy=0; iy<OutputSOA::NY; iy++)
            for(int ix=0; ix<OutputSOA::NX; ix++)
            {
                rhsu.ref(ix, iy) = txx(ix+1, iy) - txx(ix, iy) + tyx(ix, iy+1) - tyx(ix, iy);
                rhsv.ref(ix, iy) = txy(ix+1, iy) - txy(ix, iy) + tyy(ix, iy+1) - tyy(ix, iy);
                rhsw.ref(ix, iy) = txz(ix+1, iy) - txz(ix, iy) + tyz(ix, iy+1) - tyz(ix, iy);
                rhse.ref(ix, iy) = utx(ix+1, iy) - utx(ix, iy) + uty(ix, iy+1) - uty(ix, iy);
            }
    }

    virtual void _div_dz(const TempPiZSOA_ST& tzx0, const TempPiZSOA_ST& tzy0, const TempPiZSOA_ST& tzz0, const TempPiZSOA_ST& utz0,
                         const TempPiZSOA_ST& tzx1, const TempPiZSOA_ST& tzy1, const TempPiZSOA_ST& tzz1, const TempPiZSOA_ST& utz1)
    {
        for(int iy=0; iy<OutputSOA::NY; iy++)
            for(int ix=0; ix<OutputSOA::NX; ix++)
            {
                rhsu.ref(ix, iy) += tzx1(ix, iy) - tzx0(ix, iy);
                rhsv.ref(ix, iy) += tzy1(ix, iy) - tzy0(ix, iy);
                rhsw.ref(ix, iy) += tzz1(ix, iy) - tzz0(ix, iy);
                rhse.ref(ix, iy) += utz1(ix, iy) - utz0(ix, iy);
            }
    }

    inline void _tensors_next() { ringtzx.next(); ringtzy.next(); ringtzz.next(); ringutz.next(); }
};

#endif /* DIVTENSOR_CPP_H_DENX8LSZ */
