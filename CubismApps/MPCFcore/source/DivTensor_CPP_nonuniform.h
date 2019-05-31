/*
 *  DivTensor_CPP_nonuniform.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger on 7/24/17.
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef DIVTENSOR_CPP_NONUNIFORM_H_PTJBXLLM
#define DIVTENSOR_CPP_NONUNIFORM_H_PTJBXLLM

#include <cassert>
#include <cmath>

#include "DivTensor_CPP.h"
#include "FDCoefficients.h"


template <typename TInput>
class DivTensor_CPP_nonuniform : public virtual DivTensor_CPP<TInput>
{
public:
    DivTensor_CPP_nonuniform() : DivTensor_CPP<TInput>() {}

    virtual void hpc_info(float& flop_corners, int& traffic_corner,
                          float& flop_tdotu, int& traffic_tdotu,
                          float& flop_div, int& traffic_div,
                          size_t& footprint)
    {
        const int ncorners = (int)std::pow(_BLOCKSIZE_ + 1, 3);
        const int ntensors = 3 * (_BLOCKSIZE_ + 1) * _BLOCKSIZE_ * _BLOCKSIZE_;
        const int ncells = (int)std::pow(_BLOCKSIZE_, 3);

        // TODO: [fabianw@mavt.ethz.ch; Thu Jul 27 2017 09:53:08 AM (+0200)]
        // adjust this once implemented
        flop_corners = 21 * ncorners;
        traffic_corner = (8 + 3) * sizeof(Real) * ncorners;
        flop_tdotu = 8 * ntensors;
        traffic_tdotu = (9 + 1) * sizeof(Real) * ntensors;
        flop_div = 5 * 4 * ncells;
        traffic_div = (6 * 4 + 4) * sizeof(Real) * ncells;

        footprint = sizeof(*this);
    }

protected:
    ///////////////////////////////////////////////////////////////////////////
    // Gradient approximations:
    // Compute spatial derivatives at all cell corners using the eight cells
    // that share the corner.  Interpolates the derivatives computed at the
    // center of the cell faces to the corresponding corner.
    ///////////////////////////////////////////////////////////////////////////
    // 2nd-order
    virtual void _corners(const TInput& s_zm1, const TInput& s_z00, const TInput& s_zp1,
            const CoeffsFD_Block_t& xcoeffs,
            const CoeffsFD_Block_t& ycoeffs,
            const CoeffsFD_Block_t& zcoeffs,
            const int iz,
            TempSOA_ST& gradx, TempSOA_ST& grady, TempSOA_ST& gradz)
    {
        assert(CoeffsFD_Block_t::N_COEFFS == 3);
        const int* const idx = &StaticIndexArray::array::idx[0];
        const Real cz_m1 = zcoeffs.get_coeffs(-1)[iz];
        const Real cz_00 = zcoeffs.get_coeffs(0)[iz];
        const Real cz_p1 = zcoeffs.get_coeffs(1)[iz];

// TODO: [fabianw@mavt.ethz.ch; Thu Feb 22 2018 06:53:39 PM (+0100)]
// Remove redundant computations
        for(int iy=0; iy<TempSOA_ST::NY; iy++)
        {
            const int iy00 = idx[iy];
            const int iym1 = iy00-1;
            const int iyp1 = iy00+1;
            assert(iym1 >= -1 && iym1 < _BLOCKSIZE_-1);
            assert(iyp1 >= 1  && iyp1 < _BLOCKSIZE_+1);
            const Real cy_m1 = ycoeffs.get_coeffs(-1)[iy];
            const Real cy_00 = ycoeffs.get_coeffs(0)[iy];
            const Real cy_p1 = ycoeffs.get_coeffs(1)[iy];

            for(int ix=0; ix<TempSOA_ST::NX; ix++)
            {
                const int ix00 = idx[ix];
                const int ixm1 = ix00-1;
                const int ixp1 = ix00+1;
                assert(ixm1 >= -1 && ixm1 < _BLOCKSIZE_-1);
                assert(ixp1 >= 1  && ixp1 < _BLOCKSIZE_+1);
                const Real cx_m1 = xcoeffs.get_coeffs(-1)[ix];
                const Real cx_00 = xcoeffs.get_coeffs(0)[ix];
                const Real cx_p1 = xcoeffs.get_coeffs(1)[ix];

                // cell values
                const Real c00m1   = s_zm1(ix00, iy00);
                const Real cm10m1  = s_zm1(ixm1, iy00);
                const Real c0m1m1  = s_zm1(ix00, iym1);
                const Real cm1m1m1 = s_zm1(ixm1, iym1);
                const Real cp10m1  = s_zm1(ixp1, iy00);
                const Real cp1m1m1 = s_zm1(ixp1, iym1);
                const Real c0p1m1  = s_zm1(ix00, iyp1);
                const Real cm1p1m1 = s_zm1(ixm1, iyp1);

                const Real c000    = s_z00(ix00, iy00);
                const Real cm100   = s_z00(ixm1, iy00);
                const Real c0m10   = s_z00(ix00, iym1);
                const Real cm1m10  = s_z00(ixm1, iym1);
                const Real cp100   = s_z00(ixp1, iy00);
                const Real cp1m10  = s_z00(ixp1, iym1);
                const Real c0p10   = s_z00(ix00, iyp1);
                const Real cm1p10  = s_z00(ixm1, iyp1);

                const Real c00p1   = s_zp1(ix00, iy00);
                const Real cm10p1  = s_zp1(ixm1, iy00);
                const Real c0m1p1  = s_zp1(ix00, iym1);
                const Real cm1m1p1 = s_zp1(ixm1, iym1);

// TODO: [fabianw@mavt.ethz.ch; Thu Feb 22 2018 06:28:15 PM (+0100)]
// interpolation
                gradx.ref(ix, iy) = (
                          cx_m1*(cm100 + cm1m10 + cm10m1 + cm1m1m1)
                        + cx_00*(c000  + c0m10  + c00m1  + c0m1m1)
                        + cx_p1*(cp100 + cp1m10 + cp10m1 + cp1m1m1)
                        );

                grady.ref(ix, iy) = (
                          cy_m1*(c0m10 + cm1m10 + c0m1m1 + cm1m1m1)
                        + cy_00*(c000  + cm100  + c00m1  + cm10m1)
                        + cy_p1*(c0p10 + cm1p10 + c0p1m1 + cm1p1m1)
                        );

                gradz.ref(ix, iy) = (
                          cz_m1*(c00m1 + cm10m1 + c0m1m1 + cm1m1m1)
                        + cz_00*(c000  + cm100  + c0m10  + cm1m10)
                        + cz_p1*(c00p1 + cm10p1 + c0m1p1 + cm1m1p1)
                        );
            }
        }
    }

// TODO: [fabianw@mavt.ethz.ch; Fri Feb 23 2018 09:20:57 AM (+0100)]

    ///////////////////////////////////////////////////////////////////////////
    // 4th-order
    virtual void _corners(const TInput& sm2, const TInput& sm1, const TInput& s0, const TInput& sp1, const TInput& sp2,
            const CoeffsFD_Block_t& xcoeffs,
            const CoeffsFD_Block_t& ycoeffs,
            const CoeffsFD_Block_t& zcoeffs,
            const int iz,
            TempSOA_ST& gradx, TempSOA_ST& grady, TempSOA_ST& gradz)
    {
        // const Real cz = zcoeffs.a[iz];
        for(int iy=0; iy<TempSOA_ST::NY; iy++)
        {
            // const Real cy = ycoeffs.a[iy];
            for(int ix=0; ix<TempSOA_ST::NX; ix++)
            {
                // const Real cx = xcoeffs.a[ix];

                const Real c00p1   = sp1(ix,iy);
                const Real cm10p1  = sp1(ix-1,iy);
                const Real c0m1p1  = sp1(ix,iy-1);
                const Real cm1m1p1 = sp1(ix-1,iy-1);

                const Real cp100   = s0(ix+1,iy);
                const Real c000    = s0(ix,  iy);
                const Real cm100   = s0(ix-1,iy);
                const Real cm200   = s0(ix-2,iy);
                const Real cp1m10  = s0(ix+1,iy-1);
                const Real c0m10   = s0(ix,  iy-1);
                const Real cm1m10  = s0(ix-1,iy-1);
                const Real cm2m10  = s0(ix-2,iy-1);
                const Real c0p10   = s0(ix,  iy+1);
                const Real c0m20   = s0(ix,  iy-2);
                const Real cm1p10  = s0(ix-1,iy+1);
                const Real cm1m20  = s0(ix-1,iy-2);

                const Real cp10m1  = sm1(ix+1,iy);
                const Real c00m1   = sm1(ix,  iy);
                const Real cm10m1  = sm1(ix-1,iy);
                const Real cm20m1  = sm1(ix-2,iy);
                const Real cp1m1m1 = sm1(ix+1,iy-1);
                const Real c0m1m1  = sm1(ix,  iy-1);
                const Real cm1m1m1 = sm1(ix-1,iy-1);
                const Real cm2m1m1 = sm1(ix-2,iy-1);
                const Real c0p1m1  = sm1(ix,  iy+1);
                const Real c0m2m1  = sm1(ix,  iy-2);
                const Real cm1p1m1 = sm1(ix-1,iy+1);
                const Real cm1m2m1 = sm1(ix-1,iy-2);

                const Real c00m2   = sm2(ix,iy);
                const Real cm10m2  = sm2(ix-1,iy);
                const Real c0m1m2  = sm2(ix,iy-1);
                const Real cm1m1m2 = sm2(ix-1,iy-1);

                // gradx.ref(ix, iy) = cx*(-cp100+(Real)27.0*c000-(Real)27.0*cm100+cm200 -cp1m10+(Real)27.0*c0m10-(Real)27.0*cm1m10+cm2m10 -cp10m1+(Real)27.0*c00m1-(Real)27.0*cm10m1+cm20m1 -cp1m1m1+(Real)27.0*c0m1m1-(Real)27.0*cm1m1m1+cm2m1m1);
                // grady.ref(ix, iy) = cy*(-cm1p10+(Real)27.0*cm100-(Real)27.0*cm1m10+cm1m20 -c0p10+(Real)27.0*c000-(Real)27.0*c0m10+c0m20 -cm1p1m1+(Real)27.0*cm10m1-(Real)27.0*cm1m1m1+cm1m2m1 -c0p1m1+(Real)27.0*c00m1-(Real)27.0*c0m1m1+c0m2m1);
                // gradz.ref(ix, iy) = cz*(-cm1m1p1+(Real)27.0*cm1m10-(Real)27.0*cm1m1m1+cm1m1m2 -cm10p1+(Real)27.0*cm100-(Real)27.0*cm10m1+cm10m2 -c0m1p1+(Real)27.0*c0m10-(Real)27.0*c0m1m1+c0m1m2 -c00p1+(Real)27.0*c000-(Real)27.0*c00m1+c00m2);
                // gradx.ref(ix, iy) = cx*(-cp100+cm200 -cp1m10+cm2m10 -cp10m1+cm20m1 -cp1m1m1+cm2m1m1 + (Real)27.0*(c000-cm100 + c0m10-cm1m10 + c00m1-cm10m1 + c0m1m1-cm1m1m1));
                // grady.ref(ix, iy) = cy*(-cm1p10+cm1m20 -c0p10+c0m20 -cm1p1m1+cm1m2m1 -c0p1m1+c0m2m1 + (Real)27.0*(cm100-cm1m10 + c000-c0m10 + cm10m1-cm1m1m1 + c00m1-c0m1m1));
                // gradz.ref(ix, iy) = cz*(-cm1m1p1+cm1m1m2 -cm10p1+cm10m2 -c0m1p1+c0m1m2 -c00p1+c00m2 + (Real)27.0*(cm1m10-cm1m1m1 + cm100-cm10m1 + c0m10-c0m1m1 + c000-c00m1));
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Divergence approximation.  Currently 2nd-order (otherwise larger stencil
    // required)
    ///////////////////////////////////////////////////////////////////////////
    virtual void _div_dxy(const Real* const invh_x, const Real* const invh_y)
    {
        for(int iy=0; iy<OutputSOA::NY; iy++)
        {
            const Real ihy = invh_y[iy];
            for(int ix=0; ix<OutputSOA::NX; ix++)
            {
                const Real ihx = invh_x[ix];
                this->rhsu.ref(ix, iy) = ihx*(this->txx(ix+1, iy) - this->txx(ix, iy)) + ihy*(this->tyx(ix, iy+1) - this->tyx(ix, iy));
                this->rhsv.ref(ix, iy) = ihx*(this->txy(ix+1, iy) - this->txy(ix, iy)) + ihy*(this->tyy(ix, iy+1) - this->tyy(ix, iy));
                this->rhsw.ref(ix, iy) = ihx*(this->txz(ix+1, iy) - this->txz(ix, iy)) + ihy*(this->tyz(ix, iy+1) - this->tyz(ix, iy));
                this->rhse.ref(ix, iy) = ihx*(this->utx(ix+1, iy) - this->utx(ix, iy)) + ihy*(this->uty(ix, iy+1) - this->uty(ix, iy));
            }
        }
    }

    virtual void _div_dz(const TempPiZSOA_ST& tzx0, const TempPiZSOA_ST& tzy0, const TempPiZSOA_ST& tzz0, const TempPiZSOA_ST& utz0,
                         const TempPiZSOA_ST& tzx1, const TempPiZSOA_ST& tzy1, const TempPiZSOA_ST& tzz1, const TempPiZSOA_ST& utz1,
                         const Real ihz)
    {
        for(int iy=0; iy<OutputSOA::NY; iy++)
            for(int ix=0; ix<OutputSOA::NX; ix++)
            {
                this->rhsu.ref(ix, iy) += ihz*(tzx1(ix, iy) - tzx0(ix, iy));
                this->rhsv.ref(ix, iy) += ihz*(tzy1(ix, iy) - tzy0(ix, iy));
                this->rhsw.ref(ix, iy) += ihz*(tzz1(ix, iy) - tzz0(ix, iy));
                this->rhse.ref(ix, iy) += ihz*(utz1(ix, iy) - utz0(ix, iy));
            }
    }
};

#endif /* DIVTENSOR_CPP_NONUNIFORM_H_PTJBXLLM */
