/*
 *  VectorOperatorNouninform.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 08/21/2017
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef VECTOROPERATORNONUNIFORM_H_U4LE7NP8
#define VECTOROPERATORNONUNIFORM_H_U4LE7NP8

#include <string>
#include <Cubism/BlockInfo.h>
#include <Cubism/StencilInfo.h>
#include "common.h"
#include "Types.h"
using namespace cubism;

// Operator classes:
// Some of the streamers depend on a specific operator, which must be evaluated
// in advance and writes the result into a specific location.  When evaluating
// a set of operators, they must be evaluated in _increasing_ class order.  The
// following classes for operators with writes are defined:
// class 0:   no operator (default)
// class 1:   writes into dummy (scalar)              (MASK = 0x1;   INVALID = 0x0)
// class 2:   writes into tmp[iz][iy][ix][0] (scalar) (MASK = 0x2;   INVALID = 0x1)
// class 3:   writes into tmp[iz][iy][ix][1] (scalar) (MASK = 0x4;   INVALID = 0x3)
// class 4:   writes into tmp[iz][iy][ix][2] (scalar) (MASK = 0x8;   INVALID = 0x7)
// class 5:   writes into tmp[iz][iy][ix][3] (scalar) (MASK = 0x10;  INVALID = 0xf)
// class 6:   writes into tmp[iz][iy][ix][4] (scalar) (MASK = 0x20;  INVALID = 0x1f)
// class 7:   writes into tmp[iz][iy][ix][5] (scalar) (MASK = 0x40;  INVALID = 0x3f)
// class 8:   writes into tmp[iz][iy][ix][6] (scalar) (MASK = 0x80;  INVALID = 0x7f)
// class 9:   writes into tmp[iz][iy][ix][7] (scalar) (MASK = 0x100; INVALID = 0xff)
// class 10:  writes into tmp[iz][iy][ix][0-1] (two-component vector) (MASK = 0x200; INVALID = 0x1f9)
// class 11:  writes into tmp[iz][iy][ix][2-3] (two-component vector) (MASK = 0x400; INVALID = 0x3e7)
// class 12:  writes into tmp[iz][iy][ix][4-5] (two-component vector) (MASK = 0x800; INVALID = 0x79f)
// class 13:  writes into tmp[iz][iy][ix][6-7] (two-component vector) (MASK = 0x1000; INVALID = 0xe7f)
// class 14:  writes into tmp[iz][iy][ix][0-2] (three-component vector) (MASK = 0x2000; INVALID = 0x19f1)
// class 15:  writes into tmp[iz][iy][ix][3-5] (three-component vector) (MASK = 0x4000; INVALID = 0x338f)
// class 16:  writes into dummy + tmp[iz][iy][ix][0-7] (nine-component rank-2 tensor) (MASK = 0x8000; INVALID = 0x0)
#define OPMAXCLASS 17

// NOTE: Inclusion of this header must be guarded by #ifdef _NONUNIFORM_BLOCK_

// helper macros
#define __VELU(I,J,K) (lab((I),(J),(K)).ru/(lab((I),(J),(K)).alpha1rho1 + lab((I),(J),(K)).alpha2rho2))
#define __VELV(I,J,K) (lab((I),(J),(K)).rv/(lab((I),(J),(K)).alpha1rho1 + lab((I),(J),(K)).alpha2rho2))
#define __VELW(I,J,K) (lab((I),(J),(K)).rw/(lab((I),(J),(K)).alpha1rho1 + lab((I),(J),(K)).alpha2rho2))
#define __ALPHA(I,J,K) (lab((I),(J),(K)).alpha2)

#define __FD_2ND(i,c,um1,u00,up1) (c.get_coeffs(-1)[i]*um1 + c.get_coeffs(0)[i]*u00 + c.get_coeffs(1)[i]*up1)
#define __FD_4TH(i,c,um2,um1,u00,up1,up2) (c.get_coeffs(-2)[i]*um2 + c.get_coeffs(-1)[i]*um1 + c.get_coeffs(0)[i]*u00 + c.get_coeffs(1)[i]*up1 + c.get_coeffs(2)[i]*up2)

///////////////////////////////////////////////////////////////////////////////
// Finite differences based on FDCoefficients.h:
// Derivatives are evaluated at cell faces.  The cell-center value is obtained
// by averaging the two values on the faces.
///////////////////////////////////////////////////////////////////////////////
struct Operator_Vort
{
    static const std::string NAME;
    static const int CLASS = 14;
    static const unsigned int MASK  = 0x2000;
    static const unsigned int INVALID = 0x19f1;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_Vort(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
    Operator_Vort(const Operator_Vort& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
#else
    Operator_Vort(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_Vort(const Operator_Vort& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        // scratch pads used as temporary buffers
        Real sx[2];
        Real sy[TBlock::sizeX][2];
        Real sz[TBlock::sizeY][TBlock::sizeX][2];

        for(int iy=0; iy<TBlock::sizeY; iy++)
        {
            for(int ix=0; ix<TBlock::sizeX; ix++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sz[iy][ix][0] = __FD_2ND(0,cz,__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1));
                sz[iy][ix][1] = __FD_2ND(0,cz,__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1));
#else
                sz[iy][ix][0] = __FD_4TH(0,cz,__VELU(ix,iy,-2),__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1),__VELU(ix,iy,2));
                sz[iy][ix][1] = __FD_4TH(0,cz,__VELV(ix,iy,-2),__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1),__VELV(ix,iy,2));
#endif /* _SECONDORDER_FD_CORE_ */
            }
        }

        for(int iz=0; iz<TBlock::sizeZ; iz++)
        {
            for (int ix = 0; ix < TBlock::sizeX; ++ix)
            {
#ifdef _SECONDORDER_FD_CORE_
                sy[ix][0] = __FD_2ND(0,cy,__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz));
                sy[ix][1] = __FD_2ND(0,cy,__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz));
#else
                sy[ix][0] = __FD_4TH(0,cy,__VELU(ix,-2,iz),__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz),__VELU(ix,2,iz));
                sy[ix][1] = __FD_4TH(0,cy,__VELW(ix,-2,iz),__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz),__VELW(ix,2,iz));
#endif /* _SECONDORDER_FD_CORE_ */
            }

            for(int iy=0; iy<TBlock::sizeY; iy++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sx[0] = __FD_2ND(0,cx,__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz));
                sx[1] = __FD_2ND(0,cx,__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz));
#else
                sx[0] = __FD_4TH(0,cx,__VELV(-2,iy,iz),__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz),__VELV(2,iy,iz));
                sx[1] = __FD_4TH(0,cx,__VELW(-2,iy,iz),__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz),__VELW(2,iy,iz));
#endif /* _SECONDORDER_FD_CORE_ */
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
#ifdef _SECONDORDER_FD_CORE_
                    const Real vx1 = __FD_2ND(ix+1,cx,__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz));
                    const Real uy1 = __FD_2ND(iy+1,cy,__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz));
                    const Real uz1 = __FD_2ND(iz+1,cz,__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1));
                    const Real wx1 = __FD_2ND(ix+1,cx,__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz));
                    const Real wy1 = __FD_2ND(iy+1,cy,__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz));
                    const Real vz1 = __FD_2ND(iz+1,cz,__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1));
#else
                    const Real vx1 = __FD_4TH(ix+1,cx,__VELV(ix-2,iy,iz),__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz),__VELV(ix+2,iy,iz));
                    const Real uy1 = __FD_4TH(iy+1,cy,__VELU(ix,iy-2,iz),__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz),__VELU(ix,iy+2,iz));
                    const Real uz1 = __FD_4TH(iz+1,cz,__VELU(ix,iy,iz-2),__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1),__VELU(ix,iy,iz+2));
                    const Real wx1 = __FD_4TH(ix+1,cx,__VELW(ix-2,iy,iz),__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz),__VELW(ix+2,iy,iz));
                    const Real wy1 = __FD_4TH(iy+1,cy,__VELW(ix,iy-2,iz),__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz),__VELW(ix,iy+2,iz));
                    const Real vz1 = __FD_4TH(iz+1,cz,__VELV(ix,iy,iz-2),__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1),__VELV(ix,iy,iz+2));
#endif /* _SECONDORDER_FD_CORE_ */

                    o.tmp[iz][iy][ix][0] = static_cast<Real>(0.5)*( (sy[ix][1]+wy1) - (sz[iy][ix][1]+vz1) );
                    o.tmp[iz][iy][ix][1] = static_cast<Real>(0.5)*( (sz[iy][ix][0]+uz1) - (sx[1]+wx1) );
                    o.tmp[iz][iy][ix][2] = static_cast<Real>(0.5)*( (sx[0]+vx1) - (sy[ix][0]+uy1) );

                    sx[0]         = vx1;
                    sy[ix][0]     = uy1;
                    sz[iy][ix][0] = uz1;
                    sx[1]         = wx1;
                    sy[ix][1]     = wy1;
                    sz[iy][ix][1] = vz1;
                }
            }
        }
    }
};


struct Operator_IgradA2I
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_IgradA2I(): stencil(-1,-1,-1,2,2,2, false, 1,6) {}
    Operator_IgradA2I(const Operator_IgradA2I& c): stencil(-1,-1,-1,2,2,2, false, 1, 6) {}
#else
    Operator_IgradA2I(): stencil(-2,-2,-2,3,3,3, false, 1,6) {}
    Operator_IgradA2I(const Operator_IgradA2I& c): stencil(-2,-2,-2,3,3,3, false, 1, 6) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        // scratch pads
        Real sx;
        Real sy[TBlock::sizeX];
        Real sz[TBlock::sizeY][TBlock::sizeX];

        for(int iy=0; iy<TBlock::sizeY; iy++)
            for(int ix=0; ix<TBlock::sizeX; ix++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sz[iy][ix] = __FD_2ND(0,cz,__ALPHA(ix,iy,-1),__ALPHA(ix,iy,0),__ALPHA(ix,iy,1));
#else
                sz[iy][ix] = __FD_4TH(0,cz,__ALPHA(ix,iy,-2),__ALPHA(ix,iy,-1),__ALPHA(ix,iy,0),__ALPHA(ix,iy,1),__ALPHA(ix,iy,2));
#endif /* _SECONDORDER_FD_CORE_ */
            }

        for(int iz=0; iz<TBlock::sizeZ; iz++)
        {
            for (int ix = 0; ix < TBlock::sizeX; ++ix)
            {
#ifdef _SECONDORDER_FD_CORE_
                sy[ix] = __FD_2ND(0,cy,__ALPHA(ix,-1,iz),__ALPHA(ix,0,iz),__ALPHA(ix,1,iz));
#else
                sy[ix] = __FD_4TH(0,cy,__ALPHA(ix,-2,iz),__ALPHA(ix,-1,iz),__ALPHA(ix,0,iz),__ALPHA(ix,1,iz),__ALPHA(ix,2,iz));
#endif /* _SECONDORDER_FD_CORE_ */
            }

            for(int iy=0; iy<TBlock::sizeY; iy++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sx = __FD_2ND(0,cx,__ALPHA(-1,iy,iz),__ALPHA(0,iy,iz),__ALPHA(1,iy,iz));
#else
                sx = __FD_4TH(0,cx,__ALPHA(-2,iy,iz),__ALPHA(-1,iy,iz),__ALPHA(0,iy,iz),__ALPHA(1,iy,iz),__ALPHA(2,iy,iz));
#endif /* _SECONDORDER_FD_CORE_ */
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
#ifdef _SECONDORDER_FD_CORE_
                    const Real ax1 = __FD_2ND(ix+1,cx,__ALPHA(ix-1,iy,iz),__ALPHA(ix,iy,iz),__ALPHA(ix+1,iy,iz));
                    const Real ay1 = __FD_2ND(iy+1,cy,__ALPHA(ix,iy-1,iz),__ALPHA(ix,iy,iz),__ALPHA(ix,iy+1,iz));
                    const Real az1 = __FD_2ND(iz+1,cz,__ALPHA(ix,iy,iz-1),__ALPHA(ix,iy,iz),__ALPHA(ix,iy,iz+1));
#else
                    const Real ax1 = __FD_4TH(ix+1,cx,__ALPHA(ix-2,iy,iz),__ALPHA(ix-1,iy,iz),__ALPHA(ix,iy,iz),__ALPHA(ix+1,iy,iz),__ALPHA(ix+2,iy,iz));
                    const Real ay1 = __FD_4TH(iy+1,cy,__ALPHA(ix,iy-2,iz),__ALPHA(ix,iy-1,iz),__ALPHA(ix,iy,iz),__ALPHA(ix,iy+1,iz),__ALPHA(ix,iy+2,iz));
                    const Real az1 = __FD_4TH(iz+1,cz,__ALPHA(ix,iy,iz-2),__ALPHA(ix,iy,iz-1),__ALPHA(ix,iy,iz),__ALPHA(ix,iy,iz+1),__ALPHA(ix,iy,iz+2));
#endif /* _SECONDORDER_FD_CORE_ */

                    const Real ax = static_cast<Real>(0.5)*(sx+ax1);
                    const Real ay = static_cast<Real>(0.5)*(sy[ix]+ay1);
                    const Real az = static_cast<Real>(0.5)*(sz[iy][ix]+az1);
                    o(ix,iy,iz).dummy = mysqrt(ax*ax + ay*ay + az*az);

                    sx = ax1;
                    sy[ix] = ay1;
                    sz[iy][ix] = az1;
                }
            }
        }
    }
};

struct Operator_divU
{
    static const std::string NAME;
    static const int CLASS = 5;
    static const unsigned int MASK  = 0x10;
    static const unsigned int INVALID = 0xf;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_divU(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
    Operator_divU(const Operator_divU& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
#else
    Operator_divU(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_divU(const Operator_divU& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        // scratch pads
        Real sx;
        Real sy[TBlock::sizeX];
        Real sz[TBlock::sizeY][TBlock::sizeX];

        for(int iy=0; iy<TBlock::sizeY; iy++)
            for(int ix=0; ix<TBlock::sizeX; ix++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sz[iy][ix] = __FD_2ND(0,cz,__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1));
#else
                sz[iy][ix] = __FD_4TH(0,cz,__VELW(ix,iy,-2),__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1),__VELW(ix,iy,2));
#endif /* _SECONDORDER_FD_CORE_ */
            }

        for(int iz=0; iz<TBlock::sizeZ; iz++)
        {
            for (int ix = 0; ix < TBlock::sizeX; ++ix)
            {
#ifdef _SECONDORDER_FD_CORE_
                sy[ix] = __FD_2ND(0,cy,__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz));
#else
                sy[ix] = __FD_4TH(0,cy,__VELV(ix,-2,iz),__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz),__VELV(ix,2,iz));
#endif /* _SECONDORDER_FD_CORE_ */
            }

            for(int iy=0; iy<TBlock::sizeY; iy++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sx = __FD_2ND(0,cx,__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz));
#else
                sx = __FD_4TH(0,cx,__VELU(-2,iy,iz),__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz),__VELU(2,iy,iz));
#endif /* _SECONDORDER_FD_CORE_ */
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
#ifdef _SECONDORDER_FD_CORE_
                    const Real ux1 = __FD_2ND(ix+1,cx,__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz));
                    const Real vy1 = __FD_2ND(iy+1,cy,__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz));
                    const Real wz1 = __FD_2ND(iz+1,cz,__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1));
#else
                    const Real ux1 = __FD_4TH(ix+1,cx,__VELU(ix-2,iy,iz),__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz),__VELU(ix+2,iy,iz));
                    const Real vy1 = __FD_4TH(iy+1,cy,__VELV(ix,iy-2,iz),__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz),__VELV(ix,iy+2,iz));
                    const Real wz1 = __FD_4TH(iz+1,cz,__VELW(ix,iy,iz-2),__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1),__VELW(ix,iy,iz+2));
#endif /* _SECONDORDER_FD_CORE_ */

                    o.tmp[iz][iy][ix][3] = static_cast<Real>(0.5)*((sx+ux1) + (sy[ix]+vy1) + (sz[iy][ix]+wz1)); // 0.5 for averaging

                    sx = ux1;
                    sy[ix] = vy1;
                    sz[iy][ix] = wz1;
                }
            }
        }
    }
};


struct Operator_SSDeviatoric
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_SSDeviatoric(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
    Operator_SSDeviatoric(const Operator_SSDeviatoric& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
#else
    Operator_SSDeviatoric(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_SSDeviatoric(const Operator_SSDeviatoric& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        // scratch pads
        Real sx[3];
        Real sy[TBlock::sizeX][3];
        Real sz[TBlock::sizeY][TBlock::sizeX][3];

        for(int iy=0; iy<TBlock::sizeY; iy++)
            for(int ix=0; ix<TBlock::sizeX; ix++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sz[iy][ix][0] = __FD_2ND(0,cz,__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1));
                sz[iy][ix][1] = __FD_2ND(0,cz,__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1));
                sz[iy][ix][2] = __FD_2ND(0,cz,__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1));
#else
                sz[iy][ix][0] = __FD_4TH(0,cz,__VELU(ix,iy,-2),__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1),__VELU(ix,iy,2));
                sz[iy][ix][1] = __FD_4TH(0,cz,__VELV(ix,iy,-2),__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1),__VELV(ix,iy,2));
                sz[iy][ix][2] = __FD_4TH(0,cz,__VELW(ix,iy,-2),__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1),__VELW(ix,iy,2));
#endif /* _SECONDORDER_FD_CORE_ */
            }

        for(int iz=0; iz<TBlock::sizeZ; iz++)
        {
            for (int ix = 0; ix < TBlock::sizeX; ++ix)
            {
#ifdef _SECONDORDER_FD_CORE_
                sy[ix][0] = __FD_2ND(0,cy,__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz));
                sy[ix][1] = __FD_2ND(0,cy,__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz));
                sy[ix][2] = __FD_2ND(0,cy,__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz));
#else
                sy[ix][0] = __FD_4TH(0,cy,__VELU(ix,-2,iz),__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz),__VELU(ix,2,iz));
                sy[ix][1] = __FD_4TH(0,cy,__VELV(ix,-2,iz),__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz),__VELV(ix,2,iz));
                sy[ix][2] = __FD_4TH(0,cy,__VELW(ix,-2,iz),__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz),__VELW(ix,2,iz));
#endif /* _SECONDORDER_FD_CORE_ */
            }

            for(int iy=0; iy<TBlock::sizeY; iy++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sx[0] = __FD_2ND(0,cx,__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz));
                sx[1] = __FD_2ND(0,cx,__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz));
                sx[2] = __FD_2ND(0,cx,__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz));
#else
                sx[0] = __FD_4TH(0,cx,__VELU(-2,iy,iz),__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz),__VELU(2,iy,iz));
                sx[1] = __FD_4TH(0,cx,__VELV(-2,iy,iz),__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz),__VELV(2,iy,iz));
                sx[2] = __FD_4TH(0,cx,__VELW(-2,iy,iz),__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz),__VELW(2,iy,iz));
#endif /* _SECONDORDER_FD_CORE_ */
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
#ifdef _SECONDORDER_FD_CORE_
                    const Real ux1 = __FD_2ND(ix+1,cx,__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz));
                    const Real uy1 = __FD_2ND(iy+1,cy,__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz));
                    const Real uz1 = __FD_2ND(iz+1,cz,__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1));
                    const Real vx1 = __FD_2ND(ix+1,cx,__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz));
                    const Real vy1 = __FD_2ND(iy+1,cy,__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz));
                    const Real vz1 = __FD_2ND(iz+1,cz,__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1));
                    const Real wx1 = __FD_2ND(ix+1,cx,__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz));
                    const Real wy1 = __FD_2ND(iy+1,cy,__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz));
                    const Real wz1 = __FD_2ND(iz+1,cz,__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1));
#else
                    const Real ux1 = __FD_4TH(ix+1,cx,__VELU(ix-2,iy,iz),__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz),__VELU(ix+2,iy,iz));
                    const Real uy1 = __FD_4TH(iy+1,cy,__VELU(ix,iy-2,iz),__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz),__VELU(ix,iy+2,iz));
                    const Real uz1 = __FD_4TH(iz+1,cz,__VELU(ix,iy,iz-2),__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1),__VELU(ix,iy,iz+2));
                    const Real vx1 = __FD_4TH(ix+1,cx,__VELV(ix-2,iy,iz),__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz),__VELV(ix+2,iy,iz));
                    const Real vy1 = __FD_4TH(iy+1,cy,__VELV(ix,iy-2,iz),__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz),__VELV(ix,iy+2,iz));
                    const Real vz1 = __FD_4TH(iz+1,cz,__VELV(ix,iy,iz-2),__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1),__VELV(ix,iy,iz+2));
                    const Real wx1 = __FD_4TH(ix+1,cx,__VELW(ix-2,iy,iz),__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz),__VELW(ix+2,iy,iz));
                    const Real wy1 = __FD_4TH(iy+1,cy,__VELW(ix,iy-2,iz),__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz),__VELW(ix,iy+2,iz));
                    const Real wz1 = __FD_4TH(iz+1,cz,__VELW(ix,iy,iz-2),__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1),__VELW(ix,iy,iz+2));
#endif /* _SECONDORDER_FD_CORE_ */

                    const Real dudx = (sx[0]+ux1);
                    const Real dudy = (sy[ix][0]+uy1);
                    const Real dudz = (sz[ix][iy][0]+uz1);
                    const Real dvdx = (sx[1]+vx1);
                    const Real dvdy = (sy[ix][1]+vy1);
                    const Real dvdz = (sz[iy][ix][1]+vz1);
                    const Real dwdx = (sx[2]+wx1);
                    const Real dwdy = (sy[ix][2]+wy1);
                    const Real dwdz = (sz[iy][ix][2]+wz1);

                    const Real t1 = static_cast<Real>(2.0)*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
                    + (dvdx + dudy)*dudy + (dwdx + dudz)*dudz
                    + (dudy + dvdx)*dvdx + (dwdy + dvdz)*dvdz
                    + (dudz + dwdx)*dwdx + (dvdz + dwdy)*dwdy;

                    const Real t2 = (dudx + dvdy + dwdz)*(dudx + dvdy + dwdz);

                    o(ix,iy,iz).dummy = static_cast<Real>(0.25)*(static_cast<Real>(0.5)*t1 - static_cast<Real>(1.0/3.0)*t2); // 0.25 from averaging

                    sx[0] = ux1;
                    sy[ix][0] = uy1;
                    sz[iy][ix][0] = uz1;
                    sx[1] = vx1;
                    sy[ix][1] = vy1;
                    sz[iy][ix][1] = vz1;
                    sx[2] = wx1;
                    sy[ix][2] = wy1;
                    sz[iy][ix][2] = wz1;
                }
            }
        }
    }
};


struct Operator_Qcrit
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_Qcrit(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
    Operator_Qcrit(const Operator_Qcrit& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
#else
    Operator_Qcrit(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_Qcrit(const Operator_Qcrit& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        // scratch pads
        Real sx[3];
        Real sy[TBlock::sizeX][3];
        Real sz[TBlock::sizeY][TBlock::sizeX][3];

        for(int iy=0; iy<TBlock::sizeY; iy++)
            for(int ix=0; ix<TBlock::sizeX; ix++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sz[iy][ix][0] = __FD_2ND(0,cz,__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1));
                sz[iy][ix][1] = __FD_2ND(0,cz,__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1));
                sz[iy][ix][2] = __FD_2ND(0,cz,__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1));
#else
                sz[iy][ix][0] = __FD_4TH(0,cz,__VELU(ix,iy,-2),__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1),__VELU(ix,iy,2));
                sz[iy][ix][1] = __FD_4TH(0,cz,__VELV(ix,iy,-2),__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1),__VELV(ix,iy,2));
                sz[iy][ix][2] = __FD_4TH(0,cz,__VELW(ix,iy,-2),__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1),__VELW(ix,iy,2));
#endif /* _SECONDORDER_FD_CORE_ */
            }

        for(int iz=0; iz<TBlock::sizeZ; iz++)
        {
            for (int ix = 0; ix < TBlock::sizeX; ++ix)
            {
#ifdef _SECONDORDER_FD_CORE_
                sy[ix][0] = __FD_2ND(0,cy,__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz));
                sy[ix][1] = __FD_2ND(0,cy,__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz));
                sy[ix][2] = __FD_2ND(0,cy,__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz));
#else
                sy[ix][0] = __FD_4TH(0,cy,__VELU(ix,-2,iz),__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz),__VELU(ix,2,iz));
                sy[ix][1] = __FD_4TH(0,cy,__VELV(ix,-2,iz),__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz),__VELV(ix,2,iz));
                sy[ix][2] = __FD_4TH(0,cy,__VELW(ix,-2,iz),__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz),__VELW(ix,2,iz));
#endif /* _SECONDORDER_FD_CORE_ */
            }

            for(int iy=0; iy<TBlock::sizeY; iy++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sx[0] = __FD_2ND(0,cx,__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz));
                sx[1] = __FD_2ND(0,cx,__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz));
                sx[2] = __FD_2ND(0,cx,__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz));
#else
                sx[0] = __FD_4TH(0,cx,__VELU(-2,iy,iz),__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz),__VELU(2,iy,iz));
                sx[1] = __FD_4TH(0,cx,__VELV(-2,iy,iz),__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz),__VELV(2,iy,iz));
                sx[2] = __FD_4TH(0,cx,__VELW(-2,iy,iz),__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz),__VELW(2,iy,iz));
#endif /* _SECONDORDER_FD_CORE_ */
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
#ifdef _SECONDORDER_FD_CORE_
                    const Real ux1 = __FD_2ND(ix+1,cx,__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz));
                    const Real uy1 = __FD_2ND(iy+1,cy,__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz));
                    const Real uz1 = __FD_2ND(iz+1,cz,__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1));
                    const Real vx1 = __FD_2ND(ix+1,cx,__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz));
                    const Real vy1 = __FD_2ND(iy+1,cy,__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz));
                    const Real vz1 = __FD_2ND(iz+1,cz,__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1));
                    const Real wx1 = __FD_2ND(ix+1,cx,__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz));
                    const Real wy1 = __FD_2ND(iy+1,cy,__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz));
                    const Real wz1 = __FD_2ND(iz+1,cz,__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1));
#else
                    const Real ux1 = __FD_4TH(ix+1,cx,__VELU(ix-2,iy,iz),__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz),__VELU(ix+2,iy,iz));
                    const Real uy1 = __FD_4TH(iy+1,cy,__VELU(ix,iy-2,iz),__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz),__VELU(ix,iy+2,iz));
                    const Real uz1 = __FD_4TH(iz+1,cz,__VELU(ix,iy,iz-2),__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1),__VELU(ix,iy,iz+2));
                    const Real vx1 = __FD_4TH(ix+1,cx,__VELV(ix-2,iy,iz),__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz),__VELV(ix+2,iy,iz));
                    const Real vy1 = __FD_4TH(iy+1,cy,__VELV(ix,iy-2,iz),__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz),__VELV(ix,iy+2,iz));
                    const Real vz1 = __FD_4TH(iz+1,cz,__VELV(ix,iy,iz-2),__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1),__VELV(ix,iy,iz+2));
                    const Real wx1 = __FD_4TH(ix+1,cx,__VELW(ix-2,iy,iz),__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz),__VELW(ix+2,iy,iz));
                    const Real wy1 = __FD_4TH(iy+1,cy,__VELW(ix,iy-2,iz),__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz),__VELW(ix,iy+2,iz));
                    const Real wz1 = __FD_4TH(iz+1,cz,__VELW(ix,iy,iz-2),__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1),__VELW(ix,iy,iz+2));
#endif /* _SECONDORDER_FD_CORE_ */

                    const Real ux = (sx[0]+ux1);
                    const Real uy = (sy[ix][0]+uy1);
                    const Real uz = (sz[ix][iy][0]+uz1);
                    const Real vx = (sx[1]+vx1);
                    const Real vy = (sy[ix][1]+vy1);
                    const Real vz = (sz[iy][ix][1]+vz1);
                    const Real wx = (sx[2]+wx1);
                    const Real wy = (sy[ix][2]+wy1);
                    const Real wz = (sz[iy][ix][2]+wz1);

                    // rate of rotation / rate of strain
                    const Real O[9] = {
                        0.0,         0.5*(uy-vx), 0.5*(uz-wx),
                        0.5*(vx-uy), 0.0,         0.5*(vz-wy),
                        0.5*(wx-uz), 0.5*(wy-vz), 0.0};

                    const Real S[9] = {
                        ux,          0.5*(uy+vx), 0.5*(uz+wx),
                        0.5*(vx+uy), vy,          0.5*(vz+wy),
                        0.5*(wx+uz), 0.5*(wy+vz), wz};

                    // frobenius norm
                    Real IOI2 = 0.0;
                    Real ISI2 = 0.0;
                    for (int i=0; i<9; ++i)
                    {
                        IOI2 += O[i]*O[i];
                        ISI2 += S[i]*S[i];
                    }

                    o(ix,iy,iz).dummy = static_cast<Real>(0.25)*(static_cast<Real>(0.5)*(IOI2 - ISI2)); // 0.25 from averaging

                    sx[0] = ux1;
                    sy[ix][0] = uy1;
                    sz[iy][ix][0] = uz1;
                    sx[1] = vx1;
                    sy[ix][1] = vy1;
                    sz[iy][ix][1] = vz1;
                    sx[2] = wx1;
                    sy[ix][2] = wy1;
                    sz[iy][ix][2] = wz1;
                }
            }
        }
    }
};


struct Operator_gradU
{
    static const std::string NAME;
    static const int CLASS = 16;
    static const unsigned int MASK  = 0x8000;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_gradU(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
    Operator_gradU(const Operator_gradU& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
#else
    Operator_gradU(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_gradU(const Operator_gradU& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        // scratch pads
        Real sx[3];
        Real sy[TBlock::sizeX][3];
        Real sz[TBlock::sizeY][TBlock::sizeX][3];

        for(int iy=0; iy<TBlock::sizeY; iy++)
            for(int ix=0; ix<TBlock::sizeX; ix++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sz[iy][ix][0] = __FD_2ND(0,cz,__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1));
                sz[iy][ix][1] = __FD_2ND(0,cz,__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1));
                sz[iy][ix][2] = __FD_2ND(0,cz,__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1));
#else
                sz[iy][ix][0] = __FD_4TH(0,cz,__VELU(ix,iy,-2),__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1),__VELU(ix,iy,2));
                sz[iy][ix][1] = __FD_4TH(0,cz,__VELV(ix,iy,-2),__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1),__VELV(ix,iy,2));
                sz[iy][ix][2] = __FD_4TH(0,cz,__VELW(ix,iy,-2),__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1),__VELW(ix,iy,2));
#endif /* _SECONDORDER_FD_CORE_ */
            }

        for(int iz=0; iz<TBlock::sizeZ; iz++)
        {
            for (int ix = 0; ix < TBlock::sizeX; ++ix)
            {
#ifdef _SECONDORDER_FD_CORE_
                sy[ix][0] = __FD_2ND(0,cy,__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz));
                sy[ix][1] = __FD_2ND(0,cy,__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz));
                sy[ix][2] = __FD_2ND(0,cy,__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz));
#else
                sy[ix][0] = __FD_4TH(0,cy,__VELU(ix,-2,iz),__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz),__VELU(ix,2,iz));
                sy[ix][1] = __FD_4TH(0,cy,__VELV(ix,-2,iz),__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz),__VELV(ix,2,iz));
                sy[ix][2] = __FD_4TH(0,cy,__VELW(ix,-2,iz),__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz),__VELW(ix,2,iz));
#endif /* _SECONDORDER_FD_CORE_ */
            }

            for(int iy=0; iy<TBlock::sizeY; iy++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sx[0] = __FD_2ND(0,cx,__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz));
                sx[1] = __FD_2ND(0,cx,__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz));
                sx[2] = __FD_2ND(0,cx,__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz));
#else
                sx[0] = __FD_4TH(0,cx,__VELU(-2,iy,iz),__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz),__VELU(2,iy,iz));
                sx[1] = __FD_4TH(0,cx,__VELV(-2,iy,iz),__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz),__VELV(2,iy,iz));
                sx[2] = __FD_4TH(0,cx,__VELW(-2,iy,iz),__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz),__VELW(2,iy,iz));
#endif /* _SECONDORDER_FD_CORE_ */

                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
#ifdef _SECONDORDER_FD_CORE_
                    const Real ux1 = __FD_2ND(ix+1,cx,__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz));
                    const Real uy1 = __FD_2ND(iy+1,cy,__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz));
                    const Real uz1 = __FD_2ND(iz+1,cz,__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1));
                    const Real vx1 = __FD_2ND(ix+1,cx,__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz));
                    const Real vy1 = __FD_2ND(iy+1,cy,__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz));
                    const Real vz1 = __FD_2ND(iz+1,cz,__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1));
                    const Real wx1 = __FD_2ND(ix+1,cx,__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz));
                    const Real wy1 = __FD_2ND(iy+1,cy,__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz));
                    const Real wz1 = __FD_2ND(iz+1,cz,__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1));
#else
                    const Real ux1 = __FD_4TH(ix+1,cx,__VELU(ix-2,iy,iz),__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz),__VELU(ix+2,iy,iz));
                    const Real uy1 = __FD_4TH(iy+1,cy,__VELU(ix,iy-2,iz),__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz),__VELU(ix,iy+2,iz));
                    const Real uz1 = __FD_4TH(iz+1,cz,__VELU(ix,iy,iz-2),__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1),__VELU(ix,iy,iz+2));
                    const Real vx1 = __FD_4TH(ix+1,cx,__VELV(ix-2,iy,iz),__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz),__VELV(ix+2,iy,iz));
                    const Real vy1 = __FD_4TH(iy+1,cy,__VELV(ix,iy-2,iz),__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz),__VELV(ix,iy+2,iz));
                    const Real vz1 = __FD_4TH(iz+1,cz,__VELV(ix,iy,iz-2),__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1),__VELV(ix,iy,iz+2));
                    const Real wx1 = __FD_4TH(ix+1,cx,__VELW(ix-2,iy,iz),__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz),__VELW(ix+2,iy,iz));
                    const Real wy1 = __FD_4TH(iy+1,cy,__VELW(ix,iy-2,iz),__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz),__VELW(ix,iy+2,iz));
                    const Real wz1 = __FD_4TH(iz+1,cz,__VELW(ix,iy,iz-2),__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1),__VELW(ix,iy,iz+2));
#endif /* _SECONDORDER_FD_CORE_ */

                    o(ix,iy,iz).dummy =    static_cast<Real>(0.5)*(sx[0]+ux1);
                    o.tmp[iz][iy][ix][0] = static_cast<Real>(0.5)*(sy[ix][0]+uy1);
                    o.tmp[iz][iy][ix][1] = static_cast<Real>(0.5)*(sz[iy][ix][0]+uz1);
                    o.tmp[iz][iy][ix][2] = static_cast<Real>(0.5)*(sx[1]+vx1);
                    o.tmp[iz][iy][ix][3] = static_cast<Real>(0.5)*(sy[ix][1]+vy1);
                    o.tmp[iz][iy][ix][4] = static_cast<Real>(0.5)*(sz[iy][ix][1]+vz1);
                    o.tmp[iz][iy][ix][5] = static_cast<Real>(0.5)*(sx[2]+wx1);
                    o.tmp[iz][iy][ix][6] = static_cast<Real>(0.5)*(sy[ix][2]+wy1);
                    o.tmp[iz][iy][ix][7] = static_cast<Real>(0.5)*(sz[iy][ix][2]+wz1);

                    sx[0]         = ux1;
                    sy[ix][0]     = uy1;
                    sz[iy][ix][0] = uz1;
                    sx[1]         = vx1;
                    sy[ix][1]     = vy1;
                    sz[iy][ix][1] = vz1;
                    sx[2]         = wx1;
                    sy[ix][2]     = wy1;
                    sz[iy][ix][2] = wz1;
                }
            }
        }
    }
};


///////////////////////////////////////////////////////////////////////////////
// special
///////////////////////////////////////////////////////////////////////////////
template <int comp>
struct Operator_gradUij
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_gradUij(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
    Operator_gradUij(const Operator_gradUij& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
#else
    Operator_gradUij(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_gradUij(const Operator_gradUij& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const int* const idx = &StaticIndexArray::array::idx[0];
        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
// #ifdef _SECONDORDER_FD_CORE_
//                     if (0 == comp)
//                     {
//                         const Real v0 = cx[ix]*(__VELU(ix,iy,iz) - __VELU(ix-1,iy,iz));
//                         const Real v1 = cx[ix+1]*(__VELU(ix+1,iy,iz) - __VELU(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (1 == comp)
//                     {
//                         const Real v0 = cy[iy]*(__VELU(ix,iy,iz) - __VELU(ix,iy-1,iz));
//                         const Real v1 = cy[iy+1]*(__VELU(ix,iy+1,iz) - __VELU(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (2 == comp)
//                     {
//                         const Real v0 = cz[iz]*(__VELU(ix,iy,iz) - __VELU(ix,iy,iz-1));
//                         const Real v1 = cz[iz+1]*(__VELU(ix,iy,iz+1) - __VELU(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (3 == comp)
//                     {
//                         const Real v0 = cx[ix]*(__VELV(ix,iy,iz) - __VELV(ix-1,iy,iz));
//                         const Real v1 = cx[ix+1]*(__VELV(ix+1,iy,iz) - __VELV(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (4 == comp)
//                     {
//                         const Real v0 = cy[iy]*(__VELV(ix,iy,iz) - __VELV(ix,iy-1,iz));
//                         const Real v1 = cy[iy+1]*(__VELV(ix,iy+1,iz) - __VELV(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (5 == comp)
//                     {
//                         const Real v0 = cz[iz]*(__VELV(ix,iy,iz) - __VELV(ix,iy,iz-1));
//                         const Real v1 = cz[iz+1]*(__VELV(ix,iy,iz+1) - __VELV(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (6 == comp)
//                     {
//                         const Real v0 = cx[ix]*(__VELW(ix,iy,iz) - __VELW(ix-1,iy,iz));
//                         const Real v1 = cx[ix+1]*(__VELW(ix+1,iy,iz) - __VELW(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (7 == comp)
//                     {
//                         const Real v0 = cy[iy]*(__VELW(ix,iy,iz) - __VELW(ix,iy-1,iz));
//                         const Real v1 = cy[iy+1]*(__VELW(ix,iy+1,iz) - __VELW(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (8 == comp)
//                     {
//                         const Real v0 = cz[iz]*(__VELW(ix,iy,iz) - __VELW(ix,iy,iz-1));
//                         const Real v1 = cz[iz+1]*(__VELW(ix,iy,iz+1) - __VELW(ix,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
// #else
//                     if (0 == comp)
//                     {
//                         const Real v0 = cx[ix]*(-__VELU(ix+1,iy,iz) + (Real)27.0*(__VELU(ix,iy,iz) - __VELU(ix-1,iy,iz)) + __VELU(ix-2,iy,iz));
//                         const Real v1 = cx[ix+1]*(-__VELU(ix+2,iy,iz) + (Real)27.0*(__VELU(ix+1,iy,iz) - __VELU(ix,iy,iz)) + __VELU(ix-1,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (1 == comp)
//                     {
//                         const Real v0 = cy[iy]*(-__VELU(ix,iy+1,iz) + (Real)27.0*(__VELU(ix,iy,iz) - __VELU(ix,iy-1,iz)) + __VELU(ix,iy-2,iz));
//                         const Real v1 = cy[iy+1]*(-__VELU(ix,iy+2,iz) + (Real)27.0*(__VELU(ix,iy+1,iz) - __VELU(ix,iy,iz)) + __VELU(ix,iy-1,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (2 == comp)
//                     {
//                         const Real v0 = cz[iz]*(-__VELU(ix,iy,iz+1) + (Real)27.0*(__VELU(ix,iy,iz) - __VELU(ix,iy,iz-1)) + __VELU(ix,iy,iz-2));
//                         const Real v1 = cz[iz+1]*(-__VELU(ix,iy,iz+2) + (Real)27.0*(__VELU(ix,iy,iz+1) - __VELU(ix,iy,iz)) + __VELU(ix,iy,iz-1));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (3 == comp)
//                     {
//                         const Real v0 = cx[ix]*(-__VELV(ix+1,iy,iz) + (Real)27.0*(__VELV(ix,iy,iz) - __VELV(ix-1,iy,iz)) + __VELV(ix-2,iy,iz));
//                         const Real v1 = cx[ix+1]*(-__VELV(ix+2,iy,iz) + (Real)27.0*(__VELV(ix+1,iy,iz) - __VELV(ix,iy,iz)) + __VELV(ix-1,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (4 == comp)
//                     {
//                         const Real v0 = cy[iy]*(-__VELV(ix,iy+1,iz) + (Real)27.0*(__VELV(ix,iy,iz) - __VELV(ix,iy-1,iz)) + __VELV(ix,iy-2,iz));
//                         const Real v1 = cy[iy+1]*(-__VELV(ix,iy+2,iz) + (Real)27.0*(__VELV(ix,iy+1,iz) - __VELV(ix,iy,iz)) + __VELV(ix,iy-1,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (5 == comp)
//                     {
//                         const Real v0 = cz[iz]*(-__VELV(ix,iy,iz+1) + (Real)27.0*(__VELV(ix,iy,iz) - __VELV(ix,iy,iz-1)) + __VELV(ix,iy,iz-2));
//                         const Real v1 = cz[iz+1]*(-__VELV(ix,iy,iz+2) + (Real)27.0*(__VELV(ix,iy,iz+1) - __VELV(ix,iy,iz)) + __VELV(ix,iy,iz-1));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (6 == comp)
//                     {
//                         const Real v0 = cx[ix]*(-__VELW(ix+1,iy,iz) + (Real)27.0*(__VELW(ix,iy,iz) - __VELW(ix-1,iy,iz)) + __VELW(ix-2,iy,iz));
//                         const Real v1 = cx[ix+1]*(-__VELW(ix+2,iy,iz) + (Real)27.0*(__VELW(ix+1,iy,iz) - __VELW(ix,iy,iz)) + __VELW(ix-1,iy,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (7 == comp)
//                     {
//                         const Real v0 = cy[iy]*(-__VELW(ix,iy+1,iz) + (Real)27.0*(__VELW(ix,iy,iz) - __VELW(ix,iy-1,iz)) + __VELW(ix,iy-2,iz));
//                         const Real v1 = cy[iy+1]*(-__VELW(ix,iy+2,iz) + (Real)27.0*(__VELW(ix,iy+1,iz) - __VELW(ix,iy,iz)) + __VELW(ix,iy-1,iz));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
//                     else if (8 == comp)
//                     {
//                         const Real v0 = cz[iz]*(-__VELW(ix,iy,iz+1) + (Real)27.0*(__VELW(ix,iy,iz) - __VELW(ix,iy,iz-1)) + __VELW(ix,iy,iz-2));
//                         const Real v1 = cz[iz+1]*(-__VELW(ix,iy,iz+2) + (Real)27.0*(__VELW(ix,iy,iz+1) - __VELW(ix,iy,iz)) + __VELW(ix,iy,iz-1));
//                         o(ix,iy,iz).dummy = (Real)0.5*(v0+v1);
//                     }
// #endif /* _SECONDORDER_FD_CORE_ */
                }
    }
};


struct Operator_Ucontraction
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_Ucontraction(): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
    Operator_Ucontraction(const Operator_Ucontraction& c): stencil(-1,-1,-1,2,2,2, false, 5, 0,1,2,3,4) {}
#else
    Operator_Ucontraction(): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
    Operator_Ucontraction(const Operator_Ucontraction& c): stencil(-2,-2,-2,3,3,3, false, 5, 0,1,2,3,4) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        // scratch pads
        Real sx[3];
        Real sy[TBlock::sizeX][3];
        Real sz[TBlock::sizeY][TBlock::sizeX][3];

        for(int iy=0; iy<TBlock::sizeY; iy++)
            for(int ix=0; ix<TBlock::sizeX; ix++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sz[iy][ix][0] = __FD_2ND(0,cz,__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1));
                sz[iy][ix][1] = __FD_2ND(0,cz,__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1));
                sz[iy][ix][2] = __FD_2ND(0,cz,__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1));
#else
                sz[iy][ix][0] = __FD_4TH(0,cz,__VELU(ix,iy,-2),__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1),__VELU(ix,iy,2));
                sz[iy][ix][1] = __FD_4TH(0,cz,__VELV(ix,iy,-2),__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1),__VELV(ix,iy,2));
                sz[iy][ix][2] = __FD_4TH(0,cz,__VELW(ix,iy,-2),__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1),__VELW(ix,iy,2));
#endif /* _SECONDORDER_FD_CORE_ */
            }

        for(int iz=0; iz<TBlock::sizeZ; iz++)
        {
            for (int ix = 0; ix < TBlock::sizeX; ++ix)
            {
#ifdef _SECONDORDER_FD_CORE_
                sy[ix][0] = __FD_2ND(0,cy,__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz));
                sy[ix][1] = __FD_2ND(0,cy,__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz));
                sy[ix][2] = __FD_2ND(0,cy,__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz));
#else
                sy[ix][0] = __FD_4TH(0,cy,__VELU(ix,-2,iz),__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz),__VELU(ix,2,iz));
                sy[ix][1] = __FD_4TH(0,cy,__VELV(ix,-2,iz),__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz),__VELV(ix,2,iz));
                sy[ix][2] = __FD_4TH(0,cy,__VELW(ix,-2,iz),__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz),__VELW(ix,2,iz));
#endif /* _SECONDORDER_FD_CORE_ */
            }

            for(int iy=0; iy<TBlock::sizeY; iy++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sx[0] = __FD_2ND(0,cx,__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz));
                sx[1] = __FD_2ND(0,cx,__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz));
                sx[2] = __FD_2ND(0,cx,__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz));
#else
                sx[0] = __FD_4TH(0,cx,__VELU(-2,iy,iz),__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz),__VELU(2,iy,iz));
                sx[1] = __FD_4TH(0,cx,__VELV(-2,iy,iz),__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz),__VELV(2,iy,iz));
                sx[2] = __FD_4TH(0,cx,__VELW(-2,iy,iz),__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz),__VELW(2,iy,iz));
#endif /* _SECONDORDER_FD_CORE_ */

                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
#ifdef _SECONDORDER_FD_CORE_
                    const Real ux1 = __FD_2ND(ix+1,cx,__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz));
                    const Real uy1 = __FD_2ND(iy+1,cy,__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz));
                    const Real uz1 = __FD_2ND(iz+1,cz,__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1));
                    const Real vx1 = __FD_2ND(ix+1,cx,__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz));
                    const Real vy1 = __FD_2ND(iy+1,cy,__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz));
                    const Real vz1 = __FD_2ND(iz+1,cz,__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1));
                    const Real wx1 = __FD_2ND(ix+1,cx,__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz));
                    const Real wy1 = __FD_2ND(iy+1,cy,__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz));
                    const Real wz1 = __FD_2ND(iz+1,cz,__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1));
#else
                    const Real ux1 = __FD_4TH(ix+1,cx,__VELU(ix-2,iy,iz),__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz),__VELU(ix+2,iy,iz));
                    const Real uy1 = __FD_4TH(iy+1,cy,__VELU(ix,iy-2,iz),__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz),__VELU(ix,iy+2,iz));
                    const Real uz1 = __FD_4TH(iz+1,cz,__VELU(ix,iy,iz-2),__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1),__VELU(ix,iy,iz+2));
                    const Real vx1 = __FD_4TH(ix+1,cx,__VELV(ix-2,iy,iz),__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz),__VELV(ix+2,iy,iz));
                    const Real vy1 = __FD_4TH(iy+1,cy,__VELV(ix,iy-2,iz),__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz),__VELV(ix,iy+2,iz));
                    const Real vz1 = __FD_4TH(iz+1,cz,__VELV(ix,iy,iz-2),__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1),__VELV(ix,iy,iz+2));
                    const Real wx1 = __FD_4TH(ix+1,cx,__VELW(ix-2,iy,iz),__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz),__VELW(ix+2,iy,iz));
                    const Real wy1 = __FD_4TH(iy+1,cy,__VELW(ix,iy-2,iz),__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz),__VELW(ix,iy+2,iz));
                    const Real wz1 = __FD_4TH(iz+1,cz,__VELW(ix,iy,iz-2),__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1),__VELW(ix,iy,iz+2));
#endif /* _SECONDORDER_FD_CORE_ */

                    const Real dudx = (sx[0]+ux1);
                    const Real dudy = (sy[ix][0]+uy1);
                    const Real dudz = (sz[ix][iy][0]+uz1);
                    const Real dvdx = (sx[1]+vx1);
                    const Real dvdy = (sy[ix][1]+vy1);
                    const Real dvdz = (sz[iy][ix][1]+vz1);
                    const Real dwdx = (sx[2]+wx1);
                    const Real dwdy = (sy[ix][2]+wy1);
                    const Real dwdz = (sz[iy][ix][2]+wz1);

                    o(ix,iy,iz).dummy = static_cast<Real>(-0.25)*(dudx*dudx+dvdy*dvdy+dwdz*dwdz + static_cast<Real>(2.0)*(dudy*dvdx+dudz*dwdx+dvdz*dwdy)); // 0.25 from averaging

                    sx[0] = ux1;
                    sy[ix][0] = uy1;
                    sz[iy][ix][0] = uz1;
                    sx[1] = vx1;
                    sy[ix][1] = vy1;
                    sz[iy][ix][1] = vz1;
                    sx[2] = wx1;
                    sy[ix][2] = wy1;
                    sz[iy][ix][2] = wz1;
                }
            }
        }
    }
};


struct Operator_PIC
{
    static const std::string NAME;
    static const int CLASS = 1;
    static const unsigned int MASK  = 0x1;
    static const unsigned int INVALID = 0x0;
    StencilInfo stencil;

#ifdef _SECONDORDER_FD_CORE_
    Operator_PIC(): stencil(-1,-1,-1,2,2,2, false, 3, 2,3,4) {}
    Operator_PIC(const Operator_PIC& c): stencil(-1,-1,-1,2,2,2, false, 3, 2,3,4) {}
#else
    Operator_PIC(): stencil(-2,-2,-2,3,3,3, false, 3, 2,3,4) {}
    Operator_PIC(const Operator_PIC& c): stencil(-2,-2,-2,3,3,3, false, 3, 2,3,4) {}
#endif /* _SECONDORDER_FD_CORE_ */

    template<typename TLab, typename TBlock>
    void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        const double gamma = Simulation_Environment::GAMMA1;
        const double factor = (gamma-1.0)/gamma;

        const CoeffsFD_Block_t& cx = o.coeffsFD_x;
        const CoeffsFD_Block_t& cy = o.coeffsFD_y;
        const CoeffsFD_Block_t& cz = o.coeffsFD_z;

        // scratch pads
        Real sx[3];
        Real sy[TBlock::sizeX][3];
        Real sz[TBlock::sizeY][TBlock::sizeX][3];

        for(int iy=0; iy<TBlock::sizeY; iy++)
            for(int ix=0; ix<TBlock::sizeX; ix++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sz[iy][ix][0] = __FD_2ND(0,cz,__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1));
                sz[iy][ix][1] = __FD_2ND(0,cz,__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1));
                sz[iy][ix][2] = __FD_2ND(0,cz,__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1));
#else
                sz[iy][ix][0] = __FD_4TH(0,cz,__VELU(ix,iy,-2),__VELU(ix,iy,-1),__VELU(ix,iy,0),__VELU(ix,iy,1),__VELU(ix,iy,2));
                sz[iy][ix][1] = __FD_4TH(0,cz,__VELV(ix,iy,-2),__VELV(ix,iy,-1),__VELV(ix,iy,0),__VELV(ix,iy,1),__VELV(ix,iy,2));
                sz[iy][ix][2] = __FD_4TH(0,cz,__VELW(ix,iy,-2),__VELW(ix,iy,-1),__VELW(ix,iy,0),__VELW(ix,iy,1),__VELW(ix,iy,2));
#endif /* _SECONDORDER_FD_CORE_ */
            }

        for(int iz=0; iz<TBlock::sizeZ; iz++)
        {
            for (int ix = 0; ix < TBlock::sizeX; ++ix)
            {
#ifdef _SECONDORDER_FD_CORE_
                sy[ix][0] = __FD_2ND(0,cy,__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz));
                sy[ix][1] = __FD_2ND(0,cy,__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz));
                sy[ix][2] = __FD_2ND(0,cy,__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz));
#else
                sy[ix][0] = __FD_4TH(0,cy,__VELU(ix,-2,iz),__VELU(ix,-1,iz),__VELU(ix,0,iz),__VELU(ix,1,iz),__VELU(ix,2,iz));
                sy[ix][1] = __FD_4TH(0,cy,__VELV(ix,-2,iz),__VELV(ix,-1,iz),__VELV(ix,0,iz),__VELV(ix,1,iz),__VELV(ix,2,iz));
                sy[ix][2] = __FD_4TH(0,cy,__VELW(ix,-2,iz),__VELW(ix,-1,iz),__VELW(ix,0,iz),__VELW(ix,1,iz),__VELW(ix,2,iz));
#endif /* _SECONDORDER_FD_CORE_ */
            }

            for(int iy=0; iy<TBlock::sizeY; iy++)
            {
#ifdef _SECONDORDER_FD_CORE_
                sx[0] = __FD_2ND(0,cx,__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz));
                sx[1] = __FD_2ND(0,cx,__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz));
                sx[2] = __FD_2ND(0,cx,__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz));
#else
                sx[0] = __FD_4TH(0,cx,__VELU(-2,iy,iz),__VELU(-1,iy,iz),__VELU(0,iy,iz),__VELU(1,iy,iz),__VELU(2,iy,iz));
                sx[1] = __FD_4TH(0,cx,__VELV(-2,iy,iz),__VELV(-1,iy,iz),__VELV(0,iy,iz),__VELV(1,iy,iz),__VELV(2,iy,iz));
                sx[2] = __FD_4TH(0,cx,__VELW(-2,iy,iz),__VELW(-1,iy,iz),__VELW(0,iy,iz),__VELW(1,iy,iz),__VELW(2,iy,iz));
#endif /* _SECONDORDER_FD_CORE_ */

                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
#ifdef _SECONDORDER_FD_CORE_
                    const Real ux1 = __FD_2ND(ix+1,cx,__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz));
                    const Real uy1 = __FD_2ND(iy+1,cy,__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz));
                    const Real uz1 = __FD_2ND(iz+1,cz,__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1));
                    const Real vx1 = __FD_2ND(ix+1,cx,__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz));
                    const Real vy1 = __FD_2ND(iy+1,cy,__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz));
                    const Real vz1 = __FD_2ND(iz+1,cz,__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1));
                    const Real wx1 = __FD_2ND(ix+1,cx,__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz));
                    const Real wy1 = __FD_2ND(iy+1,cy,__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz));
                    const Real wz1 = __FD_2ND(iz+1,cz,__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1));
#else
                    const Real ux1 = __FD_4TH(ix+1,cx,__VELU(ix-2,iy,iz),__VELU(ix-1,iy,iz),__VELU(ix,iy,iz),__VELU(ix+1,iy,iz),__VELU(ix+2,iy,iz));
                    const Real uy1 = __FD_4TH(iy+1,cy,__VELU(ix,iy-2,iz),__VELU(ix,iy-1,iz),__VELU(ix,iy,iz),__VELU(ix,iy+1,iz),__VELU(ix,iy+2,iz));
                    const Real uz1 = __FD_4TH(iz+1,cz,__VELU(ix,iy,iz-2),__VELU(ix,iy,iz-1),__VELU(ix,iy,iz),__VELU(ix,iy,iz+1),__VELU(ix,iy,iz+2));
                    const Real vx1 = __FD_4TH(ix+1,cx,__VELV(ix-2,iy,iz),__VELV(ix-1,iy,iz),__VELV(ix,iy,iz),__VELV(ix+1,iy,iz),__VELV(ix+2,iy,iz));
                    const Real vy1 = __FD_4TH(iy+1,cy,__VELV(ix,iy-2,iz),__VELV(ix,iy-1,iz),__VELV(ix,iy,iz),__VELV(ix,iy+1,iz),__VELV(ix,iy+2,iz));
                    const Real vz1 = __FD_4TH(iz+1,cz,__VELV(ix,iy,iz-2),__VELV(ix,iy,iz-1),__VELV(ix,iy,iz),__VELV(ix,iy,iz+1),__VELV(ix,iy,iz+2));
                    const Real wx1 = __FD_4TH(ix+1,cx,__VELW(ix-2,iy,iz),__VELW(ix-1,iy,iz),__VELW(ix,iy,iz),__VELW(ix+1,iy,iz),__VELW(ix+2,iy,iz));
                    const Real wy1 = __FD_4TH(iy+1,cy,__VELW(ix,iy-2,iz),__VELW(ix,iy-1,iz),__VELW(ix,iy,iz),__VELW(ix,iy+1,iz),__VELW(ix,iy+2,iz));
                    const Real wz1 = __FD_4TH(iz+1,cz,__VELW(ix,iy,iz-2),__VELW(ix,iy,iz-1),__VELW(ix,iy,iz),__VELW(ix,iy,iz+1),__VELW(ix,iy,iz+2));
#endif /* _SECONDORDER_FD_CORE_ */

                    const Real dudx = (sx[0]+ux1);
                    const Real dudy = (sy[ix][0]+uy1);
                    const Real dudz = (sz[ix][iy][0]+uz1);
                    const Real dvdx = (sx[1]+vx1);
                    const Real dvdy = (sy[ix][1]+vy1);
                    const Real dvdz = (sz[iy][ix][1]+vz1);
                    const Real dwdx = (sx[2]+wx1);
                    const Real dwdy = (sy[ix][2]+wy1);
                    const Real dwdz = (sz[iy][ix][2]+wz1);

                    o(ix,iy,iz).dummy = static_cast<Real>(-0.25)*factor*(dudx*dudx+dvdy*dvdy+dwdz*dwdz + static_cast<Real>(2.0)*(dudy*dvdx+dudz*dwdx+dvdz*dwdy)); // 0.25 from averaging

                    sx[0] = ux1;
                    sy[ix][0] = uy1;
                    sz[iy][ix][0] = uz1;
                    sx[1] = vx1;
                    sy[ix][1] = vy1;
                    sz[iy][ix][1] = vz1;
                    sx[2] = wx1;
                    sy[ix][2] = wy1;
                    sz[iy][ix][2] = wz1;
                }
            }
        }
    }
};

// abbreviations
///////////////////////////////////////////////////////////////////////////////
typedef Operator_Vort         OVort;
typedef Operator_IgradA2I     OIgradA2I;
typedef Operator_divU         OdivU;
typedef Operator_SSDeviatoric OSSd;
typedef Operator_Qcrit        OQcrit;
typedef Operator_gradU        OgradU;
typedef Operator_gradUij<0>   OgradU_0;
typedef Operator_gradUij<1>   OgradU_1;
typedef Operator_gradUij<2>   OgradU_2;
typedef Operator_gradUij<3>   OgradU_3;
typedef Operator_gradUij<4>   OgradU_4;
typedef Operator_gradUij<5>   OgradU_5;
typedef Operator_gradUij<6>   OgradU_6;
typedef Operator_gradUij<7>   OgradU_7;
typedef Operator_gradUij<8>   OgradU_8;
typedef Operator_Ucontraction OUcon;
typedef Operator_PIC          OPIC;

#undef __VELU
#undef __VELV
#undef __VELW
#undef __ALPHA
#undef __FD_2ND
#undef __FD_4TH

#endif /* VECTOROPERATORNONUNIFORM_H_U4LE7NP8 */
