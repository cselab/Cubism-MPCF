/*
 *  Streamer.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 09/13/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef STREAMER_H_0HQKKPXE
#define STREAMER_H_0HQKKPXE

#include <cassert>
#include <cmath>
#include <string>
#include <algorithm>
#include "Types.h"
#include "NonUniform.h"

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
// class 14:  writes into tmp[iz][iy][ix][0-2] (three-component vector) (MASK = 0x2000; INVALID = 0x1df1)
// class 15:  writes into tmp[iz][iy][ix][3-5] (three-component vector) (MASK = 0x4000; INVALID = 0x3c8f)
// class 16:  writes into dummy + tmp[iz][iy][ix][0-7] (nine-component rank-2 tensor) (MASK = 0x8000; INVALID = 0x0)

struct StreamerAlpha2
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
#ifndef _CONVERTCLIPDEBUGOUT_
        output[0] = std::max((Real)ALPHAEPS, std::min((Real)(1.0-ALPHAEPS),el.alpha2));
#else
        output[0] = el.alpha2;
#endif
#else
        output[0] = el.alpha2;
#endif
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
#ifndef _CONVERTCLIPDEBUGOUT_
        return std::max((Real)ALPHAEPS, std::min((Real)(1.0-ALPHAEPS),el.alpha2));
#else
        return el.alpha2;
#endif
#else
        return el.alpha2;
#endif
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerPressure
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif

        output[0] = Simulation_Environment::pressure_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif

        return Simulation_Environment::pressure_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerDensity
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        output[0] = el.alpha1rho1 + el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
        output[0] = r1 + r2;
#endif
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        return el.alpha1rho1 + el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
        return r1 + r2;
#endif
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerEnergy
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        output[0] = el.energy;
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        return el.energy;
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerQcriterion
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 1;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        output[0] = el.dummy;
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        return el.dummy;
    }

    static const char * getAttributeName() { return "Scalar"; }
};


template <int comp>
struct StreamerVelocity
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        assert(0 <= comp && comp < 3);
        const Real * const pel = (Real*)&b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = pel[0];
        const Real r2 = pel[1];
#else
        const Real r1 = std::max((Real)0.0, pel[0]);
        const Real r2 = std::max((Real)0.0, pel[1]);
#endif
        const Real rho = r1 + r2;

        output[0] = pel[2+comp]/rho;
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        assert(0 <= comp && comp < 3);
        const Real * const pel = (Real*)&b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = pel[0];
        const Real r2 = pel[1];
#else
        const Real r1 = std::max((Real)0.0, pel[0]);
        const Real r2 = std::max((Real)0.0, pel[1]);
#endif
        const Real rho = r1 + r2;

        return pel[2+comp]/rho;
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerVelocityVector
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 3;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif
        const Real rho = r1 + r2;

        output[0] = el.ru/rho;
        output[1] = el.rv/rho;
        output[2] = el.rw/rho;
    }

    static const char * getAttributeName() { return "Vector"; }
};


struct StreamerVelocityMagnitude
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif
        const Real rho = r1 + r2;

        const Real u = el.ru/rho;
        const Real v = el.rv/rho;
        const Real w = el.rw/rho;

        output[0] = std::sqrt(u*u + v*v + w*w);
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif
        const Real rho = r1 + r2;

        const Real u = el.ru/rho;
        const Real v = el.rv/rho;
        const Real w = el.rw/rho;

        return std::sqrt(u*u + v*v + w*w);
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerGradU
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 9;
    static const int CLASS = 16;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

        output[0] = el.dummy;
        output[1] = b.tmp[iz][iy][ix][0];
        output[2] = b.tmp[iz][iy][ix][1];
        output[3] = b.tmp[iz][iy][ix][2];
        output[4] = b.tmp[iz][iy][ix][3];
        output[5] = b.tmp[iz][iy][ix][4];
        output[6] = b.tmp[iz][iy][ix][5];
        output[7] = b.tmp[iz][iy][ix][6];
        output[8] = b.tmp[iz][iy][ix][7];
    }

    static const char * getAttributeName() { return "Tensor"; }
};


template <int comp>
struct StreamerVorticity
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 14;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        assert(0 <= comp && comp < 3);
        output[0] = b.tmp[iz][iy][ix][comp];
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        assert(0 <= comp && comp < 3);
        return b.tmp[iz][iy][ix][comp];
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerVorticityVector
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 3;
    static const int CLASS = 14;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        output[0] = b.tmp[iz][iy][ix][0];
        output[1] = b.tmp[iz][iy][ix][1];
        output[2] = b.tmp[iz][iy][ix][2];
    }

    static const char * getAttributeName() { return "Vector" ; }
};


struct StreamerVorticityMagnitude
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 14;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        const Real wx = b.tmp[iz][iy][ix][0];
        const Real wy = b.tmp[iz][iy][ix][1];
        const Real wz = b.tmp[iz][iy][ix][2];
        output[0] = std::sqrt(wx*wx + wy*wy + wz*wz);
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        const Real wx = b.tmp[iz][iy][ix][0];
        const Real wy = b.tmp[iz][iy][ix][1];
        const Real wz = b.tmp[iz][iy][ix][2];
        return std::sqrt(wx*wx + wy*wy + wz*wz);
    }

    static const char * getAttributeName() { return "Scalar" ; }
};


struct StreamerAllConservative
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 9;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
#ifndef _CONVERTCLIPDEBUGOUT_
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#else
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#endif
#endif
        output[0] = r1;
        output[1] = r2;
        output[2] = el.ru;
        output[3] = el.rv;
        output[4] = el.rw;
        output[5] = el.energy;
#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
#ifndef _CONVERTCLIPDEBUGOUT_
        output[6] = std::max((Real)ALPHAEPS, std::min((Real)(1.0-ALPHAEPS),el.alpha2));
#else
        output[6] = el.alpha2;
#endif
#else
        output[6] = el.alpha2;
#endif
        output[7] = el.alpha2; // non-clipped forever
        output[8] = el.dummy;
    }

    static const char * getAttributeName() { return "Tensor"; }
};


struct StreamerAllPrimitive
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 9;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
#ifndef _CONVERTCLIPDEBUGOUT_
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#else
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#endif
#endif
        output[0] = r1;
        output[1] = r2;
        const Real rho = r1 + r2;
        output[2] = rho;
        output[3] = el.ru/rho;
        output[4] = el.rv/rho;
        output[5] = el.rw/rho;
        output[6] = el.energy;
#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
#ifndef _CONVERTCLIPDEBUGOUT_
        output[7] = std::max((Real)ALPHAEPS, std::min((Real)(1.0-ALPHAEPS),el.alpha2));
#else
        output[7] = el.alpha2;
#endif
#else
        output[7] = el.alpha2;
#endif
        output[8] = Simulation_Environment::pressure_5eq(rho, el.energy, el.ru, el.rv, el.rw, el.alpha2);
    }

    static const char * getAttributeName() { return "Tensor"; }
};


struct StreamerK
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
#ifndef _NOK_
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif

        const Real pressure = Simulation_Environment::pressure_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);
#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
        const Real alpha2   = std::max((Real)ALPHAEPS, std::min((Real)(1.0-ALPHAEPS),el.alpha2));
#else
        const Real alpha2   = el.alpha2;
#endif
        const Real alpha1   = 1.0 - alpha2;
        const Real rc1sq    = Simulation_Environment::GAMMA1 * (pressure + Simulation_Environment::PC1);
        const Real rc2sq    = Simulation_Environment::GAMMA2 * (pressure + Simulation_Environment::PC2);

        output[0] = alpha1 * alpha2 * (rc1sq - rc2sq) / (alpha1 * rc2sq + alpha2 * rc1sq);
#else
        output[0] = 0.0;
#endif /* _NOK_ */
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
#ifndef _NOK_
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif

        const Real pressure = Simulation_Environment::pressure_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);
#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
        const Real alpha2   = std::max((Real)ALPHAEPS, std::min((Real)(1.0-ALPHAEPS),el.alpha2));
#else
        const Real alpha2   = el.alpha2;
#endif
        const Real alpha1   = 1.0 - alpha2;
        const Real rc1sq    = Simulation_Environment::GAMMA1 * (pressure + Simulation_Environment::PC1);
        const Real rc2sq    = Simulation_Environment::GAMMA2 * (pressure + Simulation_Environment::PC2);

        return alpha1 * alpha2 * (rc1sq - rc2sq) / (alpha1 * rc2sq + alpha2 * rc1sq);
#else
        return 0.0;
#endif /* _NOK_ */
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerDivU
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 5;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        output[0] = b.tmp[iz][iy][ix][3];
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        return b.tmp[iz][iy][ix][3];
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerKDivU
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 5;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
#ifndef _NOK_
        typedef typename TBlock::ElementType TElement;

        const Real divu = b.tmp[iz][iy][ix][3];
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif

        const Real pressure = Simulation_Environment::pressure_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);

#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
        const Real alpha2   = std::max((Real)ALPHAEPS, std::min((Real)(1.0-ALPHAEPS),el.alpha2));
#else
        const Real alpha2   = el.alpha2;
#endif
        const Real alpha1   = 1.0 - alpha2;
        const Real rc1sq    = Simulation_Environment::GAMMA1 * (pressure + Simulation_Environment::PC1);
        const Real rc2sq    = Simulation_Environment::GAMMA2 * (pressure + Simulation_Environment::PC2);

        output[0] = divu * alpha1 * alpha2 * (rc1sq - rc2sq) / (alpha1 * rc2sq + alpha2 * rc1sq);
#else
        output[0] = 0.0;
#endif /* _NOK_ */
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
#ifndef _NOK_
        typedef typename TBlock::ElementType TElement;

        const Real divu = b.tmp[iz][iy][ix][3];
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif

        const Real pressure = Simulation_Environment::pressure_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);

#if defined(_ALPHACLIP_) || defined(_CONVERTCLIP_)
        const Real alpha2   = std::max((Real)ALPHAEPS, std::min((Real)(1.0-ALPHAEPS),el.alpha2));
#else
        const Real alpha2   = el.alpha2;
#endif
        const Real alpha1   = 1.0 - alpha2;
        const Real rc1sq    = Simulation_Environment::GAMMA1 * (pressure + Simulation_Environment::PC1);
        const Real rc2sq    = Simulation_Environment::GAMMA2 * (pressure + Simulation_Environment::PC2);

        return divu * alpha1 * alpha2 * (rc1sq - rc2sq) / (alpha1 * rc2sq + alpha2 * rc1sq);
#else
        return 0.0;
#endif /* _NOK_ */
    }

    static const char * getAttributeName() { return "Scalar"; }
};


template <int comp>
struct StreamerAoScomponent
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        assert(0 <= comp && comp < sizeof(FluidElement)/sizeof(Real));
        const Real * const pel = (Real*)&b(ix,iy,iz);
        output[0] = pel[comp];
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        assert(0 <= comp && comp < sizeof(FluidElement)/sizeof(Real));
        const Real * const pel = (Real*)&b(ix,iy,iz);
        return pel[comp];
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerMach
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif
        const Real c = Simulation_Environment::speedOfSound_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);

        const Real rInv = 1.0/(r1+r2);
        const Real u = el.ru * rInv;
        const Real v = el.rv * rInv;
        const Real w = el.rw * rInv;

        const Real IUI = std::sqrt(u*u + v*v + w*w);

        output[0] = IUI/c;
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif
        const Real c = Simulation_Environment::speedOfSound_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);

        const Real rInv = 1.0/(r1+r2);
        const Real u = el.ru * rInv;
        const Real v = el.rv * rInv;
        const Real w = el.rw * rInv;

        const Real IUI = std::sqrt(u*u + v*v + w*w);

        return IUI/c;
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerSpeedOfSound
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif
        output[0] = Simulation_Environment::speedOfSound_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

#ifndef _CONVERTCLIP_
        const Real r1 = el.alpha1rho1;
        const Real r2 = el.alpha2rho2;
#else
        const Real r1 = std::max((Real)0.0,el.alpha1rho1);
        const Real r2 = std::max((Real)0.0,el.alpha2rho2);
#endif
        return Simulation_Environment::speedOfSound_5eq(
                r1 + r2,
                el.energy,
                el.ru, el.rv, el.rw,
                el.alpha2);
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerDummy
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 1;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        output[0] = el.dummy;
    }

    template <typename TBlock, typename TReal>
    static inline void operate(TBlock& b, const TReal input[NCHANNELS], const int ix, const int iy, const int iz)
    {
        b(ix,iy,iz).dummy = input[0];
    }

    template <typename TBlock>
    static inline Real operate(const TBlock& b, const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        return el.dummy;
    }

    static const char * getAttributeName() { return "Scalar"; }
};


struct StreamerSerialization
{
    // Streamer is used for serialization purposes
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 7;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& input = b(ix,iy,iz);

        output[0] = input.alpha1rho1;
        output[1] = input.alpha2rho2;
        output[2] = input.ru;
        output[3] = input.rv;
        output[4] = input.rw;
        output[5] = input.energy;
        output[6] = input.alpha2;
    }

    template <typename TBlock, typename TReal>
    static inline void operate(TBlock& b, const TReal input[NCHANNELS], const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        TElement& el = b(ix,iy,iz);

        el.alpha1rho1 = input[0];
        el.alpha2rho2 = input[1];
        el.ru         = input[2];
        el.rv         = input[3];
        el.rw         = input[4];
        el.energy     = input[5];
        el.alpha2     = input[6];
    }

    // TODO: (fabianw@mavt.ethz.ch; Tue 18 Oct 2016 04:32:32 PM CEST) these are
    // required by zbin dumper.  Unify this to one streamer
    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal *ovalue, const int field)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& input = b(ix,iy,iz);

        switch(field) {
            case 0: *ovalue = input.alpha1rho1; break;
            case 1: *ovalue = input.alpha2rho2; break;
            case 2: *ovalue = input.ru; break;
            case 3: *ovalue = input.rv; break;
            case 4: *ovalue = input.rw; break;
            case 5: *ovalue = input.energy; break;
            case 6: *ovalue = input.alpha2; break;
            default: printf("unknown field\n"); abort(); break;
        }
    }

    template <typename TBlock, typename TReal>
    static inline void operate(TBlock& b, const TReal ivalue, const int ix, const int iy, const int iz, const int field)
    {
        typedef typename TBlock::ElementType TElement;
        TElement& input = b(ix,iy,iz);

        switch(field) {
            case 0:  input.alpha1rho1     = ivalue; break;
            case 1:  input.alpha2rho2     = ivalue; break;
            case 2:  input.ru             = ivalue; break;
            case 3:  input.rv             = ivalue; break;
            case 4:  input.rw             = ivalue; break;
            case 5:  input.energy         = ivalue; break;
            case 6:  input.alpha2         = ivalue; break;
            default: printf("unknown field\n"); abort(); break;
        }
    }

    static const char * getAttributeName() { return "Serialization"; }
};

#endif /* STREAMER_H_0HQKKPXE */
