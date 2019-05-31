/*
 *  Source_CPP.h
 *  MPCFcore
 *
 *  Created by Ursula Rasthofer
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "common.h"
#include <Cubism/BlockInfo.h>
using namespace cubism;

class Source_CPP
{
protected:

    Real a, dt;
    Real g[3], f[3];

    inline bool _is_aligned(const void * const ptr, unsigned int alignment) const
    {
        return ((size_t)ptr) % alignment == 0;
    }

public:

    Source_CPP(Real a, Real dt, const Real gaccel[3], const Real volforce[3]): a(a), dt(dt)
    {
        g[0]=gaccel[0];
        g[1]=gaccel[1];
        g[2]=gaccel[2];
        f[0]=volforce[0];
        f[1]=volforce[1];
        f[2]=volforce[2];
    }

    void compute (const Real * const srcfirst, const int srcfloats, Real * const dstfirst, const int dstfloats) const;

    static void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const size_t NCORES, const size_t NT, const size_t NBLOCKS, const float MEASUREDTIME, const bool bAwk=false)
    {}
};


class OneWayAcousticSource_CPP
{
public:

    struct SourceParameter
    {
        Real gamma, smooth;
        Real t0, c0, x0;
        Real amplitude, omega, T;
        Real sigma_t;
    };

    OneWayAcousticSource_CPP(const Real t, const Real dt, const SourceParameter& p):
        t(t), dt(dt),
        // reg0(1.0/(std::sqrt(2.0*M_PI)*6.0*p.smooth)),
        // reg1(1.0/(2.0*36.0*p.smooth*p.smooth)),
        reg0(1.0/(std::sqrt(2.0*M_PI)*p.smooth)),
        reg1(1.0/(2.0*p.smooth*p.smooth)),
        param(p)
    {}

    void compute(const Real * const srcfirst, const int srcfloats, Real * const dstfirst, const int dstfloats, const BlockInfo& info) const
    {
        assert(srcfloats==dstfloats);

        const Real f = _f(t);

        for (int iz = 0; iz < _BLOCKSIZE_; ++iz)
            for (int iy = 0; iy < _BLOCKSIZE_; ++iy)
                for (int ix = 0; ix < _BLOCKSIZE_; ++ix)
                {
                    const int i = dstfloats*(ix + _BLOCKSIZE_*(iy + _BLOCKSIZE_*iz));
                    Real pos[3];
                    info.pos(pos, ix, iy, iz);

                    // if (std::abs(pos[0] - param.x0) < 6.0*param.smooth)
                    {
                        const Real reg = reg0*std::exp(-(pos[0]-param.x0)*(pos[0]-param.x0)*reg1);

                        const Real s2 = f*reg;
                        const Real s1 = s2/param.c0;
                        const Real s3 = param.c0*param.c0/(param.gamma - 1.0)*s1;

                        dstfirst[i+0] += dt*s1; // mass 1
                        dstfirst[i+2] += dt*s2; // x-momentum
                        dstfirst[i+5] += dt*s3; // energy
                    }
                }
    }

protected:

    const Real t;
    const Real dt;
    const Real reg0;
    const Real reg1;
    SourceParameter param;

    inline bool _is_aligned(const void * const ptr, unsigned int alignment) const
    {
        return ((size_t)ptr) % alignment == 0;
    }

    // source function
    inline Real _f(const Real t) const
    {
        const Real fac0 = param.amplitude/(std::sqrt(2.0*M_PI)*param.sigma_t);
        const Real fac1 = 1.0/(2.0*param.sigma_t*param.sigma_t);
        return fac0*std::exp(-(t-param.t0)*(t-param.t0)*fac1);
        // return (t>=param.t0 && t<=param.T+param.t0) ? param.amplitude*std::sin(param.omega*(t-param.t0)) : 0.0;
    }
};
