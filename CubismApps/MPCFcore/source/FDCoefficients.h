/*
 *  FDCoefficients.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 07/24/2017
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef FDCOEFFICIENTS_H_DLHT8AXM
#define FDCOEFFICIENTS_H_DLHT8AXM

#include <iostream>
#include <cassert>
#include <cstdlib>
#include "common.h"

///////////////////////////////////////////////////////////////////////////////
// Generate static (compile time index array) using recursion
// https://stackoverflow.com/questions/2978259/programmatically-create-static-arrays-at-compile-time-in-c
// Purpose:
// The finite differences are evaluated at the cell faces.  The coefficients
// for the first face in cell i=0 are computed for x_{-1/2} while the remainig
// coefficients are computed for faces at x_{i+1/2}.  The StaticIndexArray is
// required to index the correct cell data when the FD scheme is evaluated.
// This keeps the ghosts symmetric.  If all coefficients were evaluated at
// x_{i+1/2}, 2 ghosts would be required at the left boundary, while only 1
// ghost is needed on the right boundary.
///////////////////////////////////////////////////////////////////////////////
namespace StaticIndexArray
{
    template<int... args> struct ArrayHolder {
        static const int idx[sizeof...(args)];
    };

    template<int... args>
    const int ArrayHolder<args...>::idx[sizeof...(args)] = { args... };

    template<int N, template<int> class F, int... args>
    struct generate_array_impl {
        typedef typename generate_array_impl<N-1, F, F<N>::value, args...>::staticArray staticArray;
    };

    template<template<int> class F, int... args>
    struct generate_array_impl<0, F, args...> {
        typedef ArrayHolder<0, args...> staticArray;
    };

    template<int N, template<int> class F>
    struct generate_array {
        typedef typename generate_array_impl<N-1, F>::staticArray staticArray;
    };

    template<int index> struct MetaFunc {
        enum { value = index-1 };
    };

    // assumes _BLOCKSIZE_ is the same for all spatial dimensions
    typedef generate_array<_BLOCKSIZE_+1, MetaFunc>::staticArray array;
} /* StaticIndexArray */

///////////////////////////////////////////////////////////////////////////////
// generic coefficient types
///////////////////////////////////////////////////////////////////////////////
// per block
template <int _NCOEFFS_>
struct CenteredCoeffsFD_Block_t
{
    static constexpr int N_COEFFS = _NCOEFFS_;
    Real m_coeffs[N_COEFFS][_BLOCKSIZE_+1];

    // accessing coefficients
    inline Real (&get_coeffs(const int i=0)) [_BLOCKSIZE_+1]
    {
        // get_coeffs(-1) -> coefficients for data \phi_{i-1}
        // get_coeffs( 0) -> coefficients for data \phi_{i} (center coefficient)
        // get_coeffs( 1) -> coefficients for data \phi_{i+1}
        // etc.
        assert(i >= -N_COEFFS/2 && i <= N_COEFFS/2);
        return m_coeffs[i+N_COEFFS/2];
    }

    inline const Real (&get_coeffs(const int i=0) const) [_BLOCKSIZE_+1]
    {
        assert(i >= -N_COEFFS/2 && i <= N_COEFFS/2);
        return m_coeffs[i+N_COEFFS/2];
    }
};


// low level
template <int _NCOEFFS_>
struct CenteredCoeffsFD_All_t
{
    static constexpr int N_COEFFS = _NCOEFFS_;
    Real* m_coeffs[N_COEFFS];

    CenteredCoeffsFD_All_t(const unsigned int N)
    {
        for (int i = 0; i < _NCOEFFS_; ++i)
            posix_memalign((void**)&m_coeffs[i], _ALIGNBYTES_, N*sizeof(Real));
    }

    ~CenteredCoeffsFD_All_t()
    {
        for (int i = 0; i < _NCOEFFS_; ++i)
            free(m_coeffs[i]);
    }

    // accessing pointers to full coefficient arrays (not only one block)
    inline Real* get_coeffs(const int i=0)
    {
        // get_coeffs(-1) -> coefficients for data \phi_{i-1}
        // get_coeffs( 0) -> coefficients for data \phi_{i} (center coefficient)
        // get_coeffs( 1) -> coefficients for data \phi_{i+1}
        // etc.
        assert(i >= -N_COEFFS/2 && i <= N_COEFFS/2);
        return m_coeffs[i+N_COEFFS/2];
    }

    inline const Real* get_coeffs(const int i=0) const
    {
        assert(i >= -N_COEFFS/2 && i <= N_COEFFS/2);
        return m_coeffs[i+N_COEFFS/2];
    }
};

// 2nd-order
typedef CenteredCoeffsFD_Block_t<3> CoeffsFD_Block_2nd_t;
typedef CenteredCoeffsFD_All_t<3>   CoeffsFD_All_2nd_t;

// 4th-order
typedef CenteredCoeffsFD_Block_t<5> CoeffsFD_Block_4th_t;
typedef CenteredCoeffsFD_All_t<5>   CoeffsFD_All_4th_t;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Finite difference classes:
// Coefficients are based on Bengt Fornberg "Generation of Finite Difference
// Formulas on Arbitrary Spaced Grids" (1988)
///////////////////////////////////////////////////////////////////////////////
template <typename TCoeffs>
class FDCoefficients_Base
{
public:
    FDCoefficients_Base(const unsigned int N) :
        m_N(N), m_coeffs(N)
    {}

    virtual ~FDCoefficients_Base() {}

    virtual void setup(const double* const grid_spacing, const unsigned int ncells, const double xS=0.0) = 0;

    inline const TCoeffs& get_coefficients() const { return m_coeffs; }

protected:
    const unsigned int m_N;
    TCoeffs m_coeffs;

    // delta coefficients based on the Fornberg paper
    /* M: max. order of derivative < M (e.g. max. first order -> M = 2) */
    /* N: number of grid points */
    template <size_t M, size_t N>
    void _compute_delta(
            const double x0,      /* point of evaluation */
            const double xi[N],   /* array with N grid point coordinates */
            double delta[M][N][N] /* output coefficients */
            )
    {
        // initialize to zero
        for (size_t i = 0; i < M; ++i)
            for (size_t j = 0; j < N; ++j)
                for (size_t k = 0; k < N; ++k)
                    delta[i][j][k] = 0.0;
        delta[0][0][0] = 1.0;
        double c1 = 1.0;
        for (size_t n = 1; n < N; ++n)
        {
            double c2 = 1.0;
            size_t nu = 0;
            while (true)
            {
                const double c3 = xi[n] - xi[nu];
                c2 *= c3;
                for (size_t m = 0; m < std::min(n+1, M); ++m)
                {
                    double c4 = 0.0;
                    if (0 != m)
                        c4 = m * delta[m-1][n-1][nu];
                    delta[m][n][nu] = ( (xi[n] - x0) * delta[m][n-1][nu] - c4) / c3;
                }
                if (nu >= n-1)
                    break;
                nu += 1;
            }

            for (size_t m = 0; m < std::min(n+1, M); ++m)
            {
                double c4 = 0.0;
                if (0 != m)
                    c4 = m * delta[m-1][n-1][n-1];
                delta[m][n][n] = c1/c2 * (c4 - (xi[n-1] - x0) * delta[m][n-1][n-1]);
            }

            c1 = c2;
        }
    }
};


class FDCoefficients_2ndOrder: public FDCoefficients_Base<CoeffsFD_All_2nd_t>
{
public:
    FDCoefficients_2ndOrder(const unsigned int Nblocks) :
        FDCoefficients_Base<CoeffsFD_All_2nd_t>(Nblocks*(_BLOCKSIZE_+1)), // each block has _BLOCKSIZE_+1 coeffs
        m_Nblocks(Nblocks)
    {}

    virtual ~FDCoefficients_2ndOrder() {}

    static constexpr unsigned int HALO_S = 1;
    static constexpr unsigned int HALO_E = 1;

    virtual void setup(const double* const grid_spacing, const unsigned int ncells, const double xS)
    {
        assert(ncells == m_Nblocks*_BLOCKSIZE_);

        // get references to coefficient arrays
        Real* c_m1 = this->m_coeffs.get_coeffs(-1);
        Real* c_00 = this->m_coeffs.get_coeffs(0);
        Real* c_p1 = this->m_coeffs.get_coeffs(1);

        // compute coefficients for first derivative
        double x_eval = xS;
        for (int b = 0; b < m_Nblocks; ++b)
        {
            const double* const bspacing = grid_spacing + b*_BLOCKSIZE_;
            Real* const bc_m1 = c_m1 + b*(_BLOCKSIZE_+1);
            Real* const bc_00 = c_00 + b*(_BLOCKSIZE_+1);
            Real* const bc_p1 = c_p1 + b*(_BLOCKSIZE_+1);

            // compute coefficients at face x_{i-1/2}
            {
                const double hm1 = *(bspacing-1);
                const double h00 = *bspacing;
                const double hp1 = *(bspacing+1);
                const double xi[3] = {
                    x_eval - 0.5*hm1,
                    x_eval + 0.5*h00,
                    x_eval + h00 + 0.5*hp1
                };
                double delta[2][3][3];
                this->template _compute_delta<2,3>(x_eval, xi, delta);
                bc_m1[0] = delta[1][2][0];
                bc_00[0] = delta[1][2][1];
                bc_p1[0] = delta[1][2][2];
            }
            // compute coefficients at faces x_{i+1/2}
            for (int i = 0; i < (int)_BLOCKSIZE_; ++i)
            {
                const double* const si = bspacing + i;
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);

                x_eval += h00;
                const double xi[3] = {
                    x_eval - 0.5*hm1 - h00,
                    x_eval - 0.5*h00,
                    x_eval + 0.5*hp1
                };
                double delta[2][3][3];
                this->template _compute_delta<2,3>(x_eval, xi, delta);
                bc_m1[i+1] = delta[1][2][0];
                bc_00[i+1] = delta[1][2][1];
                bc_p1[i+1] = delta[1][2][2];
            }
        }
    }

private:
    const unsigned int m_Nblocks;
};


class FDCoefficients_4thOrder: public FDCoefficients_Base<CoeffsFD_All_4th_t>
{
public:
    FDCoefficients_4thOrder(const unsigned int Nblocks) :
        FDCoefficients_Base<CoeffsFD_All_4th_t>(Nblocks*(_BLOCKSIZE_+1)), // each block has _BLOCKSIZE_+1 coeffs
        m_Nblocks(Nblocks)
    {}

    virtual ~FDCoefficients_4thOrder() {}

    static const unsigned int HALO_S = 2;
    static const unsigned int HALO_E = 2;

    virtual void setup(const double* const grid_spacing, const unsigned int ncells, const double xS)
    {
        assert(ncells == m_Nblocks*_BLOCKSIZE_);

        // get references to coefficient arrays
        Real* c_m2 = this->m_coeffs.get_coeffs(-2);
        Real* c_m1 = this->m_coeffs.get_coeffs(-1);
        Real* c_00 = this->m_coeffs.get_coeffs(0);
        Real* c_p1 = this->m_coeffs.get_coeffs(1);
        Real* c_p2 = this->m_coeffs.get_coeffs(2);

        // compute coefficients for first derivative
        double x_eval = xS;
        for (int b = 0; b < m_Nblocks; ++b)
        {
            const double* const bspacing = grid_spacing + b*_BLOCKSIZE_;
            Real* const bc_m2 = c_m2 + b*(_BLOCKSIZE_+1);
            Real* const bc_m1 = c_m1 + b*(_BLOCKSIZE_+1);
            Real* const bc_00 = c_00 + b*(_BLOCKSIZE_+1);
            Real* const bc_p1 = c_p1 + b*(_BLOCKSIZE_+1);
            Real* const bc_p2 = c_p2 + b*(_BLOCKSIZE_+1);

            // compute coefficients at face x_{i-1/2}
            {
                const double hm2 = *(bspacing-2);
                const double hm1 = *(bspacing-1);
                const double h00 = *bspacing;
                const double hp1 = *(bspacing+1);
                const double hp2 = *(bspacing+2);
                const double xi[5] = {
                    x_eval - 0.5*hm2 -     hm1,
                    x_eval           - 0.5*hm1,
                    x_eval                     + 0.5*h00,
                    x_eval                     +     h00 + 0.5*hp1,
                    x_eval                     +     h00     + hp1 + 0.5*hp2
                };
                double delta[2][5][5];
                this->template _compute_delta<2,5>(x_eval, xi, delta);
                bc_m2[0] = delta[1][4][0];
                bc_m1[0] = delta[1][4][1];
                bc_00[0] = delta[1][4][2];
                bc_p1[0] = delta[1][4][3];
                bc_p2[0] = delta[1][4][4];
            }
            // compute coefficients at faces x_{i+1/2}
            for (int i = 0; i < (int)_BLOCKSIZE_; ++i)
            {
                const double* const si = bspacing + i;
                const double hm2 = *(si-2);
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);
                const double hp2 = *(si+2);

                x_eval += h00;
                const double xi[5] = {
                    x_eval - 0.5*hm2 -     hm1 -     h00,
                    x_eval           - 0.5*hm1 -     h00,
                    x_eval                     - 0.5*h00,
                    x_eval                               + 0.5*hp1,
                    x_eval                               +     hp1 + 0.5*hp2
                };
                double delta[2][5][5];
                this->template _compute_delta<2,5>(x_eval, xi, delta);
                bc_m2[i+1] = delta[1][4][0];
                bc_m1[i+1] = delta[1][4][1];
                bc_00[i+1] = delta[1][4][2];
                bc_p1[i+1] = delta[1][4][3];
                bc_p2[i+1] = delta[1][4][4];
            }
        }
    }

private:
    const unsigned int m_Nblocks;
};

#ifdef _SECONDORDER_FD_CORE_
typedef CoeffsFD_Block_2nd_t    CoeffsFD_Block_t; // 2nd-order
typedef CoeffsFD_All_2nd_t      CoeffsFD_All_t;
typedef FDCoefficients_2ndOrder FDCoefficientGenerator;
#else
typedef CoeffsFD_Block_4th_t CoeffsFD_Block_t; // 4th-order
typedef CoeffsFD_All_4th_t   CoeffsFD_All_t;
typedef FDCoefficients_4thOrder FDCoefficientGenerator;
#endif /* _SECONDORDER_FD_CORE_ */

#endif /* FDCOEFFICIENTS_H_DLHT8AXM */
