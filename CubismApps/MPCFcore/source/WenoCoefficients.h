/*
 *  WenoCoefficients.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger 05/07/2017
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef WENOCOEFFICIENTS_H_QX69ZBST
#define WENOCOEFFICIENTS_H_QX69ZBST

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <Cubism/MeshMap.h>
#include "common.h"
using namespace cubism;

///////////////////////////////////////////////////////////////////////////////
// coefficient types
///////////////////////////////////////////////////////////////////////////////
// per block
struct CoeffsWENO_Block_t
{
    // polynomials
    Real c[9][_BLOCKSIZE_+1];

    // ideal weights
    Real d[3][_BLOCKSIZE_+1];

    // smoothness indicators
    Real b[9][_BLOCKSIZE_+1];
};


// low level
struct Polynomial_coeffs_All_t
{
    Real* c00;
    Real* c01;
    Real* c02;

    Real* c10;
    Real* c11;
    Real* c12;

    Real* c20;
    Real* c21;
    Real* c22;

    Polynomial_coeffs_All_t(const unsigned int N)
    {
        posix_memalign((void**)&c00, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&c01, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&c02, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&c10, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&c11, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&c12, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&c20, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&c21, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&c22, _ALIGNBYTES_, N*sizeof(Real));
    }

    ~Polynomial_coeffs_All_t()
    {
        free(c00);
        free(c01);
        free(c02);
        free(c10);
        free(c11);
        free(c12);
        free(c20);
        free(c21);
        free(c22);
    }
};

struct Ideal_weights_All_t
{
    Real* d0;
    Real* d1;
    Real* d2;

    Ideal_weights_All_t(const unsigned int N)
    {
        posix_memalign((void**)&d0, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&d1, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&d2, _ALIGNBYTES_, N*sizeof(Real));
    }

    ~Ideal_weights_All_t()
    {
        free(d0);
        free(d1);
        free(d2);
    }
};

struct Smoothness_coeffs_All_t
{
    Real* b00;
    Real* b01;
    Real* b02;

    Real* b10;
    Real* b11;
    Real* b12;

    Real* b20;
    Real* b21;
    Real* b22;

    Smoothness_coeffs_All_t() :
        b00(nullptr), b01(nullptr), b02(nullptr),
        b10(nullptr), b11(nullptr), b12(nullptr),
        b20(nullptr), b21(nullptr), b22(nullptr)
    {}

    void alloc(const unsigned int N)
    {
        posix_memalign((void**)&b00, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&b01, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&b02, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&b10, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&b11, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&b12, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&b20, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&b21, _ALIGNBYTES_, N*sizeof(Real));
        posix_memalign((void**)&b22, _ALIGNBYTES_, N*sizeof(Real));
    }

    void dealloc()
    {
        free(b00);
        free(b01);
        free(b02);
        free(b10);
        free(b11);
        free(b12);
        free(b20);
        free(b21);
        free(b22);
    }
};

// complete set of WENO coefficients for the whole domain extent in one spatial
// direction (only used for the initial coefficient computation, they will then
// be distributed to the blocks and will be destroyed afterwards)
struct CoeffsWENO_All_t
{
    Smoothness_coeffs_All_t b;
    Polynomial_coeffs_All_t c;
    Ideal_weights_All_t d;

    CoeffsWENO_All_t(const unsigned int N) :
        b(), c(N), d(N)
    {}
};
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// General WENO3/WENO5 coefficient classes
///////////////////////////////////////////////////////////////////////////////
class WenoCoefficients_Base
{
public:
    WenoCoefficients_Base(const unsigned int N) :
        m_N(N), m_minus(N), m_plus(N)
    {
        // allocate one extra element because smoothness indicators for
        // minus/plus coefficients are shifted by 1 (see implementations below)
        m_smooth_indicator.alloc(N+1);
    }

    virtual ~WenoCoefficients_Base()
    {
        m_smooth_indicator.dealloc();
    }

    virtual void setup(const double* const grid_spacing, const unsigned int ncells) = 0;

    inline const CoeffsWENO_All_t& get_coefficients_minus() const { return m_minus; }
    inline const CoeffsWENO_All_t& get_coefficients_plus() const { return m_plus; }

protected:
    const unsigned int      m_N;
    CoeffsWENO_All_t        m_minus;
    CoeffsWENO_All_t        m_plus;
    Smoothness_coeffs_All_t m_smooth_indicator;
};


///////////////////////////////////////////////////////////////////////////////
// Weno5 coefficients based on Coralic & Colonius, JCP 274 (2014) 95-121
class Weno5Coefficients_Coralic : public WenoCoefficients_Base
{
public:
    Weno5Coefficients_Coralic(const unsigned int N) : WenoCoefficients_Base(N) {}
    virtual ~Weno5Coefficients_Coralic() {}

    static const unsigned int HALO_S = 3;
    static const unsigned int HALO_E = 3;

    virtual void setup(const double* const grid_spacing, const unsigned int ncells)
    {
        // compute coefficients
        assert(ncells+1 == m_N);

        // minus: coefficients correspond to faces located at x_{i - 1/2}
        // plus:  coefficients correspond to faces located at x_{i + 1/2}
        //
        // NOTE: The WENO kernels in the solver are implemented using the
        // notation  -|+, where | corresponds to a face located at x_{i + 1/2}.
        // Therefore, when the reconstruction is carried out for -|, this
        // corresponds to the PLUS coefficients in cell i, whereas |+
        // corresponds to the MINUS coefficients in cell i+1.  The coefficients
        // computed below are stored for all faces n_N = ncells+1 along some
        // direction with respect to the cell i, that is, +|-.

        for (int i = 0; i < (int)m_N; ++i)
        {
            // minus
            ///////////////////////////////////////////////////////////////////
            {
                const double* const si = grid_spacing + i;
                const double hm2 = *(si-2);
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);
                const double hp2 = *(si+2);

                m_minus.c.c00[i] = 1.0;
                m_minus.c.c01[i] = -(h00*(h00+hp1 + h00+hp1+hp2)) / ((h00+hp1)*(h00+hp1+hp2));
                m_minus.c.c02[i] =  (h00*(h00+hp1)) / ((h00+hp1+hp2)*(hp1+hp2));
                m_minus.c.c10[i] = 1.0;
                m_minus.c.c11[i] = -(h00*(h00+hp1)) / ((hm1+h00)*(hm1+h00+hp1));
                m_minus.c.c12[i] = -(hm1*h00) / ((hm1+h00+hp1)*(h00+hp1));
                m_minus.c.c20[i] = 1.0;
                m_minus.c.c21[i] =  (hm1*h00) / ((hm2+hm1)*(hm2+hm1+h00));
                m_minus.c.c22[i] = -(h00*(hm2+hm1 + hm1+h00)) / ((hm2+hm1+h00)*(hm1+h00));

                m_minus.d.d0[i] = ((hm2+hm1)*hm1) / ((hm2+hm1+h00+hp1+hp2)*(hm1+h00+hp1+hp2));
                m_minus.d.d2[i] = ((h00+hp1)*(h00+hp1+hp2)) / ((hm2+hm1+h00+hp1)*(hm2+hm1+h00+hp1+hp2));
                m_minus.d.d1[i] = 1.0 - m_minus.d.d0[i] - m_minus.d.d2[i];
            }

            // plus
            ///////////////////////////////////////////////////////////////////
            {
                const double* const si = grid_spacing + (i-1);
                const double hm2 = *(si-2);
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);
                const double hp2 = *(si+2);

                m_plus.c.c00[i] =  1.0;
                m_plus.c.c01[i] =  (h00*(h00+hp1 + hp1+hp2)) / ((h00+hp1)*(h00+hp1+hp2));
                m_plus.c.c02[i] = -(h00*hp1) / ((h00+hp1+hp2)*(hp1+hp2));
                m_plus.c.c10[i] =  1.0;
                m_plus.c.c11[i] =  (h00*hp1) / ((hm1+h00)*(hm1+h00+hp1));
                m_plus.c.c12[i] =  ((hm1+h00)*h00) / ((hm1+h00+hp1)*(h00+hp1));
                m_plus.c.c20[i] =  1.0;
                m_plus.c.c21[i] = -((hm1+h00)*h00) / ((hm2+hm1)*(hm2+hm1+h00));
                m_plus.c.c22[i] =  (h00*(hm2+hm1+h00 + hm1+h00)) / ((hm2+hm1+h00)*(hm1+h00));

                m_plus.d.d0[i] = ((hm2+hm1+h00)*(hm1+h00)) / ((hm2+hm1+h00+hp1+hp2)*(hm1+h00+hp1+hp2));
                m_plus.d.d2[i] = (hp1*(hp1+hp2)) / ((hm2+hm1+h00+hp1)*(hm2+hm1+h00+hp1+hp2));
                m_plus.d.d1[i] = 1.0 - m_plus.d.d0[i] - m_plus.d.d2[i];
            }
        }

        // smoothness indicator coefficients
        for (int i = 0; i < (int)(m_N+1); ++i)
        {
            const double* const si = grid_spacing + (i-1);
            const double hm2 = *(si-2);
            const double hm1 = *(si-1);
            const double h00 = *si;
            const double hp1 = *(si+1);
            const double hp2 = *(si+2);

            const double fac = 4.0*h00*h00;
            {
                const double c0 = h00+hp1 + hp1+hp2;
                const double c1 = (h00+hp1)*(h00+hp1+hp2);
                m_smooth_indicator.b00[i] =  fac*( (10.0*h00*h00 + h00*(h00+hp1 + hp1+hp2) + c0*c0) / (c1*c1) );
            }
            {
                const double c0 = h00+hp1+hp2;
                m_smooth_indicator.b01[i] = -fac*( (19.0*h00*h00 - h00*(hp1+hp2) + 2.0*(h00+hp1)*(h00+hp1 + hp1+hp2)) / ((h00+hp1)*(hp1+hp2)*c0*c0) );
            }
            {
                const double c0 = (h00+hp1+hp2)*(hp1+hp2);
                m_smooth_indicator.b02[i] =  fac*( (10.0*h00*h00 + h00*hp1 + hp1*hp1) / (c0*c0) );
            }

            {
                const double c0 = hm1+h00+hp1;
                m_smooth_indicator.b10[i] = -fac*( (h00*(hm1 + 20.0*h00) - (h00+hp1)*(2.0*hm1 + h00)) / ((hm1+h00)*(h00+hp1)*c0*c0) );
            }
            {
                const double c0 = (hm1+h00)*(hm1+h00+hp1);
                m_smooth_indicator.b11[i] =  fac*( (10.0*h00*h00 + h00*hp1 + hp1*hp1) / (c0*c0) );
            }
            {
                const double c0 = (hm1+h00+hp1)*(h00+hp1);
                m_smooth_indicator.b12[i] =  fac*( (10.0*h00*h00 + hm1*h00 + hm1*hm1) / (c0*c0) );
            }

            {
                const double c0 = hm2+hm1 + hm1;
                const double c1 = (hm2+hm1+h00)*(hm1+h00);
                m_smooth_indicator.b20[i] =  fac*( (12.0*h00*h00 + 3.0*h00*(hm2+hm1 + hm1) + c0*c0) / (c1*c1) );
            }
            {
                const double c0 = hm2+hm1+h00;
                m_smooth_indicator.b21[i] = -fac*( (19.0*h00*h00 - (hm2+hm1)*h00 + 2.0*(hm1+h00)*(hm2+hm1 + hm1+h00)) / ((hm2+hm1)*(hm1+h00)*c0*c0) );
            }
            {
                const double c0 = (hm2+hm1)*(hm2+hm1+h00);
                m_smooth_indicator.b22[i] =  fac*( (10.0*h00*h00 + hm1*h00 + hm1*hm1) / (c0*c0) );
            }
        }

        m_minus.b.b00 = &m_smooth_indicator.b00[1];
        m_minus.b.b01 = &m_smooth_indicator.b01[1];
        m_minus.b.b02 = &m_smooth_indicator.b02[1];
        m_minus.b.b10 = &m_smooth_indicator.b10[1];
        m_minus.b.b11 = &m_smooth_indicator.b11[1];
        m_minus.b.b12 = &m_smooth_indicator.b12[1];
        m_minus.b.b20 = &m_smooth_indicator.b20[1];
        m_minus.b.b21 = &m_smooth_indicator.b21[1];
        m_minus.b.b22 = &m_smooth_indicator.b22[1];

        m_plus.b.b00 = &m_smooth_indicator.b00[0];
        m_plus.b.b01 = &m_smooth_indicator.b01[0];
        m_plus.b.b02 = &m_smooth_indicator.b02[0];
        m_plus.b.b10 = &m_smooth_indicator.b10[0];
        m_plus.b.b11 = &m_smooth_indicator.b11[0];
        m_plus.b.b12 = &m_smooth_indicator.b12[0];
        m_plus.b.b20 = &m_smooth_indicator.b20[0];
        m_plus.b.b21 = &m_smooth_indicator.b21[0];
        m_plus.b.b22 = &m_smooth_indicator.b22[0];
    }
};
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Weno3 coefficients: Derived by hand (notebook) Wed Jun 21 2017
class Weno3Coefficients : public WenoCoefficients_Base
{
public:
    Weno3Coefficients(const unsigned int N) : WenoCoefficients_Base(N) {}
    virtual ~Weno3Coefficients() {}

    static const unsigned int HALO_S = 2;
    static const unsigned int HALO_E = 2;

    virtual void setup(const double* const grid_spacing, const unsigned int ncells)
    {
        // compute coefficients
        assert(ncells+1 == m_N);

        // minus: coefficients correspond to faces located at x_{i - 1/2}
        // plus:  coefficients correspond to faces located at x_{i + 1/2}
        //
        // NOTE: The WENO kernels in the solver are implemented using the
        // notation  -|+, where | corresponds to a face located at x_{i + 1/2}.
        // Therefore, when the reconstruction is carried out for -|, this
        // corresponds to the PLUS coefficients in cell i, whereas |+
        // corresponds to the MINUS coefficients in cell i+1.  The coefficients
        // computed below are stored for all faces n_N = ncells+1 along some
        // direction with respect to the cell i, that is, +|-.

        for (int i = 0; i < (int)m_N; ++i)
        {
            // minus
            ///////////////////////////////////////////////////////////////////
            {
                const double* const si = grid_spacing + i;
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);

                m_minus.c.c00[i] = (h00+hp1)/hp1 - h00*h00/((h00+hp1)*hp1);
                m_minus.c.c01[i] = -h00/(h00+hp1);
                m_minus.c.c10[i] = 1.0 - hm1/h00 + hm1*hm1/((hm1+h00)*h00);
                m_minus.c.c11[i] = hm1/(hm1+h00);

                m_minus.d.d0[i] = hm1/(hm1+h00+hp1);
                m_minus.d.d1[i] = 1.0 - m_minus.d.d0[i];
            }

            // plus
            ///////////////////////////////////////////////////////////////////
            {
                const double* const si = grid_spacing + (i-1);
                const double hm1 = *(si-1);
                const double h00 = *si;
                const double hp1 = *(si+1);

                m_plus.c.c00[i] = 1.0 - h00/hp1 + h00*h00/((h00+hp1)*hp1);
                m_plus.c.c01[i] = h00/(h00+hp1);
                m_plus.c.c10[i] = hm1/(hm1+h00) - 1.0;
                m_plus.c.c11[i] = h00/(hm1+h00) + 1.0;

                m_plus.d.d0[i] = (hm1+h00)/(hm1+h00+hp1);
                m_plus.d.d1[i] = 1.0 - m_plus.d.d0[i];
            }

            m_minus.c.c02[i] = 0.0;
            m_minus.c.c12[i] = 0.0;
            m_minus.c.c20[i] = 0.0;
            m_minus.c.c21[i] = 0.0;
            m_minus.c.c22[i] = 0.0;
            m_minus.d.d2[i] = 0.0;

            m_plus.c.c02[i] = 0.0;
            m_plus.c.c12[i] = 0.0;
            m_plus.c.c20[i] = 0.0;
            m_plus.c.c21[i] = 0.0;
            m_plus.c.c22[i] = 0.0;
            m_plus.d.d2[i] = 0.0;
        }

        // smoothness indicator coefficients
        for (int i = 0; i < (int)(m_N+1); ++i)
        {
            const double* const si = grid_spacing + (i-1);
            const double hm1 = *(si-1);
            const double h00 = *si;
            const double hp1 = *(si+1);

            const double fac = 4.0*h00*h00;
            {
                const double c0 = (h00/(h00+hp1) - 1.0)/hp1;
                const double c1 = 1.0/(h00+hp1);
                m_smooth_indicator.b00[i] = fac*c0*c0;
                m_smooth_indicator.b01[i] = fac*2.0*c0*c1;
                m_smooth_indicator.b02[i] = fac*c1*c1;
            }

            {
                const double c0 = (hm1/(hm1+h00) - 1.0)/h00;
                const double c1 = 1.0/(hm1+h00);
                m_smooth_indicator.b10[i] = fac*c0*c0;
                m_smooth_indicator.b11[i] = fac*2.0*c0*c1;
                m_smooth_indicator.b12[i] = fac*c1*c1;
            }

            m_smooth_indicator.b20[i] = 0.0;
            m_smooth_indicator.b21[i] = 0.0;
            m_smooth_indicator.b22[i] = 0.0;
        }

        m_minus.b.b00 = &m_smooth_indicator.b00[1];
        m_minus.b.b01 = &m_smooth_indicator.b01[1];
        m_minus.b.b02 = &m_smooth_indicator.b02[1];
        m_minus.b.b10 = &m_smooth_indicator.b10[1];
        m_minus.b.b11 = &m_smooth_indicator.b11[1];
        m_minus.b.b12 = &m_smooth_indicator.b12[1];
        m_minus.b.b20 = &m_smooth_indicator.b20[1];
        m_minus.b.b21 = &m_smooth_indicator.b21[1];
        m_minus.b.b22 = &m_smooth_indicator.b22[1];

        m_plus.b.b00 = &m_smooth_indicator.b00[0];
        m_plus.b.b01 = &m_smooth_indicator.b01[0];
        m_plus.b.b02 = &m_smooth_indicator.b02[0];
        m_plus.b.b10 = &m_smooth_indicator.b10[0];
        m_plus.b.b11 = &m_smooth_indicator.b11[0];
        m_plus.b.b12 = &m_smooth_indicator.b12[0];
        m_plus.b.b20 = &m_smooth_indicator.b20[0];
        m_plus.b.b21 = &m_smooth_indicator.b21[0];
        m_plus.b.b22 = &m_smooth_indicator.b22[0];
    }
};
///////////////////////////////////////////////////////////////////////////////

#ifdef _WENO3_
    typedef Weno3Coefficients WENOCoefficientGenerator;
#else
    typedef Weno5Coefficients_Coralic WENOCoefficientGenerator;
#endif /* _WENO3_ */

#endif /* WENOCOEFFICIENTS_H_QX69ZBST */
