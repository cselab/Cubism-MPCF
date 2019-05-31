//
//  Gradient_CPP_nonuniform.h
//  MPCFcore
//
//  Creadet by Fabian Wermelinger on 8/17/2017
//  Copyright 2017 ETH Zurich. All rights reserved.
//
#ifndef GRADIENT_CPP_NONUNIFORM_H_QTLUZZ6F
#define GRADIENT_CPP_NONUNIFORM_H_QTLUZZ6F

#include <cassert>
#include "DivTensor_CPP_nonuniform.h"
#include "SurfaceTension_CPP_nonuniform.h" // for the typedefs

template <typename TInput, typename TInputRing>
class Gradient_CPP_nonuniform :
    public DivTensor_CPP_nonuniform<TInput>
{
public:
    Gradient_CPP_nonuniform(const Real dt, const Real a=1.0, const Real sigma=1.0) :
        DivTensor_CPP_nonuniform<TInput>(),
        m_dt(dt), m_a(a), m_sigma(sigma)
    {}

    void compute_nonuniform(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
            Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts,
            const CoeffsFD_Block_t& xcoeffs, const CoeffsFD_Block_t& ycoeffs, const CoeffsFD_Block_t& zcoeffs,
            const Real* const invh_x, const Real* const invh_y, const Real* const invh_z);

protected:
    struct AssumedType { Real phix, phiy, phiz, dummy1, dummy2, dummy3, dummy4, dummy5; };

    const Real m_dt;
    const Real m_a;
    const Real m_sigma;

    TInputRing ringphi; //slices for the primitive values
    RingTempSOA_ST ringGradax, ringGraday, ringGradaz; //slices for the corner-gradients

    virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);

    virtual void _tensor_xface(const TempSOA_ST& nx0, const TempSOA_ST& nx1);
    virtual void _tensor_yface(const TempSOA_ST& ny0, const TempSOA_ST& ny1);
    virtual void _tensor_zface(const TempSOA_ST& nz, TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz);

    virtual void _average_xy();
    virtual void _average_z(const TempPiZSOA_ST& tzz0, const TempPiZSOA_ST& tzz1);

    virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);

    inline void _input_next() { ringphi.next(); }
    inline void _grad_next() { ringGradax.next(); ringGraday.next(); ringGradaz.next(); }
};


template <typename TInput, typename TInputRing>
void Gradient_CPP_nonuniform<TInput,TInputRing>::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    TInput& phi = this->ringphi.ref();

#ifdef _SECONDORDER_FD_CORE_
    for(int sy=0; sy<_BLOCKSIZE_+2; sy++)
        for(int sx=0; sx<_BLOCKSIZE_+2; sx++)
        {
            const int dx = sx-1;
            const int dy = sy-1;
#else
    for(int sy=0; sy<_BLOCKSIZE_+4; sy++)
        for(int sx=0; sx<_BLOCKSIZE_+4; sx++)
        {
            const int dx = sx-2;
            const int dy = sy-2;
#endif /* _SECONDORDER_FD_CORE_ */

            AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));
            phi.ref(dx, dy) = pt.phix;
            assert(!std::isnan(phi.ref(dx, dy)));
        }
}


template <typename TInput, typename TInputRing>
void Gradient_CPP_nonuniform<TInput,TInputRing>::_tensor_xface(const TempSOA_ST& nx0, const TempSOA_ST& nx1)
{
    for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiXSOA_ST::NX; ix++)
        {
            const Real nx = (nx0(ix,iy)+nx0(ix,iy+1)+nx1(ix,iy)+nx1(ix,iy+1));
            this->txx.ref(ix, iy) = nx;
            this->txy.ref(ix, iy) = 0.0;
            this->txz.ref(ix, iy) = 0.0;
        }
}

template <typename TInput, typename TInputRing>
void Gradient_CPP_nonuniform<TInput,TInputRing>::_tensor_yface(const TempSOA_ST& ny0, const TempSOA_ST& ny1)
{
    for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiYSOA_ST::NX; ix++)
        {
            const Real ny = (ny0(ix,iy)+ny0(ix+1,iy)+ny1(ix,iy)+ny1(ix+1,iy));
            this->tyx.ref(ix, iy) = 0.0;
            this->tyy.ref(ix, iy) = ny;
            this->tyz.ref(ix, iy) = 0.0;
        }
}

template <typename TInput, typename TInputRing>
void Gradient_CPP_nonuniform<TInput,TInputRing>::_tensor_zface(const TempSOA_ST& nz0,
        TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz)
{
    for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiZSOA_ST::NX; ix++)
        {
            const Real nz = (nz0(ix,iy)+nz0(ix+1,iy)+nz0(ix,iy+1)+nz0(ix+1,iy+1));
            tzx.ref(ix, iy) = 0.0;
            tzy.ref(ix, iy) = 0.0;
            tzz.ref(ix, iy) = nz;
        }
}

template <typename TInput, typename TInputRing>
void Gradient_CPP_nonuniform<TInput,TInputRing>::_average_xy()
{
    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
        {
            this->rhsu.ref(ix, iy) = this->txx(ix+1, iy) + this->txx(ix, iy);
            this->rhsv.ref(ix, iy) = this->tyy(ix, iy+1) + this->tyy(ix, iy);
        }
}

template <typename TInput, typename TInputRing>
void Gradient_CPP_nonuniform<TInput,TInputRing>::_average_z(const TempPiZSOA_ST& tzz0, const TempPiZSOA_ST& tzz1)
{
    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
            this->rhsw.ref(ix, iy) = tzz1(ix, iy) + tzz0(ix, iy);
}

template <typename TInput, typename TInputRing>
void Gradient_CPP_nonuniform<TInput,TInputRing>::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
{
// TODO: [fabianw@mavt.ethz.ch; Thu Feb 22 2018 07:14:53 PM (+0100)]
// no longer correct, will be 1/8 (1/2 comes from averaging to cell center) after interpolation
    const Real factor = this->m_sigma*m_dt/32.0;

    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
        {
            assert(!std::isnan(this->rhsu(ix, iy)));
            assert(!std::isnan(this->rhsv(ix, iy)));
            assert(!std::isnan(this->rhsw(ix, iy)));

            AssumedType& rhs = *(AssumedType*)(gptfirst + gptfloats*(ix + iy*rowgpts));

            rhs.phix  = this->m_a*rhs.phix + factor*this->rhsu(ix, iy);
            rhs.phiy  = this->m_a*rhs.phiy + factor*this->rhsv(ix, iy);
            rhs.phiz  = this->m_a*rhs.phiz + factor*this->rhsw(ix, iy);
        }
}

template <typename TInput, typename TInputRing>
void Gradient_CPP_nonuniform<TInput,TInputRing>::compute_nonuniform(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
        Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts,
        const CoeffsFD_Block_t& xcoeffs, const CoeffsFD_Block_t& ycoeffs, const CoeffsFD_Block_t& zcoeffs,
        const Real* const invh_x, const Real* const invh_y, const Real* const invh_z)
{
#ifdef _SECONDORDER_FD_CORE_
    _convert(srcfirst+0*srcfloats*slicesrcs, srcfloats, rowsrcs);
    this->_input_next();
    _convert(srcfirst+1*srcfloats*slicesrcs, srcfloats, rowsrcs);
    this->_input_next();
    _convert(srcfirst+2*srcfloats*slicesrcs, srcfloats, rowsrcs);

    this->_corners(this->ringphi(-2), this->ringphi(-1), this->ringphi(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
#else
    for(int islice=0; islice<4; islice++)
    {
        _convert(srcfirst+islice*srcfloats*slicesrcs, srcfloats, rowsrcs);
        this->_input_next();
    }
    _convert(srcfirst + 4*srcfloats*slicesrcs, srcfloats, rowsrcs);

    this->_corners(this->ringphi(-4), this->ringphi(-3), this->ringphi(-2), this->ringphi(-1), this->ringphi(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

    this->_tensor_zface(this->ringGradaz(), this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

    for(int islice=0; islice<_BLOCKSIZE_; islice++)
    {
        this->_tensors_next();
        this->_grad_next();

#ifdef _SECONDORDER_FD_CORE_
        this->_corners(this->ringphi(-2), this->ringphi(-1), this->ringphi(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
#else
        this->_corners(this->ringphi(-4), this->ringphi(-3), this->ringphi(-2), this->ringphi(-1), this->ringphi(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

        this->_tensor_xface(this->ringGradax(0), this->ringGradax(1));
        this->_tensor_yface(this->ringGraday(0), this->ringGraday(1));
        this->_tensor_zface(this->ringGradaz(), this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

        _average_xy();
        _average_z(this->ringtzz(-1), this->ringtzz(0));

        _copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);

        if (islice < _BLOCKSIZE_-1)
        {
            this->_input_next();
#ifdef _SECONDORDER_FD_CORE_
            _convert(srcfirst + (islice+3)*srcfloats*slicesrcs, srcfloats, rowsrcs);
#else
            _convert(srcfirst + (islice+5)*srcfloats*slicesrcs, srcfloats, rowsrcs);
#endif /* _SECONDORDER_FD_CORE_ */
        }
    }
}

#ifdef _SECONDORDER_FD_CORE_
typedef Gradient_CPP_nonuniform<InputSOA_ST,RingInputSOA2nd_ST> Gradient_CPP_NU; // NonUniform type def
#else
typedef Gradient_CPP_nonuniform<InputSOA_ST2,RingInputSOA4th_ST> Gradient_CPP_NU; // NonUniform type def
#endif /* _SECONDORDER_FD_CORE_ */

#endif /* GRADIENT_CPP_NONUNIFORM_H_QTLUZZ6F */
