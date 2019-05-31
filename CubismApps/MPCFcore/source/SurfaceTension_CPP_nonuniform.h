//
//  SurfaceTension_CPP_nonuniform.h
//  MPCFcore
//
//  Creadet by Fabian Wermelinger on 7/27/2017
//  Copyright 2017 ETH Zurich. All rights reserved.
//
#ifndef SURFACETENSION_CPP_NONUNIFORM_H_WGYVUDMO
#define SURFACETENSION_CPP_NONUNIFORM_H_WGYVUDMO

#include <cassert>

#include "DivTensor_CPP.h"
#include "DivTensor_CPP_nonuniform.h"
#include "SurfaceTension_CPP.h"

// 2nd-order types
typedef RingSOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1, 3> RingInputSOA2nd_ST;

// 4th-order types
typedef SOA2D<-2, _BLOCKSIZE_+2, -2, _BLOCKSIZE_+2> InputSOA_ST2;
typedef RingSOA2D<-2, _BLOCKSIZE_+2, -2, _BLOCKSIZE_+2, 5> RingInputSOA4th_ST;

template <typename TInput, typename TInputRing>
class SurfaceTension_CPP_nonuniform :
    public virtual DivTensor_CPP<TInput>,
    public virtual SurfaceTension_CPP_uniform<TInput,TInputRing>,
    public DivTensor_CPP_nonuniform<TInput>
{
public:
    SurfaceTension_CPP_nonuniform(const Real dt, const Real sigma=1.0, const Real a=1.0) :
        DivTensor_CPP<TInput>(),
        SurfaceTension_CPP_uniform<TInput,TInputRing>(0.0, 1.0, sigma, a), // dtinvh = 0 explicitly
        DivTensor_CPP_nonuniform<TInput>(),
        m_dt(dt)
    {}

    void compute_nonuniform(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
            Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts,
            const CoeffsFD_Block_t& xcoeffs, const CoeffsFD_Block_t& ycoeffs, const CoeffsFD_Block_t& zcoeffs,
            const Real* const invh_x, const Real* const invh_y, const Real* const invh_z);

protected:
    const Real m_dt;

    virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
    virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);
};


template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_nonuniform<TInput,TInputRing>::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    TInput& u = this->ringu.ref();
    TInput& v = this->ringv.ref();
    TInput& w = this->ringw.ref();
    TInput& l = this->ringaux.ref();

    typedef typename SurfaceTension_CPP_uniform<TInput,TInputRing>::AssumedType AssumedType;

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

#ifdef _CONVERTCLIP_
            const Real a1r1 = std::max(static_cast<Real>(0.0),pt.a1r1);
            const Real a2r2 = std::max(static_cast<Real>(0.0),pt.a2r2);
            const Real alpha2 = std::max(static_cast<Real>(ALPHAEPS), std::min(pt.A2, static_cast<Real>(1.0-ALPHAEPS)));
#else
            const Real a1r1 = pt.a1r1;
            const Real a2r2 = pt.a2r2;
#ifdef _ALPHACLIP_
            const Real alpha2 = std::max(static_cast<Real>(ALPHAEPS), std::min(pt.A2, static_cast<Real>(1.0-ALPHAEPS)));
#else
            const Real alpha2 = pt.A2;
#endif
#endif
            const Real rInv = static_cast<Real>(1.0) / (a1r1 + a2r2);

            u.ref(dx, dy) = pt.ru*rInv;
            v.ref(dx, dy) = pt.rv*rInv;
            w.ref(dx, dy) = pt.rw*rInv;

            // get liquid void fraction (see Perigaud and Saurel JCP 2005)
            l.ref(dx, dy) = alpha2;

            assert(!std::isnan(u.ref(dx, dy)));
            assert(!std::isnan(v.ref(dx, dy)));
            assert(!std::isnan(w.ref(dx, dy)));
            assert(!std::isnan(l.ref(dx, dy)));
        }
}

template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_nonuniform<TInput,TInputRing>::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
{
// TODO: [fabianw@mavt.ethz.ch; Thu Feb 22 2018 07:05:33 PM (+0100)]
// no longer correct, will be 1/4 after interpolation
    const Real factor = this->m_sigma*m_dt/16.0;

    typedef typename SurfaceTension_CPP_uniform<TInput,TInputRing>::AssumedType AssumedType;

    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
        {
            assert(!std::isnan(this->rhsu(ix, iy)));
            assert(!std::isnan(this->rhsv(ix, iy)));
            assert(!std::isnan(this->rhsw(ix, iy)));
            assert(!std::isnan(this->rhse(ix, iy)));

            AssumedType& rhs = *(AssumedType*)(gptfirst + gptfloats*(ix + iy*rowgpts));

            rhs.ru  = this->m_a*rhs.ru + factor*this->rhsu(ix, iy);
            rhs.rv  = this->m_a*rhs.rv + factor*this->rhsv(ix, iy);
            rhs.rw  = this->m_a*rhs.rw + factor*this->rhsw(ix, iy);
// TODO: [fabianw@mavt.ethz.ch; Fri Feb 23 2018 11:53:21 AM (+0100)]
// 0.5 is only correct for uniform grid.  Needs interpolation in udot
            rhs.E   = this->m_a*rhs.E  + 0.5*factor*this->rhse(ix, iy);
        }
}

template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_nonuniform<TInput,TInputRing>::compute_nonuniform(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
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

    this->_corners(this->ringaux(-2), this->ringaux(-1), this->ringaux(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
#else
    for(int islice=0; islice<4; islice++)
    {
        _convert(srcfirst+islice*srcfloats*slicesrcs, srcfloats, rowsrcs);
        this->_input_next();
    }
    _convert(srcfirst + 4*srcfloats*slicesrcs, srcfloats, rowsrcs);

    this->_corners(this->ringaux(-4), this->ringaux(-3), this->ringaux(-2), this->ringaux(-1), this->ringaux(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

    this->_tensor_zface(this->ringGradax(), this->ringGraday(), this->ringGradaz(), this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

#ifdef _SECONDORDER_FD_CORE_
    this->_udot_tz(this->ringu(-2), this->ringv(-2), this->ringw(-2), this->ringu(-1), this->ringv(-1), this->ringw(-1), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());
#else
    this->_udot_tz(this->ringu(-3), this->ringv(-3), this->ringw(-3), this->ringu(-2), this->ringv(-2), this->ringw(-2), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

    for(int islice=0; islice<_BLOCKSIZE_; islice++)
    {
        this->_tensors_next();
        this->_grad_next();

#ifdef _SECONDORDER_FD_CORE_
        this->_corners(this->ringaux(-2), this->ringaux(-1), this->ringaux(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
#else
        this->_corners(this->ringaux(-4), this->ringaux(-3), this->ringaux(-2), this->ringaux(-1), this->ringaux(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

        this->_tensor_xface(this->ringGradax(0), this->ringGradax(1), this->ringGraday(0), this->ringGraday(1), this->ringGradaz(0), this->ringGradaz(1));
        this->_tensor_yface(this->ringGradax(0), this->ringGradax(1), this->ringGraday(0), this->ringGraday(1), this->ringGradaz(0), this->ringGradaz(1));
        this->_tensor_zface(this->ringGradax(), this->ringGraday(), this->ringGradaz(), this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

#ifdef _SECONDORDER_FD_CORE_
        this->_udot_tx(this->ringu(-1), this->ringv(-1), this->ringw(-1));
        this->_udot_ty(this->ringu(-1), this->ringv(-1), this->ringw(-1));
        this->_udot_tz(this->ringu(-1), this->ringv(-1), this->ringw(-1), this->ringu(0), this->ringv(0), this->ringw(0), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());
#else
        this->_udot_tx(this->ringu(-2), this->ringv(-2), this->ringw(-2));
        this->_udot_ty(this->ringu(-2), this->ringv(-2), this->ringw(-2));
        this->_udot_tz(this->ringu(-2), this->ringv(-2), this->ringw(-2), this->ringu(-1), this->ringv(-1), this->ringw(-1), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

        const Real ihz = invh_z[islice];
        this->_div_dxy(invh_x, invh_y);
        this->_div_dz(this->ringtzx(-1), this->ringtzy(-1), this->ringtzz(-1), this->ringutz(-1), this->ringtzx(0), this->ringtzy(0), this->ringtzz(0), this->ringutz(0), ihz);

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
typedef SurfaceTension_CPP_nonuniform<InputSOA_ST,RingInputSOA2nd_ST> SurfaceTension_CPP_NU; // NonUniform type def
#else
typedef SurfaceTension_CPP_nonuniform<InputSOA_ST2,RingInputSOA4th_ST> SurfaceTension_CPP_NU; // NonUniform type def
#endif /* _SECONDORDER_FD_CORE_ */

#endif /* SURFACETENSION_CPP_NONUNIFORM_H_WGYVUDMO */
