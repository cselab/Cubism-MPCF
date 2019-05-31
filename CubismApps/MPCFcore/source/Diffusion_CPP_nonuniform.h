/*
 *  Diffusion_CPP_nonuniform.h
 *  MPCFcore
 *
 *  Created by Fabian Wermelinger on 7/24/17.
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef DIFFUSION_CPP_NONUNIFORM_H_R7TDKV6A
#define DIFFUSION_CPP_NONUNIFORM_H_R7TDKV6A

#include <cmath>

#include "DivTensor_CPP.h"
#include "SurfaceTension_CPP_nonuniform.h"
#include "Diffusion_CPP.h"


template <typename TInput, typename TInputRing>
class Diffusion_CPP_nonuniform :
    public virtual DivTensor_CPP<TInput>,
    public virtual SurfaceTension_CPP_uniform<TInput,TInputRing>,
    public virtual Diffusion_CPP_uniform<TInput,TInputRing>,
    public virtual SurfaceTension_CPP_nonuniform<TInput,TInputRing> // phew!
{
public:
    Diffusion_CPP_nonuniform(const Real dt, const Real mu1, const Real mu2, const Real a=1.0, const Real sigma=1.0) :
        DivTensor_CPP<TInput>(),
        SurfaceTension_CPP_uniform<TInput,TInputRing>(0.0, 1.0, sigma, a), // set dtinvh=0 explicitely
        Diffusion_CPP_uniform<TInput,TInputRing>(mu1, mu2, 0.0, 1.0, a, sigma),
        SurfaceTension_CPP_nonuniform<TInput,TInputRing>(dt, sigma, a)
    {}

    void compute_nonuniform(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
            Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts,
            const CoeffsFD_Block_t& xcoeffs, const CoeffsFD_Block_t& ycoeffs, const CoeffsFD_Block_t& zcoeffs,
            const Real* const invh_x, const Real* const invh_y, const Real* const invh_z);

protected:
    virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
};


template <typename TInput, typename TInputRing>
void Diffusion_CPP_nonuniform<TInput,TInputRing>::compute_nonuniform(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
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

    this->_corners(this->ringu(-2), this->ringu(-1), this->ringu(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
    this->_corners(this->ringv(-2), this->ringv(-1), this->ringv(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradvx.ref(), this->ringGradvy.ref(), this->ringGradvz.ref());
    this->_corners(this->ringw(-2), this->ringw(-1), this->ringw(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradwx.ref(), this->ringGradwy.ref(), this->ringGradwz.ref());
#else
    for(int islice=0; islice<4; islice++)
    {
        _convert(srcfirst+islice*srcfloats*slicesrcs, srcfloats, rowsrcs);
        this->_input_next();
    }
    _convert(srcfirst + 4*srcfloats*slicesrcs, srcfloats, rowsrcs);

    this->_corners(this->ringu(-4), this->ringu(-3), this->ringu(-2), this->ringu(-1), this->ringu(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
    this->_corners(this->ringv(-4), this->ringv(-3), this->ringv(-2), this->ringv(-1), this->ringv(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradvx.ref(), this->ringGradvy.ref(), this->ringGradvz.ref());
    this->_corners(this->ringw(-4), this->ringw(-3), this->ringw(-2), this->ringw(-1), this->ringw(0),
            xcoeffs, ycoeffs, zcoeffs,
            0, this->ringGradwx.ref(), this->ringGradwy.ref(), this->ringGradwz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

    this->_tensor_zface(this->ringGradax(), this->ringGradaz(), this->ringGradvy(), this->ringGradvz(),
            this->ringGradwx(), this->ringGradwy(), this->ringGradwz(),
            this->ringaux(-1), this->ringaux(),
            this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

#ifdef _SECONDORDER_FD_CORE_
    this->_udot_tz(this->ringu(-2), this->ringv(-2), this->ringw(-2),  this->ringu(-1), this->ringv(-1), this->ringw(-1), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());
#else
    this->_udot_tz(this->ringu(-3), this->ringv(-3), this->ringw(-3),  this->ringu(-2), this->ringv(-2), this->ringw(-2), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

    for(int islice=0; islice<_BLOCKSIZE_; islice++)
    {
        this->_tensors_next();
        Diffusion_CPP_uniform<TInput,TInputRing>::_grad_next();

#ifdef _SECONDORDER_FD_CORE_
        this->_corners(this->ringu(-2), this->ringu(-1), this->ringu(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
        this->_corners(this->ringv(-2), this->ringv(-1), this->ringv(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradvx.ref(), this->ringGradvy.ref(), this->ringGradvz.ref());
        this->_corners(this->ringw(-2), this->ringw(-1), this->ringw(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradwx.ref(), this->ringGradwy.ref(), this->ringGradwz.ref());
#else
        this->_corners(this->ringu(-4), this->ringu(-3), this->ringu(-2), this->ringu(-1), this->ringu(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
        this->_corners(this->ringv(-4), this->ringv(-3), this->ringv(-2), this->ringv(-1), this->ringv(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradvx.ref(), this->ringGradvy.ref(), this->ringGradvz.ref());
        this->_corners(this->ringw(-4), this->ringw(-3), this->ringw(-2), this->ringw(-1), this->ringw(0),
                xcoeffs, ycoeffs, zcoeffs,
                islice+1, this->ringGradwx.ref(), this->ringGradwy.ref(), this->ringGradwz.ref());
#endif /* _SECONDORDER_FD_CORE_ */

        this->_tensor_xface(this->ringGradax(-1), this->ringGraday(-1), this->ringGradaz(-1),
                this->ringGradax(), this->ringGraday(), this->ringGradaz(),
                this->ringGradvx(-1), this->ringGradvy(-1), this->ringGradvx(), this->ringGradvy(),
                this->ringGradwx(-1), this->ringGradwz(-1), this->ringGradwx(), this->ringGradwz(),
                this->ringaux(-1));

        this->_tensor_yface(this->ringGradax(-1), this->ringGraday(-1), this->ringGradax(), this->ringGraday(),
                this->ringGradvx(-1), this->ringGradvy(-1), this->ringGradvz(-1),
                this->ringGradvx(), this->ringGradvy(), this->ringGradvz(),
                this->ringGradwy(-1), this->ringGradwz(-1), this->ringGradwy(), this->ringGradwz(),
                this->ringaux(-1));

        this->_tensor_zface(this->ringGradax(), this->ringGradaz(),
                this->ringGradvy(), this->ringGradvz(),
                this->ringGradwx(), this->ringGradwy(), this->ringGradwz(),
                this->ringaux(-1), this->ringaux(),
                this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

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

        SurfaceTension_CPP_nonuniform<TInput,TInputRing>::_copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);

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


template <typename TInput, typename TInputRing>
void Diffusion_CPP_nonuniform<TInput,TInputRing>::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    TInput& u = this->ringu.ref();
    TInput& v = this->ringv.ref();
    TInput& w = this->ringw.ref();
    TInput& m = this->ringaux.ref(); // inherited and used for viscosity here

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
            const Real alpha1 = static_cast<Real>(1.0) - alpha2;
            const Real rInv = static_cast<Real>(1.0) / (a1r1 + a2r2);

            // velocity
            u.ref(dx, dy) = pt.ru*rInv;
            v.ref(dx, dy) = pt.rv*rInv;
            w.ref(dx, dy) = pt.rw*rInv;

            // viscosity
            m.ref(dx, dy) = alpha1*this->m_mu1 + alpha2*this->m_mu2;

            assert(!std::isnan(u.ref(dx, dy)));
            assert(!std::isnan(v.ref(dx, dy)));
            assert(!std::isnan(w.ref(dx, dy)));
            assert(!std::isnan(m.ref(dx, dy)));
        }
}

#ifdef _SECONDORDER_FD_CORE_
typedef Diffusion_CPP_nonuniform<InputSOA_ST,RingInputSOA2nd_ST> Diffusion_CPP_NU; // NonUniform type def
#else
typedef Diffusion_CPP_nonuniform<InputSOA_ST2,RingInputSOA4th_ST> Diffusion_CPP_NU; // NonUniform type def
#endif /* _SECONDORDER_FD_CORE_ */

#endif /* DIFFUSION_CPP_NONUNIFORM_H_R7TDKV6A */
