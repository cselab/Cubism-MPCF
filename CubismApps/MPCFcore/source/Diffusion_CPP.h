/*
 *  Diffusion_CPP.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/27/12.
 *  Revision by Fabian Wermelinger on 7/26/2017
 *  Copyright 2012/2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef DIFFUSION_CPP_H_EMCEKDQP
#define DIFFUSION_CPP_H_EMCEKDQP

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "DivTensor_CPP.h"
#include "SurfaceTension_CPP.h"

#ifdef _TEST_1D_DIFFUSION_
#warning "WARNING: 1D DIFFUSION TEST IS ENABLED!"
#endif /* _TEST_1D_DIFFUSION_ */


template <typename TInput, typename TInputRing>
class Diffusion_CPP_uniform :
    public virtual DivTensor_CPP<TInput>,
    public virtual SurfaceTension_CPP_uniform<TInput,TInputRing>
{
public:
    Diffusion_CPP_uniform(const Real mu1, const Real mu2,
            const Real dtinvh, const Real h,
            const Real a=1.0, const Real sigma=1.0) :
        DivTensor_CPP<TInput>(),
        SurfaceTension_CPP_uniform<TInput,TInputRing>(dtinvh, h, sigma, a),
        m_mu1(mu1), m_mu2(mu2)
    {}

    void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
            Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts);

    void hpc_info(float& flop_convert, int& traffic_convert,
            float& flop_corners, int& traffic_corner,
            float& flop_tensorface, int& traffic_tensorface,
            float& flop_tdotu, int& traffic_tdotu,
            float& flop_div, int& traffic_div,
            float& flop_copyback, int& traffic_copyback,
            size_t& footprint);

protected:
    typedef Real RealTemp;

    const RealTemp m_mu1, m_mu2;

    // slices for velocity gradients, ringGradu = ringGrada inherited
    RingTempSOA_ST ringGradvx, ringGradvy, ringGradvz, ringGradwx, ringGradwy, ringGradwz;

    virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);

    virtual void _tensor_xface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
                               const TempSOA_ST& ny0, const TempSOA_ST& ny1,
                               const TempSOA_ST& nz0, const TempSOA_ST& nz1)
    { _oh_no(); }

    virtual void _tensor_yface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
                               const TempSOA_ST& ny0, const TempSOA_ST& ny1,
                               const TempSOA_ST& nz0, const TempSOA_ST& nz1)
    { _oh_no(); }

    virtual void _tensor_zface(const TempSOA_ST& nx, const TempSOA_ST& ny, const TempSOA_ST& nz,
                               TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz)
    { _oh_no(); }

    virtual void _tensor_xface(const TempSOA_ST& ux0, const TempSOA_ST& uy0, const TempSOA_ST& uz0,
            const TempSOA_ST& ux1, const TempSOA_ST& uy1, const TempSOA_ST& uz1,
            const TempSOA_ST& vx0, const TempSOA_ST& vy0, const TempSOA_ST& vx1, const TempSOA_ST& vy1,
            const TempSOA_ST& wx0, const TempSOA_ST& wz0, const TempSOA_ST& wx1, const TempSOA_ST& wz1,
            const TInput& mu0);

    virtual void _tensor_yface(const TempSOA_ST& ux0, const TempSOA_ST& uy0, const TempSOA_ST& ux1, const TempSOA_ST& uy1,
            const TempSOA_ST& vx0, const TempSOA_ST& vy0, const TempSOA_ST& vz0,
            const TempSOA_ST& vx1, const TempSOA_ST& vy1, const TempSOA_ST& vz1,
            const TempSOA_ST& wy0, const TempSOA_ST& wz0, const TempSOA_ST& wy1, const TempSOA_ST& wz1,
            const TInput& mu0);

    virtual void _tensor_zface(const TempSOA_ST& ux, const TempSOA_ST& uz,
            const TempSOA_ST& vy, const TempSOA_ST& vz,
            const TempSOA_ST& wx, const TempSOA_ST& wy, const TempSOA_ST& wz,
            const TInput& mu0, const TInput& mu1,
            TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz);

    inline void _grad_next()
    {
        // velocity gradients
        // Grada is reused from SurfaceTension_CPP for Gradu
        this->ringGradax.next(); this->ringGraday.next(); this->ringGradaz.next();
        ringGradvx.next(); ringGradvy.next(); ringGradvz.next();
        ringGradwx.next(); ringGradwy.next(); ringGradwz.next();
    }

private:
    inline void _oh_no()
    {
        fprintf(stderr, "You should not be calling this method here!\n");
        abort();
    }
};


template <typename TInput, typename TInputRing>
void Diffusion_CPP_uniform<TInput,TInputRing>::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    TInput& u = this->ringu.ref();
    TInput& v = this->ringv.ref();
    TInput& w = this->ringw.ref();
    TInput& m = this->ringaux.ref(); // inherited and used for viscosity here

    typedef typename SurfaceTension_CPP_uniform<TInput,TInputRing>::AssumedType AssumedType;

    for(int sy=0; sy<_BLOCKSIZE_+2; sy++)
        for(int sx=0; sx<_BLOCKSIZE_+2; sx++)
        {
            AssumedType pt = *(AssumedType*)(gptfirst + gptfloats*(sx + sy*rowgpts));

            const int dx = sx-1;
            const int dy = sy-1;

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
            m.ref(dx, dy) = alpha1*m_mu1 + alpha2*m_mu2;

            assert(!std::isnan(u.ref(dx, dy)));
            assert(!std::isnan(v.ref(dx, dy)));
            assert(!std::isnan(w.ref(dx, dy)));
            assert(!std::isnan(m.ref(dx, dy)));
        }
}


///////////////////////////////////////////////////////////////////////////////
// tensor implementation based on
// \tau_{ij} = \mu (\partial_i u_j + \partial_j u_i - 2/3\partial_k u_k \delta_{ij})
///////////////////////////////////////////////////////////////////////////////
template <typename TInput, typename TInputRing>
void Diffusion_CPP_uniform<TInput,TInputRing>::_tensor_xface(const TempSOA_ST& ux0, const TempSOA_ST& uy0, const TempSOA_ST& uz0,
        const TempSOA_ST& ux1, const TempSOA_ST& uy1, const TempSOA_ST& uz1,
        const TempSOA_ST& vx0, const TempSOA_ST& vy0, const TempSOA_ST& vx1, const TempSOA_ST& vy1,
        const TempSOA_ST& wx0, const TempSOA_ST& wz0, const TempSOA_ST& wx1, const TempSOA_ST& wz1,
        const TInput& mu0)
{
    const RealTemp factor = static_cast<RealTemp>(2.0/3.0);

    for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiXSOA_ST::NX; ix++)
        {
            const RealTemp dudx = (ux0(ix,iy)+ux0(ix,iy+1)+ux1(ix,iy)+ux1(ix,iy+1));
            const RealTemp dudy = (uy0(ix,iy)+uy0(ix,iy+1)+uy1(ix,iy)+uy1(ix,iy+1));
            const RealTemp dudz = (uz0(ix,iy)+uz0(ix,iy+1)+uz1(ix,iy)+uz1(ix,iy+1));
            const RealTemp dvdx = (vx0(ix,iy)+vx0(ix,iy+1)+vx1(ix,iy)+vx1(ix,iy+1));
            const RealTemp dvdy = (vy0(ix,iy)+vy0(ix,iy+1)+vy1(ix,iy)+vy1(ix,iy+1));
            const RealTemp dwdx = (wx0(ix,iy)+wx0(ix,iy+1)+wx1(ix,iy)+wx1(ix,iy+1));
            const RealTemp dwdz = (wz0(ix,iy)+wz0(ix,iy+1)+wz1(ix,iy)+wz1(ix,iy+1));

            const RealTemp mu = static_cast<RealTemp>(0.5) * (mu0(ix-1,iy) + mu0(ix,iy));

#ifdef _TEST_1D_DIFFUSION_
            this->txx.ref(ix, iy) = mu * dudx;
            this->txy.ref(ix, iy) = 0.0;
            this->txz.ref(ix, iy) = 0.0;
#else
            this->txx.ref(ix, iy) = mu * (static_cast<RealTemp>(2.0)*dudx - factor * (dudx + dvdy + dwdz));
            this->txy.ref(ix, iy) = mu * (dvdx + dudy);
            this->txz.ref(ix, iy) = mu * (dwdx + dudz);
#endif /* _TEST_1D_DIFFUSION_ */
        }
}

template <typename TInput, typename TInputRing>
void Diffusion_CPP_uniform<TInput,TInputRing>::_tensor_yface(const TempSOA_ST& ux0, const TempSOA_ST& uy0, const TempSOA_ST& ux1, const TempSOA_ST& uy1,
        const TempSOA_ST& vx0, const TempSOA_ST& vy0, const TempSOA_ST& vz0,
        const TempSOA_ST& vx1, const TempSOA_ST& vy1, const TempSOA_ST& vz1,
        const TempSOA_ST& wy0, const TempSOA_ST& wz0, const TempSOA_ST& wy1, const TempSOA_ST& wz1,
        const TInput& mu0)
{
    const RealTemp factor = static_cast<RealTemp>(2.0/3.0);

    for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiYSOA_ST::NX; ix++)
        {
            const RealTemp dudx = (ux0(ix,iy)+ux0(ix+1,iy)+ux1(ix,iy)+ux1(ix+1,iy));
            const RealTemp dudy = (uy0(ix,iy)+uy0(ix+1,iy)+uy1(ix,iy)+uy1(ix+1,iy));
            const RealTemp dvdx = (vx0(ix,iy)+vx0(ix+1,iy)+vx1(ix,iy)+vx1(ix+1,iy));
            const RealTemp dvdy = (vy0(ix,iy)+vy0(ix+1,iy)+vy1(ix,iy)+vy1(ix+1,iy));
            const RealTemp dvdz = (vz0(ix,iy)+vz0(ix+1,iy)+vz1(ix,iy)+vz1(ix+1,iy));
            const RealTemp dwdy = (wy0(ix,iy)+wy0(ix+1,iy)+wy1(ix,iy)+wy1(ix+1,iy));
            const RealTemp dwdz = (wz0(ix,iy)+wz0(ix+1,iy)+wz1(ix,iy)+wz1(ix+1,iy));

            const RealTemp mu = static_cast<RealTemp>(0.5) * (mu0(ix,iy-1) + mu0(ix,iy));

#ifdef _TEST_1D_DIFFUSION_
            this->tyx.ref(ix, iy) = 0.0;
            this->tyy.ref(ix, iy) = mu * dvdy;
            this->tyz.ref(ix, iy) = 0.0;
#else
            this->tyx.ref(ix, iy) = mu * (dudy + dvdx);
            this->tyy.ref(ix, iy) = mu * (static_cast<RealTemp>(2.0)*dvdy - factor*(dudx + dvdy + dwdz));
            this->tyz.ref(ix, iy) = mu * (dwdy + dvdz);
#endif /* _TEST_1D_DIFFUSION_ */
        }
}

template <typename TInput, typename TInputRing>
void Diffusion_CPP_uniform<TInput,TInputRing>::_tensor_zface(const TempSOA_ST& ux, const TempSOA_ST& uz,
        const TempSOA_ST& vy, const TempSOA_ST& vz,
        const TempSOA_ST& wx, const TempSOA_ST& wy, const TempSOA_ST& wz,
        const TInput& mu0, const TInput& mu1,
        TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz)
{
    const RealTemp factor = static_cast<RealTemp>(2.0/3.0);

    for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiZSOA_ST::NX; ix++)
        {
            const RealTemp dudx = (ux(ix,iy)+ux(ix+1,iy)+ux(ix,iy+1)+ux(ix+1,iy+1));
            const RealTemp dudz = (uz(ix,iy)+uz(ix+1,iy)+uz(ix,iy+1)+uz(ix+1,iy+1));
            const RealTemp dvdy = (vy(ix,iy)+vy(ix+1,iy)+vy(ix,iy+1)+vy(ix+1,iy+1));
            const RealTemp dvdz = (vz(ix,iy)+vz(ix+1,iy)+vz(ix,iy+1)+vz(ix+1,iy+1));
            const RealTemp dwdx = (wx(ix,iy)+wx(ix+1,iy)+wx(ix,iy+1)+wx(ix+1,iy+1));
            const RealTemp dwdy = (wy(ix,iy)+wy(ix+1,iy)+wy(ix,iy+1)+wy(ix+1,iy+1));
            const RealTemp dwdz = (wz(ix,iy)+wz(ix+1,iy)+wz(ix,iy+1)+wz(ix+1,iy+1));

            const RealTemp mu = static_cast<RealTemp>(0.5) * (mu0(ix,iy) + mu1(ix,iy));

#ifdef _TEST_1D_DIFFUSION_
            tzx.ref(ix, iy) = 0.0;
            tzy.ref(ix, iy) = 0.0;
            tzz.ref(ix, iy) = mu * dwdz;
#else
            tzx.ref(ix, iy) = mu * (dudz + dwdx);
            tzy.ref(ix, iy) = mu * (dvdz + dwdy);
            tzz.ref(ix, iy) = mu * (static_cast<RealTemp>(2.0)*dwdz - factor*(dudx + dvdy + dwdz));
#endif /* _TEST_1D_DIFFUSION_ */
        }
}

template <typename TInput, typename TInputRing>
void Diffusion_CPP_uniform<TInput,TInputRing>::compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
        Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
{
    _convert(srcfirst, srcfloats, rowsrcs);
    this->_input_next();

    _convert(srcfirst + srcfloats*slicesrcs, srcfloats, rowsrcs);

    this->_corners(this->ringu(-1), this->ringu(), this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
    this->_corners(this->ringv(-1), this->ringv(), ringGradvx.ref(), ringGradvy.ref(), ringGradvz.ref());
    this->_corners(this->ringw(-1), this->ringw(), ringGradwx.ref(), ringGradwy.ref(), ringGradwz.ref());

    _tensor_zface(this->ringGradax(), this->ringGradaz(), ringGradvy(), ringGradvz(),
            ringGradwx(), ringGradwy(), ringGradwz(),
            this->ringaux(-1), this->ringaux(),
            this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

    this->_udot_tz(this->ringu(-1), this->ringv(-1), this->ringw(-1),  this->ringu(), this->ringv(), this->ringw(), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());

    for(int islice=0; islice<_BLOCKSIZE_; islice++)
    {
        this->_tensors_next();
        this->_input_next();
        _grad_next();

        _convert(srcfirst + (islice+2)*srcfloats*slicesrcs, srcfloats, rowsrcs);

        this->_corners(this->ringu(-1), this->ringu(0), this->ringGradax.ref(), this->ringGraday.ref(), this->ringGradaz.ref());
        this->_corners(this->ringv(-1), this->ringv(0), ringGradvx.ref(), ringGradvy.ref(), ringGradvz.ref());
        this->_corners(this->ringw(-1), this->ringw(0), ringGradwx.ref(), ringGradwy.ref(), ringGradwz.ref());

        _tensor_xface(this->ringGradax(-1), this->ringGraday(-1), this->ringGradaz(-1),
                this->ringGradax(), this->ringGraday(), this->ringGradaz(),
                ringGradvx(-1), ringGradvy(-1), ringGradvx(), ringGradvy(),
                ringGradwx(-1), ringGradwz(-1), ringGradwx(), ringGradwz(),
                this->ringaux(-1));

        _tensor_yface(this->ringGradax(-1), this->ringGraday(-1), this->ringGradax(), this->ringGraday(),
                ringGradvx(-1), ringGradvy(-1), ringGradvz(-1),
                ringGradvx(), ringGradvy(), ringGradvz(),
                ringGradwy(-1), ringGradwz(-1), ringGradwy(), ringGradwz(),
                this->ringaux(-1));

        _tensor_zface(this->ringGradax(), this->ringGradaz(),
                ringGradvy(), ringGradvz(),
                ringGradwx(), ringGradwy(), ringGradwz(),
                this->ringaux(-1), this->ringaux(),
                this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

        this->_udot_tx(this->ringu(-1), this->ringv(-1), this->ringw(-1));
        this->_udot_ty(this->ringu(-1), this->ringv(-1), this->ringw(-1));
        this->_udot_tz(this->ringu(-1), this->ringv(-1), this->ringw(-1), this->ringu(0), this->ringv(0), this->ringw(0), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());

        this->_div_dxy();
        this->_div_dz(this->ringtzx(-1), this->ringtzy(-1), this->ringtzz(-1), this->ringutz(-1), this->ringtzx(0), this->ringtzy(0), this->ringtzz(0), this->ringutz(0));

        this->_copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
    }
}

template <typename TInput, typename TInputRing>
void Diffusion_CPP_uniform<TInput,TInputRing>::hpc_info(float& flop_convert, int& traffic_convert,
        float& flop_corners, int& traffic_corner,
        float& flop_tensorface, int& traffic_tensorface,
        float& flop_tdotu, int& traffic_tdotu,
        float& flop_div, int& traffic_div,
        float& flop_copyback, int& traffic_copyback,
        size_t& footprint)
{
    SurfaceTension_CPP_uniform<TInput,TInputRing>::hpc_info(flop_convert, traffic_convert, flop_corners, traffic_corner, flop_tensorface, traffic_tensorface,
            flop_tdotu, traffic_tdotu, flop_div, traffic_div, flop_copyback, traffic_copyback, footprint);

    const int ninputs = (int)powf(_BLOCKSIZE_ + 2, 3);
    const int ntensors = 3 * (_BLOCKSIZE_ + 1) * _BLOCKSIZE_ * _BLOCKSIZE_;

    flop_convert = 13 * ninputs;
    traffic_convert = (4 + 4) * sizeof(Real) * ninputs;

    flop_corners *= 3;
    traffic_corner *= 3;

    flop_tensorface =  32 * ntensors;
    traffic_tensorface = (30 + 3) * sizeof(Real) * ntensors;

    footprint = sizeof(*this);
}

typedef Diffusion_CPP_uniform<InputSOA_ST,RingInputSOA_ST> Diffusion_CPP;

#endif /* DIFFUSION_CPP_H_EMCEKDQP */
