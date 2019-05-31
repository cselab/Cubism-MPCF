//
//  SurfaceTension_CPP.h
//  MPCFcore
//
//  Created by Babak Hejazialhosseini on 8/2/11.
//  Revision by Fabian Wermelinger on 7/26/2017
//  Copyright 2011/2017 ETH Zurich. All rights reserved.
//
#ifndef SURFACETENSION_CPP_H_7BRNCV2R
#define SURFACETENSION_CPP_H_7BRNCV2R
#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>

#include "DivTensor_CPP.h"

typedef SOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1, float> InputSOAf_ST;
typedef SOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1> InputSOA_ST;
typedef RingSOA2D<-1, _BLOCKSIZE_+1, -1, _BLOCKSIZE_+1, 2> RingInputSOA_ST;

template <typename TInput, typename TInputRing>
class SurfaceTension_CPP_uniform : public virtual DivTensor_CPP<TInput>
{
public:
    SurfaceTension_CPP_uniform(const Real dtinvh, const Real h, const Real sigma=1.0, const Real a=1.0) :
        DivTensor_CPP<TInput>(),
        m_dtinvh(dtinvh), m_h(h), m_a(a), m_sigma(sigma)
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

    void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NT, const int NBLOCKS, const float MEASUREDTIME_, const bool bAwk=false);

protected:
    struct AssumedType { Real a1r1, a2r2, ru, rv, rw, E, A2, dummy; };

    const Real m_dtinvh;
    const Real m_h;
    const Real m_a;
    const Real m_sigma; // amplification factor of the rhs

    TInputRing ringu, ringv, ringw, ringaux; //slices for the primitive values
    RingTempSOA_ST ringGradax, ringGraday, ringGradaz; //slices for the corner-gradients

    virtual void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);

    virtual void _tensor_xface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
                               const TempSOA_ST& ny0, const TempSOA_ST& ny1,
                               const TempSOA_ST& nz0, const TempSOA_ST& nz1);

    virtual void _tensor_yface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
                               const TempSOA_ST& ny0, const TempSOA_ST& ny1,
                               const TempSOA_ST& nz0, const TempSOA_ST& nz1);

    virtual void _tensor_zface(const TempSOA_ST& nx, const TempSOA_ST& ny, const TempSOA_ST& nz,
                               TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz);

    virtual void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);

    inline void _input_next() { ringu.next(); ringv.next(); ringw.next(); ringaux.next(); }
    inline void _grad_next() { ringGradax.next(); ringGraday.next(); ringGradaz.next(); }
};


template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_uniform<TInput,TInputRing>::_convert(const Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    TInput& u = ringu.ref();
    TInput& v = ringv.ref();
    TInput& w = ringw.ref();
    TInput& l = ringaux.ref();

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

///////////////////////////////////////////////////////////////////////////////
// tensor implementation based on I*|n|/3 - nn^t
// Hint: for factor 1/3 see Hu & Adams 2006
///////////////////////////////////////////////////////////////////////////////
template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_uniform<TInput,TInputRing>::_tensor_xface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
        const TempSOA_ST& ny0, const TempSOA_ST& ny1,
        const TempSOA_ST& nz0, const TempSOA_ST& nz1)
{
#ifdef _SURFF3D_
    const Real factor = (Real)(1./3);
#elif _SURFF2D_
    const Real factor = (Real)(1./2);
#else
    const Real factor = (Real)(1.0);
#endif

    for(int iy=0; iy<TempPiXSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiXSOA_ST::NX; ix++)
        {
            const Real nx = (nx0(ix,iy)+nx0(ix,iy+1)+nx1(ix,iy)+nx1(ix,iy+1));
            const Real ny = (ny0(ix,iy)+ny0(ix,iy+1)+ny1(ix,iy)+ny1(ix,iy+1));
            const Real nz = (nz0(ix,iy)+nz0(ix,iy+1)+nz1(ix,iy)+nz1(ix,iy+1));

            const Real mag = mysqrt(nx*nx + ny*ny + nz*nz);
            const Real inv_mag = mag > static_cast<Real>(1.0e-9) ? ((Real)1)/mag : 0;
            assert(!std::isnan(inv_mag));

            this->txx.ref(ix, iy) = factor*mag - nx*nx*inv_mag;
            this->txy.ref(ix, iy) = -nx*ny*inv_mag;
            this->txz.ref(ix, iy) = -nx*nz*inv_mag;
        }
}

template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_uniform<TInput,TInputRing>::_tensor_yface(const TempSOA_ST& nx0, const TempSOA_ST& nx1,
        const TempSOA_ST& ny0, const TempSOA_ST& ny1,
        const TempSOA_ST& nz0, const TempSOA_ST& nz1)
{
#ifdef _SURFF3D_
    const Real factor = (Real)(1./3);
#elif _SURFF2D_
    const Real factor = (Real)(1./2);
#else
    const Real factor = (Real)(1.0);
#endif

    for(int iy=0; iy<TempPiYSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiYSOA_ST::NX; ix++)
        {
            const Real nx = (nx0(ix,iy)+nx0(ix+1,iy)+nx1(ix,iy)+nx1(ix+1,iy));
            const Real ny = (ny0(ix,iy)+ny0(ix+1,iy)+ny1(ix,iy)+ny1(ix+1,iy));
            const Real nz = (nz0(ix,iy)+nz0(ix+1,iy)+nz1(ix,iy)+nz1(ix+1,iy));

            const Real mag = mysqrt(nx*nx + ny*ny + nz*nz);
            const Real inv_mag = mag > static_cast<Real>(1.0e-9) ? ((Real)1)/mag : 0;
            assert(!std::isnan(inv_mag));

            this->tyx.ref(ix, iy) = -ny*nx*inv_mag ;
            this->tyy.ref(ix, iy) = factor*mag - ny*ny*inv_mag;
            this->tyz.ref(ix, iy) = -ny*nz*inv_mag;
        }
}

template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_uniform<TInput,TInputRing>::_tensor_zface(const TempSOA_ST& nx0, const TempSOA_ST& ny0, const TempSOA_ST& nz0,
        TempPiZSOA_ST& tzx, TempPiZSOA_ST& tzy, TempPiZSOA_ST& tzz)
{
#ifdef _SURFF3D_
    const Real factor = (Real)(1./3);
#elif _SURFF2D_
    const Real factor = (Real)(1./2);
#else
    const Real factor = (Real)(1.0);
#endif

    for(int iy=0; iy<TempPiZSOA_ST::NY; iy++)
        for(int ix=0; ix<TempPiZSOA_ST::NX; ix++)
        {
            const Real nx = (nx0(ix,iy)+nx0(ix+1,iy)+nx0(ix,iy+1)+nx0(ix+1,iy+1));
            const Real ny = (ny0(ix,iy)+ny0(ix+1,iy)+ny0(ix,iy+1)+ny0(ix+1,iy+1));
            const Real nz = (nz0(ix,iy)+nz0(ix+1,iy)+nz0(ix,iy+1)+nz0(ix+1,iy+1));

            const Real mag = mysqrt(nx*nx + ny*ny + nz*nz);
            const Real inv_mag = mag > static_cast<Real>(1.0e-9) ? ((Real)1)/mag : 0;
            assert(!std::isnan(inv_mag));

            tzx.ref(ix, iy) = -nz*nx*inv_mag ;
            tzy.ref(ix, iy) = -nz*ny*inv_mag;
            tzz.ref(ix, iy) = factor*mag - nz*nz*inv_mag;
        }
}

template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_uniform<TInput,TInputRing>::_copyback(Real * const gptfirst, const int gptfloats, const int rowgpts)
{
    const Real factor = m_sigma*m_dtinvh/(16*m_h);

    for(int iy=0; iy<OutputSOA::NY; iy++)
        for(int ix=0; ix<OutputSOA::NX; ix++)
        {
            assert(!std::isnan(this->rhsu(ix, iy)));
            assert(!std::isnan(this->rhsv(ix, iy)));
            assert(!std::isnan(this->rhsw(ix, iy)));
            assert(!std::isnan(this->rhse(ix, iy)));

            AssumedType& rhs = *(AssumedType*)(gptfirst + gptfloats*(ix + iy*rowgpts));

            rhs.ru  = m_a*rhs.ru + factor*this->rhsu(ix, iy);
            rhs.rv  = m_a*rhs.rv + factor*this->rhsv(ix, iy);
            rhs.rw  = m_a*rhs.rw + factor*this->rhsw(ix, iy);
            rhs.E   = m_a*rhs.E  + 0.5*factor*this->rhse(ix, iy);
        }
}

template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_uniform<TInput,TInputRing>::compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
        Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
{
    _convert(srcfirst, srcfloats, rowsrcs);
    _input_next();

    _convert(srcfirst + srcfloats*slicesrcs, srcfloats, rowsrcs);

    this->_corners(ringaux(-1), ringaux(0), ringGradax.ref(), ringGraday.ref(), ringGradaz.ref());
    _tensor_zface(ringGradax(), ringGraday(), ringGradaz(), this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());
    this->_udot_tz(ringu(-1), ringv(-1), ringw(-1),  ringu(0), ringv(0), ringw(0), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());

    for(int islice=0; islice<_BLOCKSIZE_; islice++)
    {
        this->_tensors_next();
        _input_next();
        _grad_next();

        _convert(srcfirst + (islice+2)*srcfloats*slicesrcs, srcfloats, rowsrcs);
        this->_corners(ringaux(-1), ringaux(0), ringGradax.ref(), ringGraday.ref(), ringGradaz.ref());

        _tensor_xface(ringGradax(0), ringGradax(1), ringGraday(0), ringGraday(1), ringGradaz(0), ringGradaz(1));
        _tensor_yface(ringGradax(0), ringGradax(1), ringGraday(0), ringGraday(1), ringGradaz(0), ringGradaz(1));
        _tensor_zface(ringGradax(), ringGraday(), ringGradaz(), this->ringtzx.ref(), this->ringtzy.ref(), this->ringtzz.ref());

        this->_udot_tx(ringu(-1), ringv(-1), ringw(-1));
        this->_udot_ty(ringu(-1), ringv(-1), ringw(-1));
        this->_udot_tz(ringu(-1), ringv(-1), ringw(-1), ringu(0), ringv(0), ringw(0), this->ringtzx(), this->ringtzy(), this->ringtzz(), this->ringutz.ref());

        this->_div_dxy();
        this->_div_dz(this->ringtzx(-1), this->ringtzy(-1), this->ringtzz(-1), this->ringutz(-1), this->ringtzx(0), this->ringtzy(0), this->ringtzz(0), this->ringutz(0));

        _copyback(dstfirst + islice*dstfloats*slicedsts, dstfloats, rowdsts);
    }
}

template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_uniform<TInput,TInputRing>::hpc_info(float& flop_convert, int& traffic_convert,
        float& flop_corners, int& traffic_corner,
        float& flop_tensorface, int& traffic_tensorface,
        float& flop_tdotu, int& traffic_tdotu,
        float& flop_div, int& traffic_div,
        float& flop_copyback, int& traffic_copyback,
        size_t& footprint)
{
    DivTensor_CPP<TInput>::hpc_info(flop_corners, traffic_corner, flop_tdotu, traffic_tdotu, flop_div, traffic_div, footprint);

    const int ninputs = (int)powf(_BLOCKSIZE_ + 2, 3);
    const int ntensors = 3 * (_BLOCKSIZE_ + 1) * _BLOCKSIZE_ * _BLOCKSIZE_;
    const int ncells = (int)powf(_BLOCKSIZE_, 3);

    flop_convert = 5 * ninputs;
    traffic_convert = (4 + 4) * sizeof(Real) * ninputs;
    flop_tensorface = 23 * ntensors;
    traffic_tensorface = (12 + 3) * sizeof(Real) * ntensors;
    flop_copyback = 3 * 4 * ncells;
    traffic_copyback = (4 + 4) * sizeof(Real) * ncells;

    footprint = sizeof(*this);
}

template <typename TInput, typename TInputRing>
void SurfaceTension_CPP_uniform<TInput,TInputRing>::printflops(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NT, const int NBLOCKS, const float MEASUREDTIME_, const bool bAwk)
{
    const float MEASUREDTIME = MEASUREDTIME_;
    const float PEAKPERF = PEAKPERF_CORE*NCORES;

    float flop_convert, flop_corners, flop_tensorface, flop_tdotu, flop_div, flop_copyback;
    int traffic_convert, traffic_corner, traffic_tensorface, traffic_tdotu, traffic_div, traffic_copyback;
    size_t footprint;

    this->hpc_info(flop_convert, traffic_convert, flop_corners, traffic_corner, flop_tensorface, traffic_tensorface,
            flop_tdotu, traffic_tdotu, flop_div, traffic_div, flop_copyback, traffic_copyback, footprint);

    double texpected_ai = 0;
    double totflop = 0;

    //compute texpected_ai
    {
        std::vector<float> ai(6), flop(6);

        ai[0] = flop_convert/traffic_convert;
        ai[1] = flop_corners/traffic_corner;
        ai[2] = flop_tensorface/traffic_tensorface;
        ai[3] = flop_tdotu/traffic_tdotu;
        ai[4] = flop_div/traffic_div;
        ai[5] = flop_copyback/traffic_copyback;

        flop[0] = flop_convert;
        flop[1] = flop_corners;
        flop[2] = flop_tensorface;
        flop[3] = flop_tdotu;
        flop[4] = flop_div;
        flop[5] = flop_copyback;

        for(int i=0; i<ai.size(); ++i)
            texpected_ai += NT * NBLOCKS * flop[i] / std::min(PEAKPERF, PEAKBAND*ai[i]);

        for(int i=0; i<ai.size(); ++i)
            totflop += NT * NBLOCKS * flop[i];
    }

    const double ai_overall = std::min((double)PEAKPERF , (totflop / texpected_ai) / PEAKBAND);

    const double inout_footprint =  NT * NBLOCKS * (5 * sizeof(Real) * powf(_BLOCKSIZE_+2, 3) + 2 * 4 * sizeof(Real) * powf(_BLOCKSIZE_, 3));

    const double oi_overall = totflop/(inout_footprint + NT * (2 + 1) * footprint);
    const double texpected_oi = totflop/std::min((double)PEAKPERF, PEAKBAND*oi_overall);

    const double perf_measured = 1e-9*totflop/MEASUREDTIME;

    printPerformanceTitle();
    printf("\tINTERMEDIATE MEMORY FOOTPRINT: %.4f MB\tTOTAL TRAFFIC: %.4f MB\n", footprint/1024./1024, (NT * NBLOCKS * footprint + inout_footprint)/1024./1024);
    printf("\tASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
    printf("\tRIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
    printf("\tDIVTENSOR THIS ONE IS %.2f GFLOP/s,\t\"per block\" %.2f FLOP/B [AI] - %.2f FLOP/B [OI]\n", perf_measured, ai_overall, oi_overall);
    printf("\tTIME PER BLOCK: %.5f ms (expected %.5f [AI] - %.5f [OI] ms)\n",  1e3*MEASUREDTIME/(NT * NBLOCKS), 1e3*texpected_ai/(NT * NBLOCKS), 1e3*texpected_oi/(NT * NBLOCKS));
    printf("\tExpected Performance is: %.2f  GFLOP/s [AI], %.2f  GFLOP/s [OI]\n", totflop*1e-9/texpected_ai, totflop*1e-9/texpected_oi);
    printf("\tEFFICIENCY: %.2f%% [AI] - %.2f%% [OI], HW-UTILIZATION: %.2f%%\n", 100.*std::min(1., texpected_ai/MEASUREDTIME), 100.*texpected_oi/MEASUREDTIME, 100*perf_measured*1e9/PEAKPERF);
    printEndLine();
}

typedef SurfaceTension_CPP_uniform<InputSOA_ST, RingInputSOA_ST> SurfaceTension_CPP;

#endif /* SURFACETENSION_CPP_H_7BRNCV2R */
