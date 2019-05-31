/*
 *  GridTypes.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger on 05/10/17
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef GRIDTYPES_H_LXM9TWTK
#define GRIDTYPES_H_LXM9TWTK

#ifdef _FLOAT_PRECISION_
typedef float Real;
typedef float DumpReal;
#else
typedef double Real;
typedef double DumpReal;
#endif


struct FluidElement
{
    typedef Real RealType;

    Real alpha1rho1, alpha2rho2, ru, rv, rw, energy, alpha2, dummy;

    void clear() { alpha1rho1 = alpha2rho2 = ru = rv = rw = energy = alpha2 = dummy = static_cast<Real>(0.0); }

    void init(const Real val) { alpha1rho1 = alpha2rho2 = ru = rv = rw = energy = alpha2 = dummy = val; }

    FluidElement& operator=(const FluidElement& gp)
    {
        this->alpha1rho1 = gp.alpha1rho1;
        this->alpha2rho2 = gp.alpha2rho2;
        this->ru = gp.ru;
        this->rv = gp.rv;
        this->rw = gp.rw;
        this->energy=gp.energy;
        this->alpha2 = gp.alpha2;
        this->dummy = gp.dummy;

        return *this;
    }
};

FluidElement operator*(const Real a, FluidElement gp);
FluidElement operator+(FluidElement gpa, FluidElement gpb);
FluidElement operator-(FluidElement gpa, FluidElement gpb);


struct FluidBlock
{
    static const int sizeX = _BLOCKSIZE_;
    static const int sizeY = _BLOCKSIZE_;
    static const int sizeZ = _BLOCKSIZE_;

    static const int gptfloats = sizeof(FluidElement)/sizeof(FluidElement::RealType);

    typedef typename FluidElement::RealType RealType;
    typedef FluidElement ElementType;
    typedef FluidElement element_type;

    FluidElement __attribute__((__aligned__(_ALIGNBYTES_))) data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];

    RealType __attribute__((__aligned__(_ALIGNBYTES_))) tmp[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][gptfloats];

    void clear_data()
    {
        const int N = sizeX*sizeY*sizeZ;
        FluidElement * const e = &data[0][0][0];
        for(int i=0; i<N; ++i) e[i].clear();
    }

    void clear_tmp()
    {
        const int N = sizeX * sizeY * sizeZ * gptfloats;

        RealType * const e = &tmp[0][0][0][0];
        for(int i=0; i<N; ++i) e[i] = 0;
    }

    void clear()
    {
        clear_data();
        clear_tmp();
    }

    inline FluidElement& operator()(int ix, int iy=0, int iz=0)
    {
        assert(ix>=0 && ix<sizeX);
        assert(iy>=0 && iy<sizeY);
        assert(iz>=0 && iz<sizeZ);

        return data[iz][iy][ix];
    }

    inline const FluidElement& operator()(int ix, int iy=0, int iz=0) const
    {
        assert(ix>=0 && ix<sizeX);
        assert(iy>=0 && iy<sizeY);
        assert(iz>=0 && iz<sizeZ);

        return data[iz][iy][ix];
    }

    template <typename TStreamer>
    inline void Write(ofstream& output) const
    {
        for(int iz=0; iz<sizeZ; iz++)
            for(int iy=0; iy<sizeY; iy++)
                for(int ix=0; ix<sizeX; ix++)
                    TStreamer::operate(data[iz][iy][ix], output);
    }

    template <typename TStreamer>
    inline void Read(ifstream& input)
    {
        for(int iz=0; iz<sizeZ; iz++)
            for(int iy=0; iy<sizeY; iy++)
                for(int ix=0; ix<sizeX; ix++)
                    TStreamer::operate(input, data[iz][iy][ix]);
    }

// TODO: VOLFRAC_EQSYS_URSULA: include again
/*
    template <typename TStreamer>
    inline void minmax(RealType minval[TStreamer::channels], RealType maxval[TStreamer::channels], TStreamer streamer = TStreamer())
    {
        enum { NCHANNELS = TStreamer::channels };

        streamer.operate(data[0][0][0], minval);
        streamer.operate(data[0][0][0], maxval);

        for(int iz=0; iz<sizeZ; iz++)
            for(int iy=0; iy<sizeY; iy++)
                for(int ix=0; ix<sizeX; ix++)
                {
                    RealType tmp[NCHANNELS];

                    streamer.operate(data[iz][iy][ix], tmp);

                    for(int ic = 0; ic < NCHANNELS; ++ic)
                        minval[ic] = std::min(minval[ic], tmp[ic]);

                    for(int ic = 0; ic < NCHANNELS; ++ic)
                        maxval[ic] = std::max(maxval[ic], tmp[ic]);
                }
    }
*/
};


struct FluidBlockNonUniform
{
    static const int sizeX = _BLOCKSIZE_;
    static const int sizeY = _BLOCKSIZE_;
    static const int sizeZ = _BLOCKSIZE_;

    static const int gptfloats = sizeof(FluidElement)/sizeof(FluidElement::RealType);

    // NOTE:
    // |-  (i-1)  +|-  (i)  +|-  (i+1)  +|
    // The WENO::minus coefficients correspond to locations - in cell i.
    // The WENO::plus  coefficients correspond to locations + in cell i.
    enum WENO { minus = 0, plus };

    typedef typename FluidElement::RealType RealType;
    typedef FluidElement ElementType;
    typedef FluidElement element_type;

    FluidElement __attribute__((__aligned__(_ALIGNBYTES_))) data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];

    RealType __attribute__((__aligned__(_ALIGNBYTES_))) tmp[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][gptfloats];

    // coefficients for non-uniform schemes
    CoeffsWENO_Block_t __attribute__((__aligned__(_ALIGNBYTES_))) coeffsWENO_x[2]; // WENO: minus/plus
    CoeffsWENO_Block_t __attribute__((__aligned__(_ALIGNBYTES_))) coeffsWENO_y[2]; // WENO: minus/plus
    CoeffsWENO_Block_t __attribute__((__aligned__(_ALIGNBYTES_))) coeffsWENO_z[2]; // WENO: minus/plus
    CoeffsFD_Block_t   __attribute__((__aligned__(_ALIGNBYTES_))) coeffsFD_x;      // finite differences
    CoeffsFD_Block_t   __attribute__((__aligned__(_ALIGNBYTES_))) coeffsFD_y;      // finite differences
    CoeffsFD_Block_t   __attribute__((__aligned__(_ALIGNBYTES_))) coeffsFD_z;      // finite differences
    RealType __attribute__((__aligned__(_ALIGNBYTES_))) invh_x[_BLOCKSIZE_]; // pre-compute inverse mesh-spacings
    RealType __attribute__((__aligned__(_ALIGNBYTES_))) invh_y[_BLOCKSIZE_]; // pre-compute inverse mesh-spacings
    RealType __attribute__((__aligned__(_ALIGNBYTES_))) invh_z[_BLOCKSIZE_]; // pre-compute inverse mesh-spacings

    void clear_data()
    {
        const int N = sizeX*sizeY*sizeZ;
        FluidElement * const e = &data[0][0][0];
        for(int i=0; i<N; ++i) e[i].clear();
    }

    void clear_tmp()
    {
        const int N = sizeX * sizeY * sizeZ * gptfloats;

        RealType * const e = &tmp[0][0][0][0];
        for(int i=0; i<N; ++i) e[i] = 0;
    }

    void clear()
    {
        clear_data();
        clear_tmp();
    }

    inline FluidElement& operator()(int ix, int iy=0, int iz=0)
    {
        assert(ix>=0 && ix<sizeX);
        assert(iy>=0 && iy<sizeY);
        assert(iz>=0 && iz<sizeZ);

        return data[iz][iy][ix];
    }

    inline const FluidElement& operator()(int ix, int iy=0, int iz=0) const
    {
        assert(ix>=0 && ix<sizeX);
        assert(iy>=0 && iy<sizeY);
        assert(iz>=0 && iz<sizeZ);

        return data[iz][iy][ix];
    }
};


// TODO: [fabianw@mavt.ethz.ch; Wed May 10 2017 06:58:52 PM (-0700)] This is
// similar to the BlockLab issue where we need a typedef for the Lab type.
#ifdef _NONUNIFORM_BLOCK_
typedef FluidBlockNonUniform Block_t; // nonuniform block
#else
typedef FluidBlock Block_t;  // uniform block
#endif /* _NONUNIFORM_BLOCK_ */

typedef Grid<Block_t, std::allocator> GridBase;
typedef GridBase Grid_t;

#endif /* GRIDTYPES_H_LXM9TWTK */
