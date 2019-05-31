// File       : RestrictBlockAverage.h
// Created    : Mon Jul 10 2017 12:41:30 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Restriction operator by simple averaging
// Copyright 2017 ETH Zurich. All Rights Reserved.
#ifndef RESTRICTBLOCKAVERAGE_H_EJPTMN1I
#define RESTRICTBLOCKAVERAGE_H_EJPTMN1I

#include <cassert>
#include <vector>

#include <Cubism/BlockInfo.h>
#include "Types.h"
#include "GridOperator.h"
using namespace cubism;

template <typename TGridIn, typename TGridOut>
class RestrictBlockAverage : public GridOperator<TGridIn,TGridOut>
{
public:
    RestrictBlockAverage(ArgumentParser& p) :
        GridOperator<TGridIn,TGridOut>(p)
    {}
    virtual ~RestrictBlockAverage() {}

    virtual void operator()(const TGridIn& grid_in, TGridOut& grid_out, const bool verbose)
    {
        // 0.) checks
        typedef typename TGridIn::BlockType TBlockIn;
        typedef typename TGridOut::BlockType TBlockOut;
        assert(TBlockIn::sizeX == TBlockOut::sizeX);
        assert(TBlockIn::sizeY == TBlockOut::sizeY);
        assert(TBlockIn::sizeZ == TBlockOut::sizeZ);
        const int NX_in = grid_in.getResidentBlocksPerDimension(0);
        const int NY_in = grid_in.getResidentBlocksPerDimension(1);
        const int NZ_in = grid_in.getResidentBlocksPerDimension(2);
        const int NX_out= grid_out.getResidentBlocksPerDimension(0);
        const int NY_out= grid_out.getResidentBlocksPerDimension(1);
        const int NZ_out= grid_out.getResidentBlocksPerDimension(2);
        assert(NX_in >= NX_out); // our goal is to coarsen
        assert(NY_in >= NY_out); // our goal is to coarsen
        assert(NZ_in >= NZ_out); // our goal is to coarsen

        // 1.) assign fine to coarse map
        assert(NX_in % NX_out == 0);
        assert(NY_in % NY_out == 0);
        assert(NZ_in % NZ_out == 0);
        const int cpdx = NX_in / NX_out;
        const int cpdy = NY_in / NY_out;
        const int cpdz = NZ_in / NZ_out;
        assert(cpdx > 0 && cpdx <= TBlockOut::sizeX);
        assert(cpdy > 0 && cpdy <= TBlockOut::sizeY);
        assert(cpdz > 0 && cpdz <= TBlockOut::sizeZ);

        assert(TBlockOut::sizeX % cpdx == 0);
        assert(TBlockOut::sizeY % cpdy == 0);
        assert(TBlockOut::sizeZ % cpdz == 0);
        const int bstrideX = TBlockOut::sizeX/cpdx;
        const int bstrideY = TBlockOut::sizeY/cpdy;
        const int bstrideZ = TBlockOut::sizeZ/cpdz;

        // 2.) average
        std::vector<BlockInfo> info_in  = grid_in.getResidentBlocksInfo();
        std::vector<BlockInfo> info_out = grid_out.getResidentBlocksInfo();

#pragma omp parallel for
        for(int i=0; i<(int)info_out.size(); i++)
        {
            BlockInfo infoout = info_out[i];
            TBlockOut& bout = *(TBlockOut*)infoout.ptrBlock;
            bout.clear(); // zero data
        }

        for(int i=0; i<(int)info_in.size(); i++)
        {
            // src
            BlockInfo infoin = info_in[i];
            TBlockIn& bin = *(TBlockIn*)infoin.ptrBlock;

            // dst
            const int bix = infoin.index[0] / cpdx;
            const int biy = infoin.index[1] / cpdy;
            const int biz = infoin.index[2] / cpdz;
            const int cox = bstrideX*(infoin.index[0]%cpdx);
            const int coy = bstrideY*(infoin.index[1]%cpdy);
            const int coz = bstrideZ*(infoin.index[2]%cpdz);
            BlockInfo infoout = info_out[bix + NX_out*(biy + NY_out*biz)];
            TBlockOut& bout = *(TBlockOut*)infoout.ptrBlock;

            for(int iz=0; iz<TBlockIn::sizeZ; iz++)
                for(int iy=0; iy<TBlockIn::sizeY; iy++)
                    for(int ix=0; ix<TBlockIn::sizeX; ix++)
                    {
                        const int cix = ix / cpdx + cox;
                        const int ciy = iy / cpdy + coy;
                        const int ciz = iz / cpdz + coz;
                        bout(cix,ciy,ciz) = bout(cix,ciy,ciz) + bin(ix,iy,iz);
                    }
        }

        const Real factor = 1.0 / (cpdx * cpdy * cpdz);
#pragma omp parallel for
        for(int i=0; i<(int)info_out.size(); i++)
        {
            BlockInfo infoout = info_out[i];
            TBlockOut& bout = *(TBlockOut*)infoout.ptrBlock;
            for(int iz=0; iz<TBlockOut::sizeZ; iz++)
                for(int iy=0; iy<TBlockOut::sizeY; iy++)
                    for(int ix=0; ix<TBlockOut::sizeX; ix++)
                        bout(ix,iy,iz) = factor * bout(ix,iy,iz);
        }
    }
};

#endif /* RESTRICTBLOCKAVERAGE_H_EJPTMN1I */
