// File       : Smoother.h
// Created    : Sun Oct 29 2017 12:36:21 PM (+0100)
// Author     : Fabian Wermelinger
// Description: Smooth data
// Copyright 2017 ETH Zurich. All Rights Reserved.
#ifndef SMOOTHER_H_2PBGUG6D
#define SMOOTHER_H_2PBGUG6D

#include <cassert>
#include <vector>

#include <Cubism/BlockInfo.h>
#include <Cubism/BlockLabMPI.h>
#include "Types.h"
#include "GridOperator.h"
#include "Prolongation/MPI_GridTransfer.h"
#include "Tests.h"
using namespace cubism;

template <typename TGridIn, typename TGridOut>
class Smoother : public GridOperator<TGridIn,TGridOut>
{
public:
    Smoother(ArgumentParser& p) :
        GridOperator<TGridIn,TGridOut>(p)
    {}
    virtual ~Smoother() {}

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
        assert(NX_in == NX_out);
        assert(NY_in == NY_out);
        assert(NZ_in == NZ_out);

        const int smooth_iter = this->m_parser("smooth_iter").asInt(0);
        typedef BlockLabMPI<Lab> LabMPI;

        // copy over
        std::vector<BlockInfo> info_in  = grid_in.getResidentBlocksInfo();
        std::vector<BlockInfo> info_out = grid_out.getResidentBlocksInfo();
        assert(info_in.size() == info_out.size());

#pragma omp parallel for
        for(int i=0; i<(int)info_out.size(); i++)
        {
            BlockInfo infoout = info_out[i];
            TBlockOut& bout = *(TBlockOut*)infoout.ptrBlock;
            bout.clear(); // zero data
        }

#pragma omp parallel for
        for(int i=0; i<(int)info_in.size(); i++)
        {
            // src
            BlockInfo infoin = info_in[i];
            TBlockIn& bin = *(TBlockIn*)infoin.ptrBlock;

            // dst
            BlockInfo infoout = info_out[i];
            TBlockOut& bout = *(TBlockOut*)infoout.ptrBlock;

            for(int iz=0; iz<TBlockIn::sizeZ; iz++)
                for(int iy=0; iy<TBlockIn::sizeY; iy++)
                    for(int ix=0; ix<TBlockIn::sizeX; ix++)
                        bout(ix,iy,iz) = bin(ix,iy,iz);
        }

        // smooth out grid
        for (int i = 0; i < smooth_iter; ++i)
        {
            if (verbose)
                std::cout << "smoothing grid: iteration " << i+1 << std::endl;
            grid_smoother smoother;
            process<LabMPI>(smoother, grid_out, 0, 0);
        }
    }
};

#endif /* SMOOTHER_H_2PBGUG6D */
