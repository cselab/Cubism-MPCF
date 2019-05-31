// File       : common.h
// Created    : Sun Aug 12 2018 04:33:51 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Common stuff
// Copyright 2018 ETH Zurich. All Rights Reserved.
#ifndef COMMON_H_QLCQRKJP
#define COMMON_H_QLCQRKJP

#ifdef CUBISM_TEST_DOUBLE
using MyReal = double;
#else
using MyReal = float;
#endif /* CUBISM_TEST_DOUBLE */

constexpr int BLOCK_SIZE = 16;

#include <cassert>
#include <mpi.h>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>

#include <Cubism/BlockInfo.h>

using namespace cubism;

template <typename TReal, size_t _AOSmembers=1>
struct Block
{
    typedef TReal ElementType;
    typedef TReal element_type;
    typedef TReal Real;
    static constexpr size_t sizeX   = BLOCK_SIZE;
    static constexpr size_t sizeY   = BLOCK_SIZE;
    static constexpr size_t sizeZ   = BLOCK_SIZE;
    static constexpr size_t members = _AOSmembers;

    inline void clear() {
        memset(&m_data[0][0][0][0], 0, BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE * _AOSmembers * sizeof(TReal));
    }

    inline const TReal (&operator()(const size_t ix, const size_t iy, const size_t iz) const)[_AOSmembers]
    {
        assert(ix<BLOCK_SIZE);
        assert(iy<BLOCK_SIZE);
        assert(iz<BLOCK_SIZE);
        return this->m_data[iz][iy][ix];
    }

    inline TReal (&operator()(const size_t ix, const size_t iy, const size_t iz))[_AOSmembers]
    {
        assert(ix<BLOCK_SIZE);
        assert(iy<BLOCK_SIZE);
        assert(iz<BLOCK_SIZE);
        return this->m_data[iz][iy][ix];
    }

    TReal m_data[BLOCK_SIZE][BLOCK_SIZE][BLOCK_SIZE][_AOSmembers];
};

template <size_t _comp=0>
struct Streamer
{
    static const int NCHANNELS = 1;

    template <typename TBlock, typename T>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, T output[NCHANNELS])
    {
        assert(_comp<TBlock::members);
        output[0] = b(ix,iy,iz)[_comp];
    }

    static const char * getAttributeName() { return "Scalar"; }
};

template <typename TGrid>
void set_grid_ic(TGrid* grid, const int myrank=0)
{
    using TBlock = typename TGrid::BlockType;
    std::vector<BlockInfo> infos = grid->getResidentBlocksInfo();

    const int NX = grid->getResidentBlocksPerDimension(0);
    const int NY = grid->getResidentBlocksPerDimension(1);
    const int NZ = grid->getResidentBlocksPerDimension(2);

    const double blocksize = TBlock::members*TBlock::sizeX*TBlock::sizeY*TBlock::sizeZ;
    const double offset = myrank * NX*NY*NZ * blocksize;

#pragma omp parallel for
    for (size_t i = 0; i < infos.size(); ++i)
    {
        BlockInfo info = infos[i];
        TBlock& b = *(TBlock*)info.ptrBlock;
        for (size_t iz = 0; iz < TBlock::sizeZ; ++iz)
            for (size_t iy = 0; iy < TBlock::sizeY; ++iy)
                for (size_t ix = 0; ix < TBlock::sizeX; ++ix)
                    for (size_t m = 0; m < TBlock::members; ++m)
                        b(ix,iy,iz)[m] = offset + info.blockID*blocksize + m + TBlock::members*(ix + (double)TBlock::sizeX*(iy + TBlock::sizeY*iz));
    }
}

#endif /* COMMON_H_QLCQRKJP */
