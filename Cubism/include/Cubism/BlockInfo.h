/*
 *  BlockInfo.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cstdlib>
#include "MeshMap.h"

CUBISM_NAMESPACE_BEGIN

struct BlockInfo
{
    long long blockID;
    void * ptrBlock;
    bool special;
    int index[3];

    double origin[3];
    double h, h_gridpoint;
    double uniform_grid_spacing[3];
    double block_extent[3];

    double* ptr_grid_spacing[3];

    bool bUniform[3];

///////////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline void pos(T p[2], int ix, int iy) const
    {
        const int I[2] = {ix, iy};
        for (int j = 0; j < 2; ++j)
        {
            T delta = 0.0;
            if (bUniform[j])
                delta = uniform_grid_spacing[j]*(I[j]+0.5);
            else
            {
                const double* const dh = ptr_grid_spacing[j];
                for (int i = 0; i < I[j]; ++i)
                    delta += dh[i];
                delta += 0.5*dh[I[j]];
            }
            p[j] = origin[j] + delta;
        }
    }

    template <typename T>
    inline void pos(T p[3], int ix, int iy, int iz) const
    {
        const int I[3] = {ix, iy, iz};
        for (int j = 0; j < 3; ++j)
        {
            T delta = 0.0;
            if (bUniform[j])
                delta = uniform_grid_spacing[j]*(I[j]+0.5);
            else
            {
                const double* const dh = ptr_grid_spacing[j];
                for (int i = 0; i < I[j]; ++i)
                    delta += dh[i];
                delta += 0.5*dh[I[j]];
            }
            p[j] = origin[j] + delta;
        }
    }

    // Return grid spacing dx in all dimensions for cell {ix, iy, iz}
    template <typename T> // 2D
    inline void spacing(T dx[2], int ix, int iy) const
    {
        const int I[2] = {ix, iy};
        for (int i = 0; i < 2; ++i)
        {
            if (bUniform[i])
                dx[i] = uniform_grid_spacing[i];
            else
                dx[i] = ptr_grid_spacing[i][I[i]];
        }
    }

    template <typename T> // 3D
    inline void spacing(T dx[3], int ix, int iy, int iz) const
    {
        const int I[3] = {ix, iy, iz};
        for (int i = 0; i < 3; ++i)
        {
            if (bUniform[i])
                dx[i] = uniform_grid_spacing[i];
            else
                dx[i] = ptr_grid_spacing[i][I[i]];
        }
    }
///////////////////////////////////////////////////////////////////////////////

    BlockInfo(long long ID, const int idx[3], const double _pos[3], const double _spacing, double h_gridpoint_, void * ptr=NULL, const bool _special=false):
    blockID(ID), ptrBlock(ptr), special(_special)
    {
        h = _spacing;
        h_gridpoint = h_gridpoint_;

        for (int i = 0; i < 3; ++i)
        {
            ptr_grid_spacing[i] = NULL;

            index[i] = idx[i];

            origin[i] = _pos[i];

            uniform_grid_spacing[i] = h_gridpoint_;

            block_extent[i] = h;

            bUniform[i] = true;
        }
    }

    template <typename TBlock>
    BlockInfo(long long ID, const int idx[3], MeshMap<TBlock>* const mapX, MeshMap<TBlock>* const mapY, MeshMap<TBlock>* const mapZ, void * ptr=NULL, const bool _special=false):
    blockID(ID), ptrBlock(ptr), special(_special)
    {
        // TODO: [fabianw@mavt.ethz.ch; Wed May 03 2017 05:06:58 PM (-0700)]
        // ugly but keep this for the moment to ensure that nothing breaks
        // down as this has propagated into uniform mesh applications.
        // WARNING: THIS CAN CAUSE UNEXPECTED RESULTS (if used with a
        // nonuniform grid)!
        h = -1.0;
        h_gridpoint = -1.0;

        MeshMap<TBlock>* const ptr_map[3] = {mapX, mapY, mapZ};
        for (int i = 0; i < 3; ++i)
        {
            index[i] = idx[i];

            origin[i] = ptr_map[i]->block_origin(idx[i]);

            block_extent[i] = ptr_map[i]->block_width(idx[i]);

            ptr_grid_spacing[i] = ptr_map[i]->get_grid_spacing(idx[i]);

            bUniform[i] = ptr_map[i]->uniform();

            const double* const dh = ptr_grid_spacing[i];
            if (bUniform[i]) uniform_grid_spacing[i] = dh[0];
            else             uniform_grid_spacing[i] = -1.0;
        }
    }

    BlockInfo():blockID(-1), ptrBlock(NULL) {}
};

CUBISM_NAMESPACE_END
