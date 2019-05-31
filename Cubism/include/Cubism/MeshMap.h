/*
 *  MeshMap.h
 *  Cubism
 *
 *  Created by Fabian Wermelinger on 05/03/17.
 *  Copyright 2017 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
#ifndef MESHMAP_H_UAYWTJDH
#define MESHMAP_H_UAYWTJDH

#include <cassert>
#include <cstddef>
#include <string>

#include "Cubism/Common.h"

CUBISM_NAMESPACE_BEGIN

// Mesh kernel base class
class MeshDensity
{
public:
    const bool uniform;
    MeshDensity(const bool _uniform) : uniform(_uniform) {}
    virtual ~MeshDensity() {}

    virtual void compute_spacing(const double xS, const double xE, const unsigned int ncells, double* const ary,
            const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL) const = 0;
    virtual std::string name() const = 0;
};

// Default uniform density kenrel
class UniformDensity : public MeshDensity
{
public:
    UniformDensity() : MeshDensity(true) {}

    virtual void compute_spacing(const double xS, const double xE, const unsigned int ncells, double* const ary,
            const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL) const
    {
        const double h = (xE - xS) / ncells;
        for (unsigned int i = 0; i < ncells; ++i)
            ary[i] = h;

        // ghost cells are given by ghost start (ghostS) and ghost end
        // (ghostE) and count the number of ghosts on either side (inclusive).
        // For example, for a symmetric 6-point stencil -> ghostS = 3 and
        // ghostE = 3.  ghost_spacing must provide valid memory for it.
        if (ghost_spacing)
            for (unsigned int i = 0; i < ghostS+ghostE; ++i)
                ghost_spacing[i] = h;
    }

    virtual std::string name() const { return std::string("UniformDensity"); }
};


template <typename TBlock>
class MeshMap
{
public:
    MeshMap(const double xS, const double xE, const unsigned int Nblocks) :
        m_xS(xS), m_xE(xE), m_extent(xE-xS), m_Nblocks(Nblocks),
        m_Ncells(Nblocks*TBlock::sizeX), // assumes uniform cells in all directions!
        m_uniform(true), m_initialized(false)
    {}

    ~MeshMap()
    {
        if (m_initialized)
        {
            delete[] m_grid_spacing;
            delete[] m_block_spacing;
        }
    }

    void init(const MeshDensity* const kernel, const unsigned int ghostS=0, const unsigned int ghostE=0, double* const ghost_spacing=NULL)
    {
        _alloc();

        kernel->compute_spacing(m_xS, m_xE, m_Ncells, m_grid_spacing, ghostS, ghostE, ghost_spacing);

        assert(m_Nblocks > 0);
        for (std::size_t i = 0; i < m_Nblocks; ++i)
        {
            double delta_block = 0.0;
            for (std::size_t j = 0; j < TBlock::sizeX; ++j)
                delta_block += m_grid_spacing[i*TBlock::sizeX + j];
            m_block_spacing[i] = delta_block;
        }

        m_uniform = kernel->uniform;
        m_kernel_name = kernel->name();
        m_initialized = true;
    }

    inline double start() const { return m_xS; }
    inline double end() const { return m_xE; }
    inline double extent() const { return m_extent; }
    inline unsigned int nblocks() const { return m_Nblocks; }
    inline unsigned int ncells() const { return m_Ncells; }
    inline bool uniform() const { return m_uniform; }
    inline std::string kernel_name() const { return m_kernel_name; }

    inline double cell_width(const unsigned int ix) const
    {
        assert(m_initialized && ix >= 0 && ix < m_Ncells);
        return m_grid_spacing[ix];
    }

    inline double block_width(const unsigned int bix) const
    {
        assert(m_initialized && bix >= 0 && bix < m_Nblocks);
        return m_block_spacing[bix];
    }

    inline double block_origin(const unsigned int bix) const
    {
        assert(m_initialized && bix >= 0 && bix < m_Nblocks);
        double offset = m_xS;
        for (unsigned int i = 0; i < bix; ++i)
            offset += m_block_spacing[i];
        return offset;
    }

    inline double* get_grid_spacing(const unsigned int bix)
    {
        assert(m_initialized && bix >= 0 && bix < m_Nblocks);
        return &m_grid_spacing[bix*TBlock::sizeX];
    }

    inline const double* get_grid_spacing(const unsigned int bix) const
    {
        assert(m_initialized && bix >= 0 && bix < m_Nblocks);
        return &m_grid_spacing[bix*TBlock::sizeX];
    }

    inline double* data_grid_spacing() { return m_grid_spacing; }
    inline const double* data_grid_spacing() const { return m_grid_spacing; }

private:
    const double m_xS;
    const double m_xE;
    const double m_extent;
    const unsigned int m_Nblocks;
    const unsigned int m_Ncells;

    bool m_uniform;
    bool m_initialized;
    double* m_grid_spacing;
    double* m_block_spacing;
    std::string m_kernel_name;

    inline void _alloc()
    {
        m_grid_spacing = new double[m_Ncells];
        m_block_spacing= new double[m_Nblocks];
    }
};

CUBISM_NAMESPACE_END

#endif /* MESHMAP_H_UAYWTJDH */
