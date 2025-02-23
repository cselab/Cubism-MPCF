/*
 *  NonUniform.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger 05/08/2017
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef NONUNIFORM_H_BDLG0RET
#define NONUNIFORM_H_BDLG0RET

#include <iostream>
#include <vector>
#include <cstdlib>

#include <Cubism/BlockInfo.h>
#include <Cubism/MeshMap.h>
#include "common.h"
#include "WenoCoefficients.h"
#include "FDCoefficients.h"
using namespace cubism;


///////////////////////////////////////////////////////////////////////////////
template <typename TBlock, typename TWENO, typename TFD>
class NonUniformScheme
{
public:
    NonUniformScheme(const double xS, const double xE,
            const double yS, const double yE,
            const double zS, const double zE,
            const unsigned int nBlocksX, const unsigned int nBlocksY, const unsigned int nBlocksZ) :
        m_h_min(HUGE_VAL),
        m_initialized(false),
        m_map_x(xS,xE,nBlocksX),
        m_map_y(yS,yE,nBlocksY),
        m_map_z(zS,zE,nBlocksZ)
    {}
    ~NonUniformScheme() {}


    typedef MeshMap<TBlock> TMeshMap;

    template <int _S, int _E>
    class HaloVector
    {
    public:
        static const int START = _S;
        static const int END   = _E;

        inline double operator()(const int i) const { return m_data[i+_S]; }
        inline double& operator()(const int i) { return m_data[i+_S]; }
        inline void clear() { std::vector<double>().swap(m_data); }
        inline void fill(TMeshMap& mmap, const double* const halos)
        {
            m_data.resize(mmap.ncells() + _S + _E);
            m_data.insert(m_data.begin(), &halos[0], &halos[_S]);
            m_data.insert(m_data.begin()+_S, mmap.data_grid_spacing(), mmap.data_grid_spacing()+mmap.ncells());
            m_data.insert(m_data.begin()+_S+mmap.ncells(), &halos[_S], &halos[_S+_E]);
        }

    private:
        std::vector<double> m_data;
    };

    typedef HaloVector<3,3> TVector;


    void init(const MeshDensity* const kernel_x, const MeshDensity* const kernel_y, const MeshDensity* const kernel_z)
    {
        assert(TWENO::HALO_S <= TVector::START);
        assert(TWENO::HALO_E <= TVector::END);
        assert(TFD::HALO_S <= TVector::START);
        assert(TFD::HALO_E <= TVector::END);

        double ghosts[TVector::START + TVector::END];

        m_map_x.init(kernel_x, TVector::START, TVector::END, &ghosts[0]);
        m_all_delta_x.fill(m_map_x, &ghosts[0]);

        m_map_y.init(kernel_y, TVector::START, TVector::END, &ghosts[0]);
        m_all_delta_y.fill(m_map_y, &ghosts[0]);

        m_map_z.init(kernel_z, TVector::START, TVector::END, &ghosts[0]);
        m_all_delta_z.fill(m_map_z, &ghosts[0]);

        for (int i = 0; i < m_map_x.ncells(); ++i)
            if (m_map_x.cell_width(i) < m_h_min)
                m_h_min = m_map_x.cell_width(i);

        for (int i = 0; i < m_map_y.ncells(); ++i)
            if (m_map_y.cell_width(i) < m_h_min)
                m_h_min = m_map_y.cell_width(i);

        for (int i = 0; i < m_map_z.ncells(); ++i)
            if (m_map_z.cell_width(i) < m_h_min)
                m_h_min = m_map_z.cell_width(i);

        m_initialized = true;
    }


    void setup_coefficients(std::vector<BlockInfo>& infos)
    {
        if (!m_initialized)
        {
            std::cerr << "ERROR: NonUniformScheme: Not initialized." << std::endl;
            abort();
        }

        // non-uniform WENO scheme
        TWENO weno_x(m_map_x.ncells()+1);
        TWENO weno_y(m_map_y.ncells()+1);
        TWENO weno_z(m_map_z.ncells()+1);
        weno_x.setup(&m_all_delta_x(0), m_map_x.ncells());
        weno_y.setup(&m_all_delta_y(0), m_map_y.ncells());
        weno_z.setup(&m_all_delta_z(0), m_map_z.ncells());

        // non-uniform finite differences
        TFD fd_x(m_map_x.nblocks());
        TFD fd_y(m_map_y.nblocks());
        TFD fd_z(m_map_z.nblocks());
        fd_x.setup(&m_all_delta_x(0), m_map_x.ncells(), m_map_x.start());
        fd_y.setup(&m_all_delta_y(0), m_map_y.ncells(), m_map_y.start());
        fd_z.setup(&m_all_delta_z(0), m_map_z.ncells(), m_map_z.start());

#pragma omp parallel for
        for(int i=0; i<(int)infos.size(); ++i)
        {
            BlockInfo info = infos[i];
            TBlock& b = *(TBlock*)info.ptrBlock;

            {
                const int index = info.index[0];
                const unsigned int offset0 = TBlock::sizeX * index;
                const unsigned int offset1 = (TBlock::sizeX+1) * index;
                _set_block_weno<_BLOCKSIZE_>(weno_x, b.coeffsWENO_x[TBlock::WENO::minus], b.coeffsWENO_x[TBlock::WENO::plus], offset0);
                _set_block_fd<_BLOCKSIZE_>(fd_x, b.coeffsFD_x, offset1);
                _set_block_invh<_BLOCKSIZE_>(m_map_x.get_grid_spacing(index), &b.invh_x[0]);
            }
            {
                const int index = info.index[1];
                const unsigned int offset0 = TBlock::sizeY * index;
                const unsigned int offset1 = (TBlock::sizeY+1) * index;
                _set_block_weno<_BLOCKSIZE_>(weno_y, b.coeffsWENO_y[TBlock::WENO::minus], b.coeffsWENO_y[TBlock::WENO::plus], offset0);
                _set_block_fd<_BLOCKSIZE_>(fd_y, b.coeffsFD_y, offset1);
                _set_block_invh<_BLOCKSIZE_>(m_map_y.get_grid_spacing(index), &b.invh_y[0]);
            }
            {
                const int index = info.index[2];
                const unsigned int offset0 = TBlock::sizeZ * index;
                const unsigned int offset1 = (TBlock::sizeZ+1) * index;
                _set_block_weno<_BLOCKSIZE_>(weno_z, b.coeffsWENO_z[TBlock::WENO::minus], b.coeffsWENO_z[TBlock::WENO::plus], offset0);
                _set_block_fd<_BLOCKSIZE_>(fd_z, b.coeffsFD_z, offset1);
                _set_block_invh<_BLOCKSIZE_>(m_map_z.get_grid_spacing(index), &b.invh_z[0]);
            }
        }

        m_all_delta_x.clear();
        m_all_delta_y.clear();
        m_all_delta_z.clear();
    }

    inline const TMeshMap& get_map_x() const { return m_map_x; }
    inline const TMeshMap& get_map_y() const { return m_map_y; }
    inline const TMeshMap& get_map_z() const { return m_map_z; }
    inline TMeshMap& get_map_x() { return m_map_x; }
    inline TMeshMap& get_map_y() { return m_map_y; }
    inline TMeshMap& get_map_z() { return m_map_z; }

    inline double minimum_cell_width() const
    {
        if (!m_initialized)
        {
            std::cerr << "ERROR: NonUniform.h: minimum_cell_width() can not return m_h_min, not initialized." << std::endl;
            abort();
        }
        return m_h_min;
    }

    void print_mesh_statistics(const bool verb=true)
    {
        if (verb)
        {
            _compute_mesh_stats("x-direction", m_map_x.kernel_name(), m_map_x.data_grid_spacing(), m_map_x.ncells());
            _compute_mesh_stats("y-direction", m_map_y.kernel_name(), m_map_y.data_grid_spacing(), m_map_y.ncells());
            _compute_mesh_stats("z-direction", m_map_z.kernel_name(), m_map_z.data_grid_spacing(), m_map_z.ncells());
        }
    }

private:
    double m_h_min;
    bool m_initialized;
    TMeshMap m_map_x;
    TMeshMap m_map_y;
    TMeshMap m_map_z;

    TVector m_all_delta_x;
    TVector m_all_delta_y;
    TVector m_all_delta_z;

    template <int _BSIZE>
    void _set_block_weno(const TWENO& weno, CoeffsWENO_Block_t& bminus, CoeffsWENO_Block_t& bplus, const unsigned int offset=0)
    {
        // set WENO coefficients
        const CoeffsWENO_All_t& wminus = weno.get_coefficients_minus();
        const CoeffsWENO_All_t& wplus  = weno.get_coefficients_plus();
        for (unsigned int i = 0; i < _BSIZE+1; ++i)
        {
            // coefficients map to cell faces
            //
            // minus
            ///////////////////////////////////////////////////////////////////
            // polynomials
            bminus.c[0][i] = wminus.c.c00[offset+i];
            bminus.c[1][i] = wminus.c.c01[offset+i];
            bminus.c[2][i] = wminus.c.c02[offset+i];

            bminus.c[3][i] = wminus.c.c10[offset+i];
            bminus.c[4][i] = wminus.c.c11[offset+i];
            bminus.c[5][i] = wminus.c.c12[offset+i];

            bminus.c[6][i] = wminus.c.c20[offset+i];
            bminus.c[7][i] = wminus.c.c21[offset+i];
            bminus.c[8][i] = wminus.c.c22[offset+i];

            // ideal weights
            bminus.d[0][i] = wminus.d.d0[offset+i];
            bminus.d[1][i] = wminus.d.d1[offset+i];
            bminus.d[2][i] = wminus.d.d2[offset+i];

            // smoothness indicators
            bminus.b[0][i] = wminus.b.b00[offset+i];
            bminus.b[1][i] = wminus.b.b01[offset+i];
            bminus.b[2][i] = wminus.b.b02[offset+i];

            bminus.b[3][i] = wminus.b.b10[offset+i];
            bminus.b[4][i] = wminus.b.b11[offset+i];
            bminus.b[5][i] = wminus.b.b12[offset+i];

            bminus.b[6][i] = wminus.b.b20[offset+i];
            bminus.b[7][i] = wminus.b.b21[offset+i];
            bminus.b[8][i] = wminus.b.b22[offset+i];

            // plus
            ///////////////////////////////////////////////////////////////////
            // polynomials
            bplus.c[0][i] = wplus.c.c00[offset+i];
            bplus.c[1][i] = wplus.c.c01[offset+i];
            bplus.c[2][i] = wplus.c.c02[offset+i];

            bplus.c[3][i] = wplus.c.c10[offset+i];
            bplus.c[4][i] = wplus.c.c11[offset+i];
            bplus.c[5][i] = wplus.c.c12[offset+i];

            bplus.c[6][i] = wplus.c.c20[offset+i];
            bplus.c[7][i] = wplus.c.c21[offset+i];
            bplus.c[8][i] = wplus.c.c22[offset+i];

            // ideal weights
            bplus.d[0][i] = wplus.d.d0[offset+i];
            bplus.d[1][i] = wplus.d.d1[offset+i];
            bplus.d[2][i] = wplus.d.d2[offset+i];

            // smoothness indicators
            bplus.b[0][i] = wplus.b.b00[offset+i];
            bplus.b[1][i] = wplus.b.b01[offset+i];
            bplus.b[2][i] = wplus.b.b02[offset+i];

            bplus.b[3][i] = wplus.b.b10[offset+i];
            bplus.b[4][i] = wplus.b.b11[offset+i];
            bplus.b[5][i] = wplus.b.b12[offset+i];

            bplus.b[6][i] = wplus.b.b20[offset+i];
            bplus.b[7][i] = wplus.b.b21[offset+i];
            bplus.b[8][i] = wplus.b.b22[offset+i];
        }
    }

    template <int _BSIZE>
    void _set_block_fd(const TFD& fd, CoeffsFD_Block_t& bcoeffs, const unsigned int offset=0)
    {
        const CoeffsFD_All_t& fdcoeffs_all = fd.get_coefficients();
        for (int c = -CoeffsFD_All_t::N_COEFFS/2; c <= CoeffsFD_All_t::N_COEFFS/2; ++c)
        {
            for (unsigned int i = 0; i < _BSIZE+1; ++i)
                bcoeffs.get_coeffs(c)[i] = fdcoeffs_all.get_coeffs(c)[offset+i];
        }
    }

    template <int _BSIZE>
    void _set_block_invh(const double* const grid_spacing, Real* const invh)
    {
        for (int i = 0; i < _BSIZE; ++i)
            invh[i] = 1.0/grid_spacing[i];
    }

    void _compute_mesh_stats(const std::string header, const std::string name, const double* const data, const unsigned int N)
    {
        std::cout << name << " statistics: " << header << std::endl;
        {
            double mean = 0;
            double var = 0;
            double skew = 0;
            double kurt = 0;
            double min =  HUGE_VAL;
            double max = -HUGE_VAL;

            for (unsigned int i = 0; i < N; ++i)
            {
                if (data[i] < min) min = data[i];
                if (data[i] > max) max = data[i];
            }

            int k = 0;
            double M2 = 0, M3 = 0, M4 = 0;
            for (unsigned int i = 0; i < N; ++i)
            {
                const int k1 = k;
                ++k;
                const double delta = data[i] - mean;
                const double delta_k = delta / k;
                const double delta_k2 = delta_k * delta_k;
                const double term1 = delta * delta_k * k1;
                mean += delta_k;
                M4 += term1 * delta_k2 * (k*k - 3*k + 3) + 6 * delta_k2 * M2 - 4 * delta_k * M3;
                M3 += term1 * delta_k * (k - 2) - 3 * delta_k * M2;
                M2 += term1;
            }
            assert(k > 1);
            var  = M2 / (k - 1);
            skew = std::sqrt(k) * M3 / std::pow(M2, 1.5);
            kurt = k * M4 / (M2 * M2) - 3;
            std::cout << "\tMesh spacing:  mean=" << std::scientific << mean;
            std::cout << "; std=" << std::scientific << std::sqrt(var);
            std::cout << "; skew=" << std::scientific << skew;
            std::cout << "; kurt=" << std::scientific << kurt;
            std::cout << "; min=" << std::scientific << min;
            std::cout << "; max=" << std::scientific << max << std::endl;
        }
        {
            double mean = 0;
            double var = 0;
            double skew = 0;
            double kurt = 0;
            double min =  HUGE_VAL;
            double max = -HUGE_VAL;

            for (unsigned int i = 1; i < N; ++i)
            {
                const double r = data[i]/data[i-1];
                if (r < min) min = r;
                if (r > max) max = r;
            }

            int k = 0;
            double M2 = 0, M3 = 0, M4 = 0;
            for (unsigned int i = 1; i < N; ++i)
            {
                const double r = data[i]/data[i-1];
                const int k1 = k;
                ++k;
                const double delta = r - mean;
                const double delta_k = delta / k;
                const double delta_k2 = delta_k * delta_k;
                const double term1 = delta * delta_k * k1;
                mean += delta_k;
                M4 += term1 * delta_k2 * (k*k - 3*k + 3) + 6 * delta_k2 * M2 - 4 * delta_k * M3;
                M3 += term1 * delta_k * (k - 2) - 3 * delta_k * M2;
                M2 += term1;
            }
            assert(k > 1);
            var  = M2 / (k - 1);
            skew = std::sqrt(k) * M3 / std::pow(M2, 1.5);
            kurt = k * M4 / (M2 * M2) - 3;
            std::cout << "\tGrowth factor: mean=" << std::scientific << mean;
            std::cout << "; std=" << std::scientific << std::sqrt(var);
            std::cout << "; skew=" << std::scientific << skew;
            std::cout << "; kurt=" << std::scientific << kurt;
            std::cout << "; min=" << std::scientific << min;
            std::cout << "; max=" << std::scientific << max << std::endl;
        }
    }
};

#endif /* NONUNIFORM_H_BDLG0RET */
