/*
 *  GridMPI.h
 *  FacesMPI
 *
 *  Created by Diego Rossinelli on 10/21/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <vector>
#include <map>
#include <mpi.h>

#include "BlockInfo.h"
#include "StencilInfo.h"
#include "SynchronizerMPI.h"

CUBISM_NAMESPACE_BEGIN

template < typename TGrid >
class GridMPI : public TGrid
{
public:
    typedef typename TGrid::Real Real;
private:
    size_t timestamp;

protected:
    typedef SynchronizerMPI<Real> SynchronizerMPIType;
    friend class SynchronizerMPI<Real>;

    int myrank, mypeindex[3], pesize[3];
    int periodic[3];
    int mybpd[3], myblockstotalsize, blocksize[3];

    std::vector<BlockInfo> cached_blockinfo;

    std::map<StencilInfo, SynchronizerMPIType *> SynchronizerMPIs;

    MPI_Comm worldcomm;
    MPI_Comm cartcomm;

    // Subdomain handled by this node.
    double subdomain_low[3];
    double subdomain_high[3];
public:

    typedef typename TGrid::BlockType Block;
    typedef typename TGrid::BlockType BlockType;

    GridMPI(const int npeX, const int npeY, const int npeZ,
            const int nX, const int nY=1, const int nZ=1,
            const double _maxextent = 1, const MPI_Comm comm = MPI_COMM_WORLD):
      TGrid(nX, nY, nZ, _maxextent), timestamp(0), worldcomm(comm)
    {
        blocksize[0] = Block::sizeX;
        blocksize[1] = Block::sizeY;
        blocksize[2] = Block::sizeZ;

        mybpd[0] = nX;
        mybpd[1] = nY;
        mybpd[2] = nZ;
        myblockstotalsize = nX*nY*nZ;

        periodic[0] = true;
        periodic[1] = true;
        periodic[2] = true;

        pesize[0] = npeX;
        pesize[1] = npeY;
        pesize[2] = npeZ;

        int world_size;
        MPI_Comm_size(worldcomm, &world_size);
        assert(npeX*npeY*npeZ == world_size);

        MPI_Cart_create(worldcomm, 3, pesize, periodic, true, &cartcomm);
        MPI_Comm_rank(cartcomm, &myrank);
        MPI_Cart_coords(cartcomm, myrank, 3, mypeindex);

        const std::vector<BlockInfo> vInfo = TGrid::getBlocksInfo();

        // Doesn't make sense to export `h_gridpoint` and `h_block` as a member
        // variable + getter, as they are not fixed values in case of
        // non-uniform grids.
        const double h_gridpoint = _maxextent / (double)std::max(
                getBlocksPerDimension(0) * blocksize[0],
                std::max(getBlocksPerDimension(1) * blocksize[1],
                         getBlocksPerDimension(2) * blocksize[2]));
        const double h_block[3] = {
            blocksize[0] * h_gridpoint,
            blocksize[1] * h_gridpoint,
            blocksize[2] * h_gridpoint,
        };

        // setup uniform (global) meshmaps
        // discard single process mappings
        for (int i = 0; i < 3; ++i)
            delete this->m_mesh_maps[i];
        std::vector<MeshMap<Block>*> clearme;
        this->m_mesh_maps.swap(clearme);

        // global number of blocks and extents
        const int nBlocks[3] = {
            mybpd[0]*pesize[0],
            mybpd[1]*pesize[1],
            mybpd[2]*pesize[2]
        };
        const double extents[3] = {
            h_block[0]*nBlocks[0],
            h_block[1]*nBlocks[1],
            h_block[2]*nBlocks[2]
        };
        for (int i = 0; i < 3; ++i)
        {
            MeshMap<Block>* m = new MeshMap<Block>(0.0, extents[i], nBlocks[i]);
            UniformDensity uniform;
            m->init(&uniform); // uniform only for this constructor
            this->m_mesh_maps.push_back(m);
        }

        // This subdomain box is used by the coupling framework. Please don't
        // rearrange the formula for subdomain_high, as it has to exactly match
        // subdomain_low of the neighbouring nodes. This way the roundings are
        // guaranteed to be done in the same way. Well, at least without
        // -ffast-math. (October 2017, kicici)
        subdomain_low[0] = mypeindex[0] * mybpd[0] * h_block[0];
        subdomain_low[1] = mypeindex[1] * mybpd[1] * h_block[1];
        subdomain_low[2] = mypeindex[2] * mybpd[2] * h_block[2];
        subdomain_high[0] = (mypeindex[0] + 1) * mybpd[0] * h_block[0];
        subdomain_high[1] = (mypeindex[1] + 1) * mybpd[1] * h_block[1];
        subdomain_high[2] = (mypeindex[2] + 1) * mybpd[2] * h_block[2];

        for(size_t i=0; i<vInfo.size(); ++i)
        {
            BlockInfo info = vInfo[i];

            info.h_gridpoint = h_gridpoint;
            info.h = h_block[0];// only for blocksize[0]=blocksize[1]=blocksize[2]

            for(int j=0; j<3; ++j)
            {
                info.index[j] += mypeindex[j]*mybpd[j];

                info.origin[j] = this->m_mesh_maps[j]->block_origin(info.index[j]);

                info.uniform_grid_spacing[j] = h_gridpoint;

                info.block_extent[j] = this->m_mesh_maps[j]->block_width(info.index[j]);

                info.ptr_grid_spacing[j] = this->m_mesh_maps[j]->get_grid_spacing(info.index[j]);
            }

            cached_blockinfo.push_back(info);
        }
    }


    GridMPI(MeshMap<Block>* const mapX, MeshMap<Block>* const mapY, MeshMap<Block>* const mapZ,
            const int npeX, const int npeY, const int npeZ,
            const int nX=0, const int nY=0, const int nZ=0,
            const MPI_Comm comm = MPI_COMM_WORLD): TGrid(mapX,mapY,mapZ,nX,nY,nZ), timestamp(0), worldcomm(comm)
    {
        assert(mapX->nblocks() == static_cast<size_t>(npeX*nX));
        assert(mapY->nblocks() == static_cast<size_t>(npeY*nY));
        assert(mapZ->nblocks() == static_cast<size_t>(npeZ*nZ));

        blocksize[0] = Block::sizeX;
        blocksize[1] = Block::sizeY;
        blocksize[2] = Block::sizeZ;

        mybpd[0] = nX;
        mybpd[1] = nY;
        mybpd[2] = nZ;
        myblockstotalsize = nX*nY*nZ;

        periodic[0] = true;
        periodic[1] = true;
        periodic[2] = true;

        pesize[0] = npeX;
        pesize[1] = npeY;
        pesize[2] = npeZ;

        int world_size;
        MPI_Comm_size(worldcomm, &world_size);
        assert(npeX*npeY*npeZ == world_size);

        MPI_Cart_create(worldcomm, 3, pesize, periodic, true, &cartcomm);
        MPI_Comm_rank(cartcomm, &myrank);
        MPI_Cart_coords(cartcomm, myrank, 3, mypeindex);

        const std::vector<BlockInfo> vInfo = TGrid::getBlocksInfo();

        MeshMap<Block>* const ptr_map[3] = {mapX, mapY, mapZ};
        for(size_t i=0; i<vInfo.size(); ++i)
        {
            BlockInfo info = vInfo[i];

            for(int j=0; j<3; ++j)
            {
                info.index[j] += mypeindex[j]*mybpd[j];

                info.origin[j] = ptr_map[j]->block_origin(info.index[j]);

                info.block_extent[j] = ptr_map[j]->block_width(info.index[j]);

                info.ptr_grid_spacing[j] = ptr_map[j]->get_grid_spacing(info.index[j]);
            }

            cached_blockinfo.push_back(info);
        }
    }


    ~GridMPI()
    {
        for (auto it = SynchronizerMPIs.begin(); it != SynchronizerMPIs.end(); ++it)
            delete it->second;

        SynchronizerMPIs.clear();
        MPI_Comm_free(&cartcomm);
    }

    std::vector<BlockInfo>& getBlocksInfo() override
    {
        return cached_blockinfo;
    }

    const std::vector<BlockInfo>& getBlocksInfo() const override
    {
        return cached_blockinfo;
    }

    std::vector<BlockInfo>& getResidentBlocksInfo()
    {
        return TGrid::getBlocksInfo();
    }

    const std::vector<BlockInfo>& getResidentBlocksInfo() const
    {
        return TGrid::getBlocksInfo();
    }

    virtual bool avail(int ix, int iy=0, int iz=0) const override
    {
        //return true;
        const int originX = mypeindex[0]*mybpd[0];
        const int originY = mypeindex[1]*mybpd[1];
        const int originZ = mypeindex[2]*mybpd[2];

        const int nX = pesize[0]*mybpd[0];
        const int nY = pesize[1]*mybpd[1];
        const int nZ = pesize[2]*mybpd[2];

        ix = (ix + nX) % nX;
        iy = (iy + nY) % nY;
        iz = (iz + nZ) % nZ;

        const bool xinside = (ix>= originX && ix<nX);
        const bool yinside = (iy>= originY && iy<nY);
        const bool zinside = (iz>= originZ && iz<nZ);

        assert(TGrid::avail(ix-originX, iy-originY, iz-originZ));
        return xinside && yinside && zinside;
    }

    Block& operator()(int ix, int iy=0, int iz=0) const override
    {
        //assuming ix,iy,iz to be global
        const int originX = mypeindex[0]*mybpd[0];
        const int originY = mypeindex[1]*mybpd[1];
        const int originZ = mypeindex[2]*mybpd[2];

        const int nX = pesize[0]*mybpd[0];
        const int nY = pesize[1]*mybpd[1];
        const int nZ = pesize[2]*mybpd[2];

        ix = (ix + nX) % nX;
        iy = (iy + nY) % nY;
        iz = (iz + nZ) % nZ;

        assert(ix>= originX && ix<nX);
        assert(iy>= originY && iy<nY);
        assert(iz>= originZ && iz<nZ);

        return TGrid::operator()(ix-originX, iy-originY, iz-originZ);
    }

    template<typename Processing>
    SynchronizerMPIType& sync(Processing& p)
    {
        const StencilInfo stencil = p.stencil;
        assert(stencil.isvalid());

        SynchronizerMPIType * queryresult = NULL;

        typename std::map<StencilInfo, SynchronizerMPIType*>::iterator itSynchronizerMPI = SynchronizerMPIs.find(stencil);

        if (itSynchronizerMPI == SynchronizerMPIs.end())
        {
            queryresult = new SynchronizerMPIType(stencil, getBlocksInfo(), cartcomm, mybpd, blocksize);

            SynchronizerMPIs[stencil] = queryresult;
        }
        else  queryresult = itSynchronizerMPI->second;

        queryresult->sync(sizeof(typename Block::element_type)/sizeof(Real), sizeof(Real)>4 ? MPI_DOUBLE : MPI_FLOAT, timestamp);

        timestamp = (timestamp + 1) % 32768;

        return *queryresult;
    }

    template<typename Processing>
    const SynchronizerMPIType& get_SynchronizerMPI(Processing& p) const
    {
        assert((SynchronizerMPIs.find(p.stencil) != SynchronizerMPIs.end()));

        return *SynchronizerMPIs.find(p.stencil)->second;
    }

    int getResidentBlocksPerDimension(int idim) const
    {
        assert(idim>=0 && idim<3);
        return mybpd[idim];
    }

    int getBlocksPerDimension(int idim) const override
    {
        assert(idim>=0 && idim<3);
        return mybpd[idim]*pesize[idim];
    }

    void peindex(int _mypeindex[3]) const
    {
        for(int i=0; i<3; ++i)
            _mypeindex[i] = mypeindex[i];
    }

    size_t getTimeStamp() const
    {
        return timestamp;
    }

    MPI_Comm getCartComm() const
    {
        return cartcomm;
    }

    MPI_Comm getWorldComm() const
    {
        return worldcomm;
    }

    void getSubdomainLow(double low[3]) const {
        low[0] = subdomain_low[0];
        low[1] = subdomain_low[1];
        low[2] = subdomain_low[2];
    }

    void getSubdomainHigh(double high[3]) const {
        high[0] = subdomain_high[0];
        high[1] = subdomain_high[1];
        high[2] = subdomain_high[2];
    }
};

CUBISM_NAMESPACE_END
