//
//  HDF5SliceDumper.h
//  Cubism
//
//  Created by Fabian Wermelinger 09/27/2016
//  Copyright 2016 ETH Zurich. All rights reserved.
//
#ifndef HDF5SLICEDUMPER_H_QI4Y9HO7
#define HDF5SLICEDUMPER_H_QI4Y9HO7

#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "HDF5Dumper.h"
#include "ArgumentParser.h"

CUBISM_NAMESPACE_BEGIN

///////////////////////////////////////////////////////////////////////////////
// helpers
namespace SliceTypes
{
    template <typename TGrid>
    class Slice
    {
        // To generalize, check all occurences of BS and replace with
        // the appropriate X/Y/Z.
        static_assert(TGrid::BlockType::sizeX == TGrid::BlockType::sizeY
                      && TGrid::BlockType::sizeX == TGrid::BlockType::sizeZ,
                      "Only cubic block type implemented so far.");
        static constexpr int BS = TGrid::BlockType::sizeX;

    public:
        template <typename TSlice>
        static std::vector<TSlice> getEntities(ArgumentParser& parser, TGrid& grid)
        {
            const size_t nSlices = parser("nslices").asInt(0);
            std::vector<TSlice> slices;
            for (size_t i = 0; i < nSlices; ++i)
            {
                const size_t id = i+1;
                std::ostringstream identifier;
                identifier << "slice" << i+1;

                // fetch direction
                const std::string sDir = identifier.str() + "_direction";
                const int dir = parser(sDir).asInt(2); // default z-direction
                assert(dir >= 0 && dir < 3);

                // fetch location
                const std::string sIndex = identifier.str() + "_index";
                const std::string sFrac  = identifier.str() + "_fraction";
                int idx = -1;
                double frac = -1.0;
                if (parser.check(sIndex))
                    idx = parser(sIndex).asInt();
                else if (parser.check(sFrac))
                    frac = parser(sFrac).asDouble();
                else
                    frac = 0.5; // default
                slices.emplace_back(&grid, id, dir, idx, frac);
            }
            return slices;
        }


    public:
        typedef TGrid GridType;

        Slice(TGrid* grid, const int id,
                const int dir,
                const int idx,
                const double frac) :
            m_grid(grid), m_id(id), m_dir(dir)
        {
            assert(m_grid != NULL);
            assert(m_dir >= 0 && m_dir < 3);

            int Dim[3];
            Dim[0] = m_grid->getBlocksPerDimension(0)*TBlock::sizeX;
            Dim[1] = m_grid->getBlocksPerDimension(1)*TBlock::sizeY;
            Dim[2] = m_grid->getBlocksPerDimension(2)*TBlock::sizeZ;

            if (frac >= 0.0)
            {
                const int idx_low = static_cast<int>(Dim[m_dir] * frac);
                m_idx = (frac < 1.0) ? idx_low : Dim[m_dir]-1;
            }
            else if (idx >= 0)
                m_idx = idx;
            else
            {
                std::cerr << "Slice: WARNING: Ill defined slice" << m_id << "... Invalidating" << std::endl;
                m_valid = false;
                return;
            }
            assert(m_idx >= 0 && m_idx < Dim[m_dir]);
            m_valid = true;

            // define slice layout
            if (m_dir == 0) // x-normal slices
            {
                // zy-plane
                m_coord_idx[0] = 2;
                m_coord_idx[1] = 1;
            }
            else if (m_dir == 1) // y-normal slices
            {
                // zx-plane
                m_coord_idx[0] = 2;
                m_coord_idx[1] = 0;
            }
            else if (m_dir == 2) // z-normal slices
            {
                // xy-plane
                m_coord_idx[0] = 0;
                m_coord_idx[1] = 1;
            }
            m_width  = Dim[ m_coord_idx[0] ];
            m_height = Dim[ m_coord_idx[1] ];
            m_localWidth  = m_width;
            m_localHeight = m_height;
            m_max_size = static_cast<unsigned long>(m_width * m_height);

            std::vector<BlockInfo> bInfo_local = m_grid->getBlocksInfo();
            for (size_t i = 0; i < bInfo_local.size(); ++i)
            {
                const int start = bInfo_local[i].index[m_dir] * BS;
                if (start <= m_idx && m_idx < (start + BS))
                    m_intersecting_blocks.push_back(bInfo_local[i]);
            }

            if (m_intersecting_blocks.empty())
                m_valid = false;
        }

        Slice(const Slice& c) = default;

        inline int id() const { return m_id; }
        inline int coord_idx(const size_t i) const { assert(i<2); return m_coord_idx[i]; }
        inline int width() const { return m_width; }
        inline int height() const { return m_height; }
        inline int localWidth() const { return m_localWidth; }
        inline int localHeight() const { return m_localHeight; }
        inline unsigned long max_size() const { return m_max_size; }
        inline bool valid() const { return m_valid; }
        inline const std::vector<BlockInfo>& getBlocksInfo() const { return m_intersecting_blocks; }
        inline TGrid* getGrid() const { return m_grid; }
        inline std::string name() const
        {
            std::ostringstream out;
            out << "slice" << m_id;
            return out.str();
        }

        void show(const std::string prefix="") const
        {
            std::cout << prefix << "slice" << m_id << ":" << std::endl;
            std::cout << prefix << "ID               = " << m_id << std::endl;
            std::cout << prefix << "DIR              = " << m_dir << std::endl;
            std::cout << prefix << "IDX              = " << m_idx<< std::endl;
            std::cout << prefix << "COORD[0]         = " << m_coord_idx[0] << std::endl;
            std::cout << prefix << "COORD[1]         = " << m_coord_idx[1] << std::endl;
            std::cout << prefix << "WIDTH            = " << m_width << std::endl;
            std::cout << prefix << "HEIGHT           = " << m_height << std::endl;
            std::cout << prefix << "VALID            = " << m_valid << std::endl;
            std::cout << prefix << "NUMBER OF BLOCKS = " << m_intersecting_blocks.size() << std::endl;
        }

        template <typename TStreamer, typename hdf5Real>
        inline void extract(hdf5Real* const data) const
        {
            if (0 == m_dir)
                _YZ<TStreamer>(data);
            else if (1 == m_dir)
                _XZ<TStreamer>(data);
            else if (2 == m_dir)
                _YX<TStreamer>(data);
        }

    protected:
        typedef typename TGrid::BlockType TBlock;

        TGrid * m_grid;
        const int m_id;
        const int m_dir;
        int m_idx;
        int m_coord_idx[2];
        int m_width, m_height;
        int m_localWidth, m_localHeight;
        unsigned long m_max_size;
        bool m_valid;
        std::vector<BlockInfo> m_intersecting_blocks;

    private:
        template <typename TStreamer, typename hdf5Real>
        void _YZ(hdf5Real * const data) const
        {
            const int ix = m_idx % BS;
            const unsigned int NCHANNELS = TStreamer::NCHANNELS;

#pragma omp parallel for
            for(int i = 0; i < (int)m_intersecting_blocks.size(); ++i)
            {
                const BlockInfo& info = m_intersecting_blocks[i];
                const int idx[3] = {info.index[0], info.index[1], info.index[2]}; // local info
                TBlock& b = *(TBlock*)info.ptrBlock;

                for(unsigned int iz=0; iz<TBlock::sizeZ; ++iz)
                    for(unsigned int iy=0; iy<TBlock::sizeY; ++iy)
                    {
                        hdf5Real output[NCHANNELS];
                        for(unsigned int k=0; k<NCHANNELS; ++k)
                            output[k] = 0;

                        TStreamer::operate(b, ix, iy, iz, (hdf5Real*)output);

                        const unsigned int gy = idx[1]*TBlock::sizeY + iy;
                        const unsigned int gz = idx[2]*TBlock::sizeZ + iz;

                        hdf5Real * const ptr = data + NCHANNELS*(gz + m_localWidth * gy);

                        for(unsigned int k=0; k<NCHANNELS; ++k)
                            ptr[k] = output[k];
                    }
            }
        }

        template <typename TStreamer, typename hdf5Real>
        void _XZ(hdf5Real * const data) const
        {
            const int iy = m_idx % BS;
            const unsigned int NCHANNELS = TStreamer::NCHANNELS;

#pragma omp parallel for
            for(int i = 0; i < (int)m_intersecting_blocks.size(); ++i)
            {
                const BlockInfo& info = m_intersecting_blocks[i];
                const int idx[3] = {info.index[0], info.index[1], info.index[2]}; // local info
                TBlock& b = *(TBlock*)info.ptrBlock;

                for(unsigned int iz=0; iz<TBlock::sizeZ; ++iz)
                    for(unsigned int ix=0; ix<TBlock::sizeX; ++ix)
                    {
                        hdf5Real output[NCHANNELS];
                        for(unsigned int k=0; k<NCHANNELS; ++k)
                            output[k] = 0;

                        TStreamer::operate(b, ix, iy, iz, (hdf5Real*)output);

                        const unsigned int gx = idx[0]*TBlock::sizeX + ix;
                        const unsigned int gz = idx[2]*TBlock::sizeZ + iz;

                        hdf5Real * const ptr = data + NCHANNELS*(gz + m_localWidth * gx);

                        for(unsigned int k=0; k<NCHANNELS; ++k)
                            ptr[k] = output[k];
                    }
            }
        }

        template <typename TStreamer, typename hdf5Real>
        void _YX(hdf5Real * const data) const
        {
            const int iz = m_idx % BS;
            const unsigned int NCHANNELS = TStreamer::NCHANNELS;

#pragma omp parallel for
            for(int i = 0; i < (int)m_intersecting_blocks.size(); ++i)
            {
                const BlockInfo& info = m_intersecting_blocks[i];
                const int idx[3] = {info.index[0], info.index[1], info.index[2]}; // local info
                TBlock& b = *(TBlock*)info.ptrBlock;

                for(unsigned int iy=0; iy<TBlock::sizeY; ++iy)
                    for(unsigned int ix=0; ix<TBlock::sizeX; ++ix)
                    {
                        hdf5Real output[NCHANNELS];
                        for(unsigned int k=0; k<NCHANNELS; ++k)
                            output[k] = 0;

                        TStreamer::operate(b, ix, iy, iz, (hdf5Real*)output);

                        const unsigned int gx = idx[0]*TBlock::sizeX + ix;
                        const unsigned int gy = idx[1]*TBlock::sizeY + iy;

                        hdf5Real * const ptr = data + NCHANNELS*(gx + m_localWidth * gy);

                        for(unsigned int k=0; k<NCHANNELS; ++k)
                            ptr[k] = output[k];
                    }
            }
        }
    };
}

///////////////////////////////////////////////////////////////////////////////
// Dumpers
//
// The following requirements for the data TStreamer are required:
// TStreamer::NCHANNELS        : Number of data elements (1=Scalar, 3=Vector, 9=Tensor)
// TStreamer::operate          : Data access methods for read and write
// TStreamer::getAttributeName : Attribute name of the date ("Scalar", "Vector", "Tensor")
template<typename TStreamer, typename hdf5Real, typename TSlice>
void DumpSliceHDF5(const TSlice& slice,
                   const int stepID,
                   const typename TSlice::GridType::Real t,
                   const std::string &fname,
                   const std::string &dpath = ".",
                   const bool bXMF = true)
{
#ifdef CUBISM_USE_HDF
    typedef typename TSlice::GridType::BlockType B;

    static const unsigned int NCHANNELS = TStreamer::NCHANNELS;
    const unsigned int width = slice.width();
    const unsigned int height = slice.height();

    std::cout << "Allocating " << (width * height * NCHANNELS * sizeof(hdf5Real))/(1024.*1024.) << " MB of HDF5 slice data" << std::endl;;

    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    ///////////////////////////////////////////////////////////////////////////
    // startup file
    // fname is the base filepath tail without file type extension and
    // additional identifiers
    std::ostringstream filename;
    std::ostringstream fullpath;
    filename << fname << "_slice" << slice.id();
    fullpath << dpath << "/" << filename.str();

    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fcreate((fullpath.str()+".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    status = H5Pclose(fapl_id); if(status<0) H5Eprint1(stdout);

    ///////////////////////////////////////////////////////////////////////////
    // write mesh
    std::vector<int> mesh_dims;
    std::vector<std::string> dset_name;
    dset_name.push_back("/vwidth");
    dset_name.push_back("/vheight");

    for (size_t i = 0; i < 2; ++i)
    {
        const MeshMap<B>& m = slice.getGrid()->getMeshMap( slice.coord_idx(i) );
        std::vector<double> vertices(m.ncells()+1, m.start());
        mesh_dims.push_back(vertices.size());

        for (size_t j = 0; j < m.ncells(); ++j)
            vertices[j+1] = vertices[j] + m.cell_width(j);

        hsize_t dim[1] = {vertices.size()};
        fspace_id = H5Screate_simple(1, dim, NULL);
#ifndef CUBISM_ON_FERMI
        dataset_id = H5Dcreate(file_id, dset_name[i].c_str(), H5T_NATIVE_DOUBLE, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
        dataset_id = H5Dcreate2(file_id, dset_name[i].c_str(), H5T_NATIVE_DOUBLE, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vertices.data());
        status = H5Sclose(fspace_id);
        status = H5Dclose(dataset_id);
    }

    ///////////////////////////////////////////////////////////////////////////
    // write data
    hdf5Real * array_all = NULL;

    if (slice.valid())
        array_all = new hdf5Real[width * height * NCHANNELS];

    hsize_t count[3]  = { // local count
        static_cast<hsize_t>(height),
        static_cast<hsize_t>(width),
        static_cast<hsize_t>(NCHANNELS)};
    hsize_t dims[3]  = { // global dimension
        static_cast<hsize_t>(height),
        static_cast<hsize_t>(width),
        static_cast<hsize_t>(NCHANNELS)};
    hsize_t offset[3] = {0, 0, 0};

    if (slice.valid())
        slice.template extract<TStreamer>(array_all);

    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    fspace_id = H5Screate_simple(3, dims, NULL);
#ifndef CUBISM_ON_FERMI
    dataset_id = H5Dcreate(file_id, "data", get_hdf5_type<hdf5Real>(), fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    dataset_id = H5Dcreate2(file_id, "data", get_hdf5_type<hdf5Real>(), fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    fspace_id = H5Dget_space(dataset_id);
    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    mspace_id = H5Screate_simple(3, count, NULL);
    if (!slice.valid())
    {
        H5Sselect_none(fspace_id);
        H5Sselect_none(mspace_id);
    }
    status = H5Dwrite(dataset_id, get_hdf5_type<hdf5Real>(), mspace_id, fspace_id, fapl_id, array_all);
    if (status < 0) H5Eprint1(stdout);

    status = H5Sclose(mspace_id); if(status<0) H5Eprint1(stdout);
    status = H5Sclose(fspace_id); if(status<0) H5Eprint1(stdout);
    status = H5Dclose(dataset_id); if(status<0) H5Eprint1(stdout);
    status = H5Pclose(fapl_id); if(status<0) H5Eprint1(stdout);
    status = H5Fclose(file_id); if(status<0) H5Eprint1(stdout);
    H5close();

    if (slice.valid())
        delete [] array_all;

    // writing xmf wrapper
    if (bXMF)
    {
        FILE *xmf = 0;
        xmf = fopen((fullpath.str()+".xmf").c_str(), "w");
        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
        fprintf(xmf, " <Domain>\n");
        fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
        fprintf(xmf, "     <Time Value=\"%e\"/>\n\n", t);
        fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" Dimensions=\"%d %d\"/>\n\n", mesh_dims[1], mesh_dims[0]);
        fprintf(xmf, "     <Geometry GeometryType=\"VxVyVz\">\n");
        fprintf(xmf, "       <DataItem Name=\"mesh_vx\" Dimensions=\"1\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n");
        fprintf(xmf, "        %e\n", 0.0);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "       <DataItem Name=\"mesh_vy\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", mesh_dims[0]);
        fprintf(xmf, "        %s:/vwidth\n",(filename.str()+".h5").c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "       <DataItem Name=\"mesh_vz\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", mesh_dims[1]);
        fprintf(xmf, "        %s:/vheight\n",(filename.str()+".h5").c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n\n");
        fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Cell\">\n", TStreamer::getAttributeName());
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", (int)height, (int)width, (int)NCHANNELS, (int)sizeof(hdf5Real));
        fprintf(xmf, "        %s:/data\n",(filename.str()+".h5").c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
        fprintf(xmf, "   </Grid>\n");
        fprintf(xmf, " </Domain>\n");
        fprintf(xmf, "</Xdmf>\n");
        fclose(xmf);
    }

#else
#warning USE OF HDF WAS DISABLED AT COMPILE TIME
#endif
}

CUBISM_NAMESPACE_END

#endif /* HDF5SLICEDUMPER_H_QI4Y9HO7 */
