//
//  HDF5SubdomainDumper.h
//  Cubism
//
//  Created by Fabian Wermelinger 2018-08-03
//  Copyright 2018 ETH Zurich. All rights reserved.
//
#ifndef HDF5SUBDOMAINDUMPER_H_3C2DKYV4
#define HDF5SUBDOMAINDUMPER_H_3C2DKYV4

#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "HDF5Dumper.h"

CUBISM_NAMESPACE_BEGIN

///////////////////////////////////////////////////////////////////////////////
// helpers
namespace SubdomainTypes
{
    template <typename TGrid>
    class Subdomain
    {
    public:
        template <typename TSubdomain>
        static std::vector<TSubdomain> getEntities(ArgumentParser& parser, TGrid& grid)
        {
            typedef typename TGrid::BlockType B;

            // Concept:
            // 1.) Extract bounding box data from parser input subdomains
            // 2.) Compute global bounding box for each subdomain

            // defines the number of subdomains in the config file:
            // nsubdomains  1 # defines one subdomain
            const size_t nSubdomains = parser("nsubdomains").asInt(0);
            std::vector<TSubdomain> subdomains;
            for (size_t i = 0; i < nSubdomains; ++i)
            {
                // 1.)
                const size_t id = i+1;
                std::ostringstream identifier;
                identifier << "subdomain" << id;
                parser.set_strict_mode();
                // each defined subdomain requires an origin and an extent:
                // subdomain1_origin  0 1.0 2 # requires 3 coordinates, separated by white space (no newlines!)
                // subdomain1_extent  1 1 1   # if subdomain exceeds simulation domain, it will be clipped.
                const std::string input_origin = parser(identifier.str()+"_origin").asString();
                const std::string input_extent = parser(identifier.str()+"_extent").asString();
                parser.unset_strict_mode();

                std::vector<double> new_origin(0);
                {
                    std::istringstream iss(input_origin);
                    for (int coord = 0; coord < 3; ++coord)
                    {
                        assert(!iss.eof());
                        double val;
                        iss >> val;
                        new_origin.push_back(val);
                    }
                }

                std::vector<double> new_extent(0);
                {
                    std::istringstream iss(input_extent);
                    for (int coord = 0; coord < 3; ++coord)
                    {
                        assert(!iss.eof());
                        double val;
                        iss >> val;
                        new_extent.push_back(val);
                    }
                }

                // 2.)
                int idx_start[3];
                int idx_end[3];
                const double* h[3];
                double sub_start[3];
                double sub_end[3];
                for (int coord = 0; coord < 3; ++coord)
                {
                    const MeshMap<B>& mmap = grid.getMeshMap(coord);
                    const double c_start = mmap.start();
                    const double c_end   = mmap.end();
                    // check and reset if domain violation.
                    assert(new_origin[coord] < c_end);
                    assert(new_extent[coord] > 0);
                    if (new_origin[coord] < c_start)
                        new_origin[coord] = c_start; // set to domain start
                    if (new_origin[coord]+new_extent[coord] > c_end)
                        new_extent[coord] = c_end - new_origin[coord]; // clip to domain end

                    // compute bounding box
                    const unsigned int ncells = mmap.ncells();
                    const double lower_bound = new_origin[coord];
                    double c_vertex = c_start;
                    for (unsigned int cell = 0; cell < ncells; ++cell)
                    {
                        c_vertex += mmap.cell_width(cell);
                        if (lower_bound < c_vertex)
                        {
                            idx_start[coord] = cell;
                            h[coord] = mmap.data_grid_spacing() + cell;
                            sub_start[coord] = c_vertex - mmap.cell_width(cell);
                            break;
                        }
                    }

                    const double upper_bound = lower_bound + new_extent[coord];
                    c_vertex = c_end;
                    for (int cell = ncells-1; cell >= 0; --cell)
                    {
                        c_vertex -= mmap.cell_width(cell);
                        if (upper_bound > c_vertex)
                        {
                            idx_end[coord] = cell;
                            sub_end[coord] = c_vertex + mmap.cell_width(cell);
                            break;
                        }
                    }
                }
                subdomains.emplace_back(&grid, id, sub_start, sub_end, h, idx_start, idx_end);
            }
            return subdomains;
        }

    public:
        typedef TGrid GridType;

        // bb_start: cell index within which the bounding box start (lower left) lies
        // bb_end: cell index within which the bounding box end (upper right) lies
        Subdomain(TGrid* grid, const int id,
                const double start[3], const double end[3], const double* h[3],
                const int bb_start[3]=0, const int bb_end[3]=0) :
            m_grid(grid), m_id(id),
            m_bbox_start{bb_start[0], bb_start[1], bb_start[2]},
            m_bbox_end{bb_end[0], bb_end[1], bb_end[2]},
            m_subcount{0},
            m_subdim{bb_end[0]-bb_start[0]+1, bb_end[1]-bb_start[1]+1, bb_end[2]-bb_start[2]+1},
            m_valid(false),
            m_grid_spacing{h[0], h[1], h[2]},
            m_start{start[0], start[1], start[2]},
            m_end{end[0], end[1], end[2]}
            {
                assert(m_grid != NULL);
                typedef typename TGrid::BlockType TBlock;

                // process span
                std::vector<BlockInfo> infos_all = grid->getBlocksInfo();
                const BlockInfo info_first = infos_all.front();
                const BlockInfo info_last  = infos_all.back();
                const int process_start[3] = {
                    static_cast<int>(info_first.index[0] * TBlock::sizeX),
                    static_cast<int>(info_first.index[1] * TBlock::sizeY),
                    static_cast<int>(info_first.index[2] * TBlock::sizeZ)
                };
                const int process_end[3] = {
                    static_cast<int>((info_last.index[0] + 1) * TBlock::sizeX - 1),
                    static_cast<int>((info_last.index[1] + 1) * TBlock::sizeY - 1),
                    static_cast<int>((info_last.index[2] + 1) * TBlock::sizeZ - 1)
                };
                bool b_intersect[3] = { false, false, false };
                for (size_t i = 0; i < 3; ++i)
                {
                    // case 1: subdomain is fully contained in this process
                    // dimension
                    if (process_start[i] <= m_bbox_start[i] && m_bbox_end[i] <= process_end[i])
                    {
                        b_intersect[i] = true;
                        continue;
                    }

                    // case 2: subdomain is partially contained in this process
                    // dimension (distributed)
                    if (process_start[i] <= m_bbox_start[i] && m_bbox_start[i] <= process_end[i] && m_bbox_end[i] > process_end[i])
                    {
                        m_bbox_end[i] = process_end[i];
                        b_intersect[i] = true;
                    }
                    else if (m_bbox_start[i] < process_start[i] && process_end[i] < m_bbox_end[i])
                    {
                        m_bbox_start[i] = process_start[i];
                        m_bbox_end[i]   = process_end[i];
                        b_intersect[i] = true;
                    }
                    else if (m_bbox_start[i] < process_start[i] && process_start[i] <= m_bbox_end[i] && m_bbox_end[i] <= process_end[i])
                    {
                        m_bbox_start[i] = process_start[i];
                        b_intersect[i] = true;
                    }
                }

                m_valid = true;
                m_max_size = 1;
                for (size_t i = 0; i < 3; ++i)
                {
                    m_subcount[i] = m_bbox_end[i] - m_bbox_start[i] + 1;
                    m_max_size *= static_cast<unsigned long>(m_subcount[i]);
                    m_valid = m_valid && b_intersect[i];
                }

                // see which blocks are needed
                if (m_valid)
                {
                    for (size_t i = 0; i < infos_all.size(); ++i)
                    {
                        const BlockInfo info = infos_all[i];
                        const int block_start[3] = {
                            static_cast<int>(info.index[0] * TBlock::sizeX),
                            static_cast<int>(info.index[1] * TBlock::sizeY),
                            static_cast<int>(info.index[2] * TBlock::sizeZ)
                        };
                        const int block_end[3] = {
                            static_cast<int>(block_start[0] + TBlock::sizeX - 1),
                            static_cast<int>(block_start[1] + TBlock::sizeY - 1),
                            static_cast<int>(block_start[2] + TBlock::sizeZ - 1)
                        };

                        const bool b_need_X = ((block_start[0] <= m_bbox_end[0]) && (block_end[0] >= m_bbox_start[0]));
                        const bool b_need_Y = ((block_start[1] <= m_bbox_end[1]) && (block_end[1] >= m_bbox_start[1]));
                        const bool b_need_Z = ((block_start[2] <= m_bbox_end[2]) && (block_end[2] >= m_bbox_start[2]));
                        if (b_need_X && b_need_Y && b_need_Z)
                            m_intersecting_blocks.push_back( info );
                    }
                }
            }

        Subdomain(const Subdomain& c) = default;

        inline int id() const { return m_id; }
        inline const int (&bbox_start() const)[3] { return m_bbox_start; }
        inline const int (&bbox_end() const)[3] { return m_bbox_end; }
        inline const int (&count() const)[3] { return m_subcount; }
        inline const int (&dim() const)[3] { return m_subdim; }
        inline int dim(const size_t i) const { assert(i<3); return m_subdim[i]; }
        inline const double (&start() const)[3] { return m_start; }
        inline double start(const size_t i) const { assert(i<3); return m_start[i]; }
        inline const double (&end() const)[3] { return m_end; }
        inline double end(const size_t i) const { assert(i<3); return m_end[i]; }
        inline const double* grid_spacing(const size_t i) const { assert(i<3); return m_grid_spacing[i]; }
        inline unsigned long max_size() const { return m_max_size; }
        inline bool valid() const { return m_valid; }
        inline const std::vector<BlockInfo>& getBlocksInfo() const { return m_intersecting_blocks; }
        inline TGrid* getGrid() const { return m_grid; }
        inline std::string name() const
        {
            std::ostringstream out;
            out << "subdomain" << m_id;
            return out.str();
        }

        virtual void show(const std::string prefix="") const
        {
            std::cout << prefix << "subdomain" << m_id << ":" << std::endl;
            std::cout << prefix << "ID               = " << m_id << std::endl;
            std::cout << prefix << "START            = (" << m_start[0] << ", " << m_start[1] << ", " << m_start[2] << ")" << std::endl;
            std::cout << prefix << "END              = (" << m_end[0] << ", " << m_end[1] << ", " << m_end[2] << ")" << std::endl;
            std::cout << prefix << "BBOX_START       = (" << m_bbox_start[0] << ", " << m_bbox_start[1] << ", " << m_bbox_start[2] << ")" << std::endl;
            std::cout << prefix << "BBOX_END         = (" << m_bbox_end[0] << ", " << m_bbox_end[1] << ", " << m_bbox_end[2] << ")" << std::endl;
            std::cout << prefix << "DIM              = (" << m_subdim[0] << ", " << m_subdim[1] << ", " << m_subdim[2] << ")" << std::endl;
            std::cout << prefix << "SUBDIM           = (" << m_subcount[0] << ", " << m_subcount[1] << ", " << m_subcount[2] << ")" << std::endl;
            std::cout << prefix << "MAXSIZE          = " << m_max_size << std::endl;
            std::cout << prefix << "VALID            = " << m_valid << std::endl;
            std::cout << prefix << "NUMBER OF BLOCKS = " << m_intersecting_blocks.size() << std::endl;
        }

    protected:

        TGrid * m_grid;
        const int m_id;
        int m_bbox_start[3]; // local start indices of bounding box
        int m_bbox_end[3];   // local end indices of bounding box
        int m_subcount[3];   // number of elements in local subdomain
        int m_subdim[3];     // number of elements in global subdomain
        unsigned long m_max_size;
        bool m_valid;

        const double* m_grid_spacing[3];
        double m_start[3]; // lower left coordinates of smallest subdomain that contains the specified origin in config file
        double m_end[3];   // upper right coordinates of smallest subdomain that contains the specified extent in config file
        std::vector<BlockInfo> m_intersecting_blocks;
    };
}

///////////////////////////////////////////////////////////////////////////////
// Dumpers
//
// The following requirements for the data TStreamer are required:
// TStreamer::NCHANNELS        : Number of data elements (1=Scalar, 3=Vector, 9=Tensor)
// TStreamer::operate          : Data access methods for read and write
// TStreamer::getAttributeName : Attribute name of the date ("Scalar", "Vector", "Tensor")
template<typename TStreamer, typename hdf5Real, typename TSubdomain>
void DumpSubdomainHDF5(const TSubdomain& subdomain,
                       const int stepID,
                       const typename TSubdomain::GridType::Real t,
                       const std::string &fname,
                       const std::string &dpath = ".",
                       const bool bXMF = true)
{
#ifdef CUBISM_USE_HDF
    typedef typename TSubdomain::GridType::BlockType B;

    // fname is the base filepath tail without file type extension and
    // additional identifiers
    std::ostringstream filename;
    std::ostringstream fullpath;
    filename << fname << "_subdomain" << subdomain.id();
    fullpath << dpath << "/" << filename.str();

    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    ///////////////////////////////////////////////////////////////////////////
    // startup file
    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    file_id = H5Fcreate((fullpath.str()+".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    status = H5Pclose(fapl_id); if(status<0) H5Eprint1(stdout);

    ///////////////////////////////////////////////////////////////////////////
    // write mesh
    std::vector<int> mesh_dims;
    std::vector<std::string> dset_name;
    dset_name.push_back("/vx");
    dset_name.push_back("/vy");
    dset_name.push_back("/vz");
    for (size_t i = 0; i < 3; ++i)
    {
        const int nCells = subdomain.dim(i);
        const double* const h = subdomain.grid_spacing(i);
        std::vector<double> vertices(nCells+1, subdomain.start(i));
        mesh_dims.push_back(vertices.size());

        for (int j = 0; j < nCells; ++j)
            vertices[j+1] = vertices[j] + h[j];;

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
    std::vector<BlockInfo> infos_sub = subdomain.getBlocksInfo();
    static const unsigned int NCHANNELS = TStreamer::NCHANNELS;

    const unsigned int NX = subdomain.count()[0];
    const unsigned int NY = subdomain.count()[1];
    const unsigned int NZ = subdomain.count()[2];

    std::cout << "Allocating " << (subdomain.max_size() * NCHANNELS * sizeof(hdf5Real))/(1024.*1024.) << " MB of HDF5 subdomain data" << std::endl;

    hdf5Real * array_all = NULL;

    hsize_t count[4]  = { NZ, NY, NX, NCHANNELS };
    hsize_t dims[4]   = { NZ, NY, NX, NCHANNELS };
    hsize_t offset[4] = {0, 0, 0, 0};

    if (subdomain.valid())
    {
        array_all = new hdf5Real[NX * NY * NZ * NCHANNELS];

        const int bbox_start[3] = {
            subdomain.bbox_start()[0],
            subdomain.bbox_start()[1],
            subdomain.bbox_start()[2]
        };
        const int bbox_end[3]   = {
            subdomain.bbox_end()[0],
            subdomain.bbox_end()[1],
            subdomain.bbox_end()[2]
        };

#pragma omp parallel for
        for(int i=0; i<(int)infos_sub.size(); i++)
        {
            BlockInfo& info = infos_sub[i];
            const B& b = *(B*)info.ptrBlock;

            const int idx[3] = { info.index[0], info.index[1], info.index[2] };

            for(int iz=0; iz<static_cast<int>(B::sizeZ); iz++)
                for(int iy=0; iy<static_cast<int>(B::sizeY); iy++)
                    for(int ix=0; ix<static_cast<int>(B::sizeX); ix++)
                    {
                        int gx = idx[0]*B::sizeX + ix;
                        int gy = idx[1]*B::sizeY + iy;
                        int gz = idx[2]*B::sizeZ + iz;
                        const bool b_containedX = (bbox_start[0] <= gx) && (gx <= bbox_end[0]);
                        const bool b_containedY = (bbox_start[1] <= gy) && (gy <= bbox_end[1]);
                        const bool b_containedZ = (bbox_start[2] <= gz) && (gz <= bbox_end[2]);
                        if (!(b_containedX && b_containedY && b_containedZ))
                            continue;

                        hdf5Real output[NCHANNELS];
                        for(unsigned int j=0; j<NCHANNELS; ++j)
                            output[j] = 0;

                        TStreamer::operate(b, ix, iy, iz, (hdf5Real*)output);

                        gx -= bbox_start[0]; // shift to process local
                        gy -= bbox_start[1]; // shift to process local
                        gz -= bbox_start[2]; // shift to process local

                        hdf5Real * const ptr = array_all + NCHANNELS*(gx + NX * (gy + NY * gz));

                        for(unsigned int j=0; j<NCHANNELS; ++j)
                            ptr[j] = output[j];
                    }
        }
    }

    fapl_id = H5Pcreate(H5P_DATASET_XFER);

    fspace_id = H5Screate_simple(4, dims, NULL);
#ifndef CUBISM_ON_FERMI
    dataset_id = H5Dcreate(file_id, "data", get_hdf5_type<hdf5Real>(), fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
    dataset_id = H5Dcreate2(file_id, "data", get_hdf5_type<hdf5Real>(), fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    fspace_id = H5Dget_space(dataset_id);

    H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    mspace_id = H5Screate_simple(4, count, NULL);

    if (!subdomain.valid())
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

    if (subdomain.valid())
        delete [] array_all;

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
        fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" Dimensions=\"%d %d %d\"/>\n\n", mesh_dims[2], mesh_dims[1], mesh_dims[0]);
        fprintf(xmf, "     <Geometry GeometryType=\"VxVyVz\">\n");
        fprintf(xmf, "       <DataItem Name=\"mesh_vx\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", mesh_dims[0]);
        fprintf(xmf, "        %s:/vx\n",(filename.str()+".h5").c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "       <DataItem Name=\"mesh_vy\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", mesh_dims[1]);
        fprintf(xmf, "        %s:/vy\n",(filename.str()+".h5").c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "       <DataItem Name=\"mesh_vz\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", mesh_dims[2]);
        fprintf(xmf, "        %s:/vz\n",(filename.str()+".h5").c_str());
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n\n");
        fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"%s\" Center=\"Cell\">\n", TStreamer::getAttributeName());
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", (int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3], (int)sizeof(hdf5Real));
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

#endif /* HDF5SUBDOMAINDUMPER_H_3C2DKYV4 */
