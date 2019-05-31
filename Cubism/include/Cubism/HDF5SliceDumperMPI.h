//
//  HDF5SliceDumperMPI.h
//  Cubism
//
//  Created by Fabian Wermelinger 09/28/2016
//  Copyright 2016 ETH Zurich. All rights reserved.
//
#ifndef HDF5SLICEDUMPERMPI_H_ZENQHJA6
#define HDF5SLICEDUMPERMPI_H_ZENQHJA6

#include <mpi.h>
#include "HDF5SliceDumper.h"
#include "ArgumentParser.h"

CUBISM_NAMESPACE_BEGIN

///////////////////////////////////////////////////////////////////////////////
// helpers
namespace SliceTypesMPI
{
    template <typename TGrid>
    class Slice : public SliceTypes::Slice<TGrid>
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
            return SliceTypes::Slice<TGrid>::template getEntities<TSlice>(parser, grid);
        }

    public:
        typedef TGrid GridType;

        Slice(TGrid* grid, const int id,
                const int dir,
                const int idx,
                const double frac) :
            SliceTypes::Slice<TGrid>(grid, id, dir, idx, frac)
        {
            const int localDim[3] = {
                static_cast<int>(this->m_grid->getResidentBlocksPerDimension(0)*TBlock::sizeX),
                static_cast<int>(this->m_grid->getResidentBlocksPerDimension(1)*TBlock::sizeY),
                static_cast<int>(this->m_grid->getResidentBlocksPerDimension(2)*TBlock::sizeZ)
            };

            // get MPI related dimensions and offsets
            int peIdx[3];
            this->m_grid->peindex(peIdx);
            int myStart[3], myEnd[3];
            for (int i = 0; i < 3; ++i)
            {
                myStart[i] = localDim[i]*peIdx[i];
                myEnd[i]   = myStart[i] + localDim[i];
            }

            if ( !(myStart[this->m_dir] <= this->m_idx && this->m_idx < myEnd[this->m_dir]) )
                this->m_valid = false;

            // scale index to process local index and recompute intersecting
            // blocks
            this->m_idx = this->m_idx % localDim[this->m_dir];
            std::vector<BlockInfo> clean;
            std::vector<BlockInfo> bInfo_local = this->m_grid->getResidentBlocksInfo(); // local
            this->m_intersecting_blocks.swap(clean);
            for (size_t i = 0; i < bInfo_local.size(); ++i)
            {
                const int start = bInfo_local[i].index[this->m_dir] * BS;
                if (start <= this->m_idx && this->m_idx < (start + BS))
                    this->m_intersecting_blocks.push_back(bInfo_local[i]);
            }

            if (this->m_intersecting_blocks.empty())
                this->m_valid = false;

            // local dimensions and offsets
            this->m_localWidth  = localDim[ this->m_coord_idx[0] ];
            this->m_localHeight = localDim[ this->m_coord_idx[1] ];
            m_offsetWidth = peIdx[ this->m_coord_idx[0] ] * localDim[ this->m_coord_idx[0] ];
            m_offsetHeight= peIdx[ this->m_coord_idx[1] ] * localDim[ this->m_coord_idx[1] ];
        }

        Slice(const Slice& c) = default;

        inline int offsetWidth() const { return m_offsetWidth; }
        inline int offsetHeight() const { return m_offsetHeight; }

    protected:
        typedef typename SliceTypes::Slice<TGrid>::TBlock TBlock;

        int m_offsetWidth, m_offsetHeight;
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
void DumpSliceHDF5MPI(const TSlice& slice,
                      const int stepID,
                      const typename TSlice::GridType::Real t,
                      const std::string &fname,
                      const std::string &dpath = ".",
                      const bool bXMF = true)
{
#ifdef CUBISM_USE_HDF
    typedef typename TSlice::GridType::BlockType B;

    // fname is the base filepath tail without file type extension and
    // additional identifiers
    std::ostringstream filename;
    std::ostringstream fullpath;
    filename << fname << "_slice" << slice.id();
    fullpath << dpath << "/" << filename.str();

    static const unsigned int NCHANNELS = TStreamer::NCHANNELS;
    const unsigned int width = slice.localWidth();
    const unsigned int height = slice.localHeight();

    int myRank;
    MPI_Comm comm = slice.getGrid()->getCartComm();
    MPI_Comm_rank(comm, &myRank);

    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    ///////////////////////////////////////////////////////////////////////////
    // write mesh
    std::vector<int> mesh_dims;
    std::vector<std::string> dset_name;
    dset_name.push_back("/vwidth");
    dset_name.push_back("/vheight");
    if (0 == myRank)
    {
        H5open();
        fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        file_id = H5Fcreate((fullpath.str()+".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        status = H5Pclose(fapl_id);

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

        // shutdown h5 file
        status = H5Fclose(file_id);
        H5close();
    }
    MPI_Barrier(comm);

    ///////////////////////////////////////////////////////////////////////////
    // startup file
    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    status = H5Pset_fapl_mpio(fapl_id, comm, MPI_INFO_NULL); if(status<0) H5Eprint1(stdout);
    file_id = H5Fopen((fullpath.str()+".h5").c_str(), H5F_ACC_RDWR, fapl_id);
    status = H5Pclose(fapl_id); if(status<0) H5Eprint1(stdout);

    ///////////////////////////////////////////////////////////////////////////
    // write data
    if (0 == myRank)
    {
        std::cout << "Allocating " << (width * height * NCHANNELS * sizeof(hdf5Real))/(1024.*1024.) << " MB of HDF5 slice data";
    }

    hdf5Real * array_all = NULL;

    if (slice.valid())
        array_all = new hdf5Real[width * height * NCHANNELS];

    hsize_t count[3]  = { // local count
        static_cast<hsize_t>(height),
        static_cast<hsize_t>(width),
        static_cast<hsize_t>(NCHANNELS)};
    hsize_t dims[3]   = { // global dimension
        static_cast<hsize_t>(slice.height()),
        static_cast<hsize_t>(slice.width()),
        static_cast<hsize_t>(NCHANNELS)};
    hsize_t offset[3] = { // file offset
        static_cast<hsize_t>(slice.offsetHeight()),
        static_cast<hsize_t>(slice.offsetWidth()),
        0};

    if (0 == myRank)
    {
        std::cout << " (Total  " << (dims[0] * dims[1] * dims[2] * sizeof(hdf5Real))/(1024.*1024.) << " MB)" << std::endl;
    }

    if (slice.valid())
        slice.template extract<TStreamer>(array_all);

    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

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
    if (bXMF && 0 == myRank)
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
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n", (int)dims[0], (int)dims[1], (int)dims[2], (int)sizeof(hdf5Real));
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

#endif /* HDF5SLICEDUMPERMPI_H_ZENQHJA6 */
