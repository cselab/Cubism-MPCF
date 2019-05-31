//
//  HDF5SubdomainDumperMPI.h
//  Cubism
//
//  Created by Fabian Wermelinger 2018-08-03
//  Copyright 2018 ETH Zurich. All rights reserved.
//
#ifndef HDF5SUBDOMAINDUMPERMPI_H_UAFPTNPL
#define HDF5SUBDOMAINDUMPERMPI_H_UAFPTNPL

#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <mpi.h>

#include "HDF5Dumper.h"

CUBISM_NAMESPACE_BEGIN

///////////////////////////////////////////////////////////////////////////////
// helpers
namespace SubdomainTypesMPI
{
    template <typename TGrid>
    class Subdomain : public SubdomainTypes::Subdomain<TGrid>
    {
    public:
        template <typename TSubdomain>
        static std::vector<TSubdomain> getEntities(ArgumentParser& parser, TGrid& grid)
        {
            return SubdomainTypes::Subdomain<TGrid>::template getEntities<TSubdomain>(parser, grid);
        }

    public:
        typedef TGrid GridType;

        // bb_start: cell index within which the bounding box start (lower left) lies
        // bb_end: cell index within which the bounding box end (upper right) lies
        Subdomain(TGrid* grid, const int id,
                const double start[3], const double end[3], const double* h[3],
                const int bb_start[3]=0, const int bb_end[3]=0) :
            SubdomainTypes::Subdomain<TGrid>(grid, id, start, end, h, bb_start, bb_end),
            m_suboffset{0}
            {
                int myrank;
                int color = static_cast<int>( this->m_valid );
                MPI_Comm comm = this->m_grid->getCartComm();
                MPI_Comm subcomm;
                MPI_Comm_rank(comm, &myrank);
                MPI_Comm_split(comm, color, myrank, &subcomm);

                int pe_coords[3];
                this->m_grid->peindex(pe_coords);

                // compute offsets
                this->m_max_size = 1;
                unsigned long g_max_size = 0;
                if (this->m_valid)
                {
                    // 1. determine dimension and create cartesian sub-communicator
                    int pe_shift[3] = {
                        pe_coords[0],
                        pe_coords[1],
                        pe_coords[2]
                    };
                    MPI_Bcast(pe_shift, 3, MPI_INT, 0, subcomm);
                    for (int i = 0; i < 3; ++i)
                        pe_coords[i] -= pe_shift[i];

                    int pe_subdims[3];
                    for (int i = 0; i < 3; ++i)
                    {
                        MPI_Allreduce(&pe_coords[i], &pe_subdims[i], 1, MPI_INT, MPI_MAX, subcomm);
                        pe_subdims[i] += 1; // shift from index to dimension space
                    }

                    MPI_Comm subcartcomm;
                    int periodic[3] = {true};
                    MPI_Cart_create(subcomm, 3, pe_subdims, periodic, false, &subcartcomm);

                    // 2. compute file offsets using reduced 1D communicators
                    int subdims[][3] = {
                        {true, false, false},
                        {false, true, false},
                        {false, false, true}
                    };
                    for (int i = 0; i < 3; ++i)
                    {
                        MPI_Comm dimcomm;
                        MPI_Cart_sub(subcartcomm, subdims[i], &dimcomm);
                        MPI_Exscan(&this->m_subcount[i], &m_suboffset[i], 1, MPI_INT, MPI_SUM, dimcomm);
                        MPI_Comm_free(&dimcomm);
                    }

                    MPI_Comm_free(&subcartcomm);

                    // 3. reduce maximum element size of subdomain to all
                    // others in the sub-communicator
                    for (int i = 0; i < 3; ++i)
                        this->m_max_size *= static_cast<unsigned long>( this->m_subcount[i] );
                    MPI_Allreduce(&(this->m_max_size), &g_max_size, 1, MPI_UNSIGNED_LONG, MPI_MAX, subcomm);
                }
                MPI_Comm_free(&subcomm);
                // 4. update maximum size globally
                MPI_Allreduce(&g_max_size, &(this->m_max_size), 1, MPI_UNSIGNED_LONG, MPI_MAX, comm);
            }

        Subdomain(const Subdomain& c) = default;

        inline const int (&offset() const)[3] { return m_suboffset; }

        virtual void show(const std::string prefix="") const
        {
            std::cout << prefix << "subdomain" << this->m_id << ":" << std::endl;
            std::cout << prefix << "ID               = " << this->m_id << std::endl;
            std::cout << prefix << "START            = (" << this->m_start[0] << ", " << this->m_start[1] << ", " << this->m_start[2] << ")" << std::endl;
            std::cout << prefix << "END              = (" << this->m_end[0] << ", " << this->m_end[1] << ", " << this->m_end[2] << ")" << std::endl;
            std::cout << prefix << "BBOX_START       = (" << this->m_bbox_start[0] << ", " << this->m_bbox_start[1] << ", " << this->m_bbox_start[2] << ")" << std::endl;
            std::cout << prefix << "BBOX_END         = (" << this->m_bbox_end[0] << ", " << this->m_bbox_end[1] << ", " << this->m_bbox_end[2] << ")" << std::endl;
            std::cout << prefix << "DIM              = (" << this->m_subdim[0] << ", " << this->m_subdim[1] << ", " << this->m_subdim[2] << ")" << std::endl;
            std::cout << prefix << "SUBDIM           = (" << this->m_subcount[0] << ", " << this->m_subcount[1] << ", " << this->m_subcount[2] << ")" << std::endl;
            std::cout << prefix << "OFFSET           = (" << this->m_suboffset[0] << ", " << this->m_suboffset[1] << ", " << this->m_suboffset[2] << ")" << std::endl;
            std::cout << prefix << "MAXSIZE          = " << this->m_max_size << std::endl;
            std::cout << prefix << "VALID            = " << this->m_valid << std::endl;
            std::cout << prefix << "NUMBER OF BLOCKS = " << this->m_intersecting_blocks.size() << std::endl;
        }

    protected:
        int m_suboffset[3];  // index offset for my subdomain
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
void DumpSubdomainHDF5MPI(const TSubdomain& subdomain,
                          const int stepID,
                          const typename TSubdomain::GridType::Real t,
                          const std::string &fname,
                          const std::string &dpath = ".",
                          const bool bXMF = true)
{
#ifdef CUBISM_USE_HDF
    typedef typename TSubdomain::GridType::BlockType B;

    int rank;

    // fname is the base filepath tail without file type extension and
    // additional identifiers
    std::ostringstream filename;
    std::ostringstream fullpath;
    filename << fname << "_subdomain" << subdomain.id();
    fullpath << dpath << "/" << filename.str();

    MPI_Comm comm = subdomain.getGrid()->getCartComm();
    MPI_Comm_rank(comm, &rank);

    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;

    ///////////////////////////////////////////////////////////////////////////
    // write mesh
    std::vector<int> mesh_dims;
    std::vector<std::string> dset_name;
    dset_name.push_back("/vx");
    dset_name.push_back("/vy");
    dset_name.push_back("/vz");
    if (0 == rank)
    {
        H5open();
        fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        file_id = H5Fcreate((fullpath.str()+".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        status = H5Pclose(fapl_id);

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
    std::vector<BlockInfo> infos_sub = subdomain.getBlocksInfo();
    static const unsigned int NCHANNELS = TStreamer::NCHANNELS;

    const unsigned int NX = subdomain.count()[0];
    const unsigned int NY = subdomain.count()[1];
    const unsigned int NZ = subdomain.count()[2];
    const unsigned int DX = subdomain.dim()[0];
    const unsigned int DY = subdomain.dim()[1];
    const unsigned int DZ = subdomain.dim()[2];

    if (rank==0)
    {
        std::cout << "Allocating " << (subdomain.max_size() * NCHANNELS * sizeof(hdf5Real))/(1024.*1024.) << " MB of HDF5 subdomain data";
        std::cout << " (Total  " << (DX * DY * DZ * NCHANNELS * sizeof(hdf5Real))/(1024.*1024.) << " MB)" << std::endl;
    }

    hsize_t count[4]  = { NZ, NY, NX, NCHANNELS };
    hsize_t dims[4]   = { DZ, DY, DX, NCHANNELS };
    hsize_t offset[4] = {
        static_cast<hsize_t>(subdomain.offset()[2]),
        static_cast<hsize_t>(subdomain.offset()[1]),
        static_cast<hsize_t>(subdomain.offset()[0]),
        0
    };

    hdf5Real * array_all = NULL;

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
                        // cell local check: continue if the cell does not
                        // intersect the subdomain bounding box.
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

                        // shift the indices to subdomain index space
                        gx -= bbox_start[0];
                        gy -= bbox_start[1];
                        gz -= bbox_start[2];

                        hdf5Real * const ptr = array_all + NCHANNELS*(gx + NX * (gy + NY * gz));

                        for(unsigned int j=0; j<NCHANNELS; ++j)
                            ptr[j] = output[j];
                    }
        }
    }

    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);

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

    if (bXMF && rank==0)
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
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",(int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3], (int)sizeof(hdf5Real));
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

#endif /* HDF5SUBDOMAINDUMPERMPI_H_UAFPTNPL */
