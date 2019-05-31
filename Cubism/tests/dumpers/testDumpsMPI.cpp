// File       : testDumpsMPI.cpp
// Created    : Tue Jul 24 2018 02:05:12 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Test Cubism dumping facilities
// Copyright 2018 ETH Zurich. All Rights Reserved.
#include <vector>
#include <sstream>
#include "../common.h"

#include "Cubism/ArgumentParser.h"
#include "Cubism/Grid.h"
#include "Cubism/GridMPI.h"
#include "Cubism/MeshMap.h"

#define CUBISM_USE_HDF
#include "Cubism/HDF5Dumper.h"
#include "Cubism/HDF5Dumper_MPI.h"

#include "Cubism/HDF5SliceDumper.h"
#include "Cubism/HDF5SliceDumperMPI.h"

#include "Cubism/HDF5SubdomainDumper.h"
#include "Cubism/HDF5SubdomainDumperMPI.h"

// #include "ZBinDumper.h" // requires CubismZ
// #include "ZBinDumper_MPI.h" // requires CubismZ

#include "Cubism/PlainBinDumper_MPI.h"

using namespace cubism;
using namespace std;

// dumpers
#ifndef CUBISM_TEST_HDF5_DOUBLE_PRECISION
typedef float hdf5Real;
const string prec_string = "4byte";
#else
typedef double hdf5Real;
const string prec_string = "8byte";
#endif

using MyBlock        = Block<MyReal,1>;
using MyStreamer     = Streamer<0>;
using MyGrid         = Grid<MyBlock>;
using MyGridMPI      = GridMPI<MyGrid>;
using MySlice        = typename SliceTypes::Slice<MyGrid>;
using MySliceMPI     = typename SliceTypesMPI::Slice<MyGridMPI>;
using MySubdomain    = typename SubdomainTypes::Subdomain<MyGrid>;
using MySubdomainMPI = typename SubdomainTypesMPI::Subdomain<MyGridMPI>;
#ifdef CUBISM_TEST_NONUNIFORM
using MyMeshMap = MeshMap<MyBlock>;
using MyDensity = RandomDensity;
#endif /* CUBISM_TEST_NONUNIFORM */


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    ArgumentParser parser(argc, argv);

    // blocks
    const int bpdx = parser("bpdx").asInt(1);
    const int bpdy = parser("bpdy").asInt(bpdx);
    const int bpdz = parser("bpdz").asInt(bpdx);

    // processes
    const int ppdx = parser("ppdx").asInt(1);
    const int ppdy = parser("ppdy").asInt(ppdx);
    const int ppdz = parser("ppdz").asInt(ppdx);

#ifdef CUBISM_TEST_NONUNIFORM
    MyDensity mesh_density;
    MyMeshMap* xmap = new MyMeshMap(0, 1, ppdx * bpdx);
    MyMeshMap* ymap = new MyMeshMap(0, 1, ppdy * bpdy);
    MyMeshMap* zmap = new MyMeshMap(0, 1, ppdz * bpdz);
    xmap->init(&mesh_density);
    ymap->init(&mesh_density);
    zmap->init(&mesh_density);
    MyGridMPI* grid = new MyGridMPI(xmap, ymap, zmap, ppdx, ppdy, ppdz, bpdx, bpdy, bpdz);
#else
    MyGridMPI* grid = new MyGridMPI(ppdx, ppdy, ppdz, bpdx, bpdy, bpdz);
#endif /* CUBISM_TEST_NONUNIFORM */

    int myrank;
    const MPI_Comm comm = grid->getCartComm();
    MPI_Comm_rank(comm, &myrank);

    set_grid_ic(grid, myrank);

    ///////////////////////////////////////////////////////////////////////////
    // HDF5 full dumps
    {
        ostringstream fname;
        fname << "serial_rank" << myrank << "_hdf5_" << prec_string;
        DumpHDF5<MyStreamer, hdf5Real>(*(MyGrid*)grid, 0, 0, fname.str());
    }

    {
        ostringstream fname;
        fname << "mpi_hdf5_" << prec_string;
        DumpHDF5_MPI<MyStreamer, hdf5Real>(*grid, 0, 0, fname.str());
    }
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // HDF5 slice dumps
    {
        // set default slices
        parser("nslices").asInt(3);
        parser("slice1_direction").asInt(0);
        parser("slice1_fraction").asDouble(0.3);
        parser("slice2_direction").asInt(1);
        parser("slice2_fraction").asDouble(0.6);
        parser("slice3_direction").asInt(2);
        parser("slice3_fraction").asDouble(0.8);

        // parser.print_args();

        vector<MySlice> slices = MySlice::getEntities<MySlice>(parser, *(MyGrid*)grid);
        vector<MySliceMPI> slices_mpi = MySliceMPI::getEntities<MySliceMPI>(parser, *grid);

        for (size_t i = 0; i < slices.size(); ++i)
        {
            const MySlice& slice = slices[i];

            ostringstream fname;
            fname << "serial_rank" << myrank << "_hdf5_slice" << (i+1) << "_" << prec_string;
            DumpSliceHDF5<MyStreamer, hdf5Real>(slice, 0, 0, fname.str());
        }

        for (size_t i = 0; i < slices_mpi.size(); ++i)
        {
            const MySliceMPI& slice = slices_mpi[i];

            ostringstream fname;
            fname << "mpi_hdf5_slice" << (i+1) << "_" << prec_string;
            DumpSliceHDF5MPI<MyStreamer, hdf5Real>(slice, 0, 0, fname.str());
        }
    }
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // HDF5 subdomain dumps
    {
        // set default subdomains
        struct SubdomainDefinition
        {
            double origin[3];
            double extent[3];
            SubdomainDefinition(
                    const double ox, const double oy, const double oz,
                    const double ex, const double ey, const double ez
                    ) : origin{ox, oy, oz}, extent{ex, ey, ez}
                { }
        };
        vector<SubdomainDefinition> vSd;
        vSd.push_back( SubdomainDefinition(0.05, 0.05, 0.05, 0.9, 0.9, 0.9) );  // most of domain
        vSd.push_back( SubdomainDefinition(0.05, 0.05, 0.05, 0.15, 0.9, 0.9) ); // yz-tile
        vSd.push_back( SubdomainDefinition(0.05, 0.05, 0.05, 0.9, 0.15, 0.9) ); // xz-tile
        vSd.push_back( SubdomainDefinition(0.05, 0.05, 0.05, 0.9, 0.9, 0.15) ); // xy-tile
        vSd.push_back( SubdomainDefinition(0.05, 0.05, 0.05, 0.9, 0.15, 0.15) ); // x-stripe
        vSd.push_back( SubdomainDefinition(0.05, 0.05, 0.05, 0.15, 0.9, 0.15) ); // y-stripe
        vSd.push_back( SubdomainDefinition(0.05, 0.05, 0.05, 0.15, 0.15, 0.9) ); // z-stripe
        vSd.push_back( SubdomainDefinition(0.75, 0.75, 0.75, 0.20, 0.20, 0.20) ); // upper right small piece
        vSd.push_back( SubdomainDefinition(0.45, 0.45, 0.45, 0.10, 0.10, 0.10) ); // center piece

        double origin[3], extent[3];
        for (int i = 0; i < 3; ++i)
        {
            const double start = grid->getMeshMap(i).start();
            const double end   = grid->getMeshMap(i).end();
            origin[i] = start;
            extent[i] = end - start;
        }

        parser("nsubdomains").asInt(vSd.size());
        for (size_t i = 0; i < vSd.size(); ++i)
        {
            ostringstream so, se, sname;
            SubdomainDefinition& sd = vSd[i];
            for (int j = 0; j < 3; ++j)
            {
                so << origin[j] + sd.origin[j] * extent[j] << " ";
                se << sd.extent[j] * extent[j] << " ";
            }
            sname << "subdomain" << i+1;
            parser(sname.str() + "_origin").asString(so.str());
            parser(sname.str() + "_extent").asString(se.str());
        }

        // parser.print_args();

        vector<MySubdomain> subdomains = MySubdomain::getEntities<MySubdomain>(parser, *(MyGrid*)grid);
        vector<MySubdomainMPI> subdomains_mpi = MySubdomainMPI::getEntities<MySubdomainMPI>(parser, *grid);

        for (size_t i = 0; i < subdomains.size(); ++i)
        {
            const MySubdomain& subdomain = subdomains[i];

            ostringstream fname;
            fname << "serial_rank" << myrank << "_hdf5_subdomain" << (i+1) << "_" << prec_string;
            DumpSubdomainHDF5<MyStreamer, hdf5Real>(subdomain, 0, 0, fname.str());
        }

        for (size_t i = 0; i < subdomains_mpi.size(); ++i)
        {
            const MySubdomainMPI& subdomain = subdomains_mpi[i];

            ostringstream fname;
            fname << "mpi_hdf5_subdomain" << (i+1) << "_" << prec_string;
            DumpSubdomainHDF5MPI<MyStreamer, hdf5Real>(subdomain, 0, 0, fname.str());
        }
    }
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Plain binary dumper (MPI only)
    {
        ostringstream fname;
        fname << "mpi_plain_bin_" << prec_string;

        BlockInfo info_front = grid->getBlocksInfo().front();
        MyReal* const src = (MyReal*)info_front.ptrBlock;

        const int NX = grid->getResidentBlocksPerDimension(0);
        const int NY = grid->getResidentBlocksPerDimension(1);
        const int NZ = grid->getResidentBlocksPerDimension(2);
        const size_t bytes = NX*NY*NZ*sizeof(MyBlock);

        PlainDumpBin_MPI(comm, src, bytes, fname.str());
    }

    MPI_Finalize();
    return 0;
}
