// File       : testSynchronizerMPI.cpp
// Created    : Sun Aug 12 2018 04:25:45 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Test program for MPI stencil synchronizer
// Copyright 2018 ETH Zurich. All Rights Reserved.
#include <iostream>

#include "../common.h"
#include "Cubism/ArgumentParser.h"
#include "Cubism/Grid.h"
#include "Cubism/GridMPI.h"
#include "Cubism/MeshMap.h"
#include "Cubism/SynchronizerMPI.h"
#include "Cubism/StencilInfo.h"
#include "Cubism/Profiler.h"

using namespace cubism;
using namespace std;

using MyBlock   = Block<MyReal,32>;
using MyGrid    = Grid<MyBlock>;
using MyGridMPI = GridMPI<MyGrid>;
#ifdef CUBISM_TEST_NONUNIFORM
using MyMeshMap = MeshMap<MyBlock>;
using MyDensity = UniformDensity;
#endif /* CUBISM_TEST_NONUNIFORM */


template<typename TKernel, typename TGrid>
static void process(TKernel& kernel, TGrid& grid, Profiler& prof, const int rank)
{
    static size_t pass = 0;
    SynchronizerMPI<MyReal> & Synch = grid.sync(kernel);

    prof.push_start("Inner");
    const vector<BlockInfo> avail0 = Synch.avail_inner();
    prof.pop_stop();

    prof.push_start("Halo");
    const vector<BlockInfo> avail1 = Synch.avail_halo();
    prof.pop_stop();

    ++pass;
    if (pass%20 == 0)
        cout << "[rank=" << rank << "] Pass " << pass << ": nInner=" << avail0.size() << "; nHalo=" << avail1.size() << endl;
}

template <typename TBlock>
struct Kernel
{
    StencilInfo stencil;

    Kernel(ArgumentParser& parser)
    {
        stencil.sx = parser("sx").asInt(-1);
        stencil.sy = parser("sy").asInt(-1);
        stencil.sz = parser("sz").asInt(-1);
        stencil.ex = parser("ex").asInt(1);
        stencil.ey = parser("ey").asInt(1);
        stencil.ez = parser("ez").asInt(1);
        stencil.tensorial = parser("tensorial").asBool(false);
        for (int i = 0; i < static_cast<int>(TBlock::members); ++i)
            stencil.selcomponents.push_back(i);
    }

    Kernel(const Kernel& c) = default;
    Kernel& operator=(const Kernel& c) = default;
};

using MyKernel = Kernel<MyBlock>;


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

    const int nPasses = parser("passes").asInt(500);
    MyKernel kernel(parser);
    Profiler prof;
    for (int i = 0; i < nPasses; ++i)
        process(kernel, *grid, prof, myrank);

    cout << std::flush;
    MPI_Barrier(comm);
    if (0==myrank)
        prof.printSummary();

    MPI_Finalize();
    return 0;
}
