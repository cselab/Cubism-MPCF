// File       : main.cpp
// Created    : Mon Jul 10 2017 12:18:47 PM (+0200)
// Author     : Fabian Wermelinger
// Description: Tool for restriction / prolongation of grid
// Copyright 2017 ETH Zurich. All Rights Reserved.
#include <iostream>
#include <string>
#include <cstdlib>
#include <mpi.h>

#include "Types.h"
#include "Serialization.h"
#include "Streamer.h"

#include <Cubism/ArgumentParser.h>
#include <Cubism/Grid.h>
#include <Cubism/GridMPI.h>
#include <Cubism/ZBinDumper_MPI.h>
#ifdef _USE_HDF_
#include <Cubism/HDF5Dumper_MPI.h>
#endif /* _USE_HDF_ */

// grid operators
#include "GridOperator.h"
#include "Restriction/RestrictBlockAverage.h"
#include "Prolongation/ProlongHarten.h"
#include "Smoother/Smoother.h"

using namespace cubism;

// scalar field pass-through streamer
struct StreamerScalar
{
    static const std::string NAME;
    static const std::string EXT;
    static const int NCHANNELS = 1;
    static const int CLASS = 0;

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal output[NCHANNELS])
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);
        output[0] = el.dummy;
    }

    template <typename TBlock, typename TReal>
    static inline void operate(TBlock& b, const TReal input[NCHANNELS], const int ix, const int iy, const int iz)
    {
        typedef typename TBlock::ElementType TElement;
        TElement& el = b(ix,iy,iz);
        el.dummy = input[0];
    }

    template <typename TBlock, typename TReal>
    static inline void operate(const TBlock& b, const int ix, const int iy, const int iz, TReal *ovalue, const int field)
    {
        typedef typename TBlock::ElementType TElement;
        const TElement& el = b(ix,iy,iz);

        switch(field) {
            case 0: *ovalue = el.dummy; break;
            default: printf("unknown field\n"); abort(); break;
        }
    }

    template <typename TBlock, typename TReal>
    static inline void operate(TBlock& b, const TReal ivalue, const int ix, const int iy, const int iz, const int field)
    {
        typedef typename TBlock::ElementType TElement;
        TElement& el = b(ix,iy,iz);

        switch(field) {
            case 0:  el.dummy = ivalue; break;
            default: printf("unknown field\n"); abort(); break;
        }
    }

    static const char * getAttributeName() { return "Scalar"; }
};

const std::string StreamerScalar::NAME = "Scalar";
const std::string StreamerScalar::EXT = "-scalar";


template <typename TGrid>
class SerializerFactory
{
public:
    SerializerFactory(const bool isroot, ArgumentParser& p, Serialization<TGrid>& s) :
        m_isroot(isroot), m_parser(p), m_grid_serializer(s)
    {}

    typedef typename Serialization<TGrid>::TSerializer TSerializer;
    typedef typename Serialization<TGrid>::TDeserializer TDeserializer;

    void set_serializer()
    {
        const std::string out_streamer = m_parser("out_streamer").asString("serialization");
        TSerializer kernel = NULL;
        if (out_streamer == "serialization")
            kernel = _serializer<StreamerSerialization>();
        else if (out_streamer == "scalar")
            kernel = _serializer_h5<StreamerScalar>();
        // these are only useful if the grid has been loaded with the
        // 'StreamerSerialization' streamer, i.e., in_streamer = serialization
        else if (out_streamer == "pressure")
            kernel = _serializer_h5<StreamerPressure>();
        else if (out_streamer == "density")
            kernel = _serializer_h5<StreamerDensity>();
        else if (out_streamer == "velocity")
            kernel = _serializer_h5<StreamerVelocityVector>();
        else if (out_streamer == "energy")
            kernel = _serializer_h5<StreamerEnergy>();
        else if (out_streamer == "alpha2")
            kernel = _serializer_h5<StreamerAlpha2>();
        else if (out_streamer == "all_primitive")
            kernel = _serializer_h5<StreamerAllPrimitive>();
        else if (out_streamer == "all_conservative")
            kernel = _serializer_h5<StreamerAllConservative>();

        if (kernel == NULL)
        {
            if (m_isroot)
                std::cerr << "ERROR: NULL kernel" << std::endl;
            abort();
        }
        m_grid_serializer.set_serializer(kernel);
    }

    void set_deserializer()
    {
        const std::string in_streamer = m_parser("in_streamer").asString("serialization");
        TDeserializer kernel = NULL;
        if (in_streamer == "serialization")
            kernel = _deserializer<StreamerSerialization>();
        else if (in_streamer == "scalar")
            kernel = _deserializer<StreamerScalar>();

        if (kernel == NULL)
        {
            if (m_isroot)
                std::cerr << "ERROR: NULL kernel" << std::endl;
            abort();
        }
        m_grid_serializer.set_deserializer(kernel);
    }

private:
    const bool m_isroot;
    ArgumentParser& m_parser;
    Serialization<TGrid>& m_grid_serializer;

    template <typename TStreamer>
    inline TSerializer _serializer() const
    {
        if (m_parser("save_format").asString("zbin") == "zbin")
            return &DumpZBin_MPI<TStreamer,TGrid>;
#ifdef _USE_HDF_
        else if (m_parser("save_format").asString("zbin") == "h5")
            return &DumpHDF5_MPI<TStreamer,DumpReal,TGrid>;
#endif
        else
        {
            if (m_isroot)
                std::cerr << "ERROR: No suitable save format chosen" << std::endl;
            return NULL;
        }
    }

    template <typename TStreamer>
    inline TSerializer _serializer_h5() const
    {
#ifdef _USE_HDF_
        return &DumpHDF5_MPI<TStreamer,DumpReal,TGrid>;
#else
        if (m_isroot)
            std::cerr << "ERROR: HDF5 serializer requires compilation with -D_USE_HDF_" << std::endl;
        return NULL;
#endif
    }

    template <typename TStreamer>
    inline TDeserializer _deserializer() const
    {
        if (m_parser("restart_format").asString("zbin") == "zbin")
            return &ReadZBin_MPI<TStreamer,TGrid>;
#ifdef _USE_HDF_
        else if (m_parser("restart_format").asString("zbin") == "h5")
            return &ReadHDF5_MPI<TStreamer,DumpReal,TGrid>;
#endif
        else
        {
            if (m_isroot)
                std::cerr << "ERROR: No suitable restart format chosen" << std::endl;
            return NULL;
        }
    }
};


int main(int argc, char* argv[])
{
    int provided;
    MPI_Init_thread(&argc, (char ***)&argv, MPI_THREAD_MULTIPLE, &provided);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const bool isroot = (0==rank) ? true : false;

    ArgumentParser parser(argc, argv);
    if (parser.exist("conf"))
        parser.readFile(parser("conf").asString());

    // setup grids
    const Real extent = parser("extent").asDouble(1.0);
    const int xpesize = parser("xpesize").asInt(1);
    const int ypesize = parser("ypesize").asInt(1);
    const int zpesize = parser("zpesize").asInt(1);

    // fine blocks
    const int in_bpdx = parser("in_bpdx").asInt(1);
    const int in_bpdy = parser("in_bpdy").asInt(1);
    const int in_bpdz = parser("in_bpdz").asInt(1);

    // coarse blocks
    const int out_bpdx = parser("out_bpdx").asInt(1);
    const int out_bpdy = parser("out_bpdy").asInt(1);
    const int out_bpdz = parser("out_bpdz").asInt(1);

    typedef GridMPI< Grid<FluidBlock, std::allocator> > TGridIn;
    typedef GridMPI< Grid<FluidBlock, std::allocator> > TGridOut;

    TGridIn*  const grid_in = new TGridIn(xpesize, ypesize, zpesize, in_bpdx, in_bpdy, in_bpdz, extent, MPI_COMM_WORLD);
    TGridOut* const grid_out= new TGridOut(xpesize, ypesize, zpesize, out_bpdx, out_bpdy, out_bpdz, extent, MPI_COMM_WORLD);

    // setup I/O
    Serialization<TGridIn>  input_serializer(parser);
    SerializerFactory<TGridIn> infactory(isroot, parser, input_serializer);
    infactory.set_deserializer();

    Serialization<TGridOut> output_serializer(parser);
    SerializerFactory<TGridOut> outfactory(isroot, parser, output_serializer);
    outfactory.set_serializer();

    // get data for input grid
    int restart_id, step_id, step_lastdump, dumpcount;
    Real time, time_nextdump, time_lastdump;
    SerializationMeta meta(restart_id, step_id, step_lastdump, dumpcount, time, time_nextdump, time_lastdump);
    if (parser("in_streamer").asString("serialization") == "serialization")
        input_serializer.deserialize(*grid_in, meta, isroot);
    else
    {
        parser.set_strict_mode();
        const std::string infile = parser("in_file").asString();
        parser.unset_strict_mode();
        input_serializer.stream_from_file(*grid_in, infile, isroot);
    }

    // perform grid manipulation
    GridOperator<TGridIn,TGridOut>* gridmanip;
    const std::string operator_name = parser("operator").asString("RestrictBlockAverage");
    if (operator_name == "RestrictBlockAverage")
        gridmanip = new RestrictBlockAverage<TGridIn,TGridOut>(parser);
    else if (operator_name == "ProlongHarten")
        gridmanip = new ProlongHarten<TGridIn,TGridOut>(parser);
    else if (operator_name == "Smoother")
        gridmanip = new Smoother<TGridIn,TGridOut>(parser);
    else
    {
        if (isroot)
            std::cerr << "ERROR: Undefined operator '" << operator_name << "'" << std::endl;
        abort();
    }
    if (isroot) parser.print_args();
    (*gridmanip)(*grid_in, *grid_out, isroot);

    // save manipulations
    ++meta.restart_id;
    if (parser("out_streamer").asString("serialization") == "serialization")
        output_serializer.serialize(*grid_out, meta, isroot);
    else
        output_serializer.stream_to_file(*grid_out, 0, 0.0, parser("out_file").asString("stream_out"), isroot);

    MPI_Barrier(MPI_COMM_WORLD);

    // clean up
    delete grid_in;
    delete grid_out;
    delete gridmanip;

    MPI_Finalize();
    return 0;
}
