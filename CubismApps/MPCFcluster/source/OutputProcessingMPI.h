/*
 *  OutputProcessingMPI.h
 *  MPCFcluster
 *
 *  Created by Fabian Wermelinger 09/18/2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef OUTPUTPROCESSINGMPI_H_TQSCO9A3
#define OUTPUTPROCESSINGMPI_H_TQSCO9A3

#include <sstream>
#include <mpi.h>

#include <Cubism/GridMPI.h>
#include <Cubism/BlockLabMPI.h>

#include "OutputProcessing.h"
#include "BlockProcessor_MPI.h"

// MPI based dumper
#include <Cubism/ZBinDumper_MPI.h>
#include <Cubism/HDF5Dumper_MPI.h>
#include <Cubism/HDF5SliceDumperMPI.h>
#include <Cubism/HDF5SubdomainDumperMPI.h>
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
using namespace cubism;


#define __REGISTER_VP_ENTITY__(KEY, INFO, TYPE, HEAVY, TENTITY, TGRID, GRID, PARS, PROC, LAB) \
    this->_register(Item(#KEY"_rho",     new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Density", new SerializerWavelet<TGRID, StreamerDensity>(&PARS), &GRID, &PARS, NULL, StreamerDensity::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_Ux",      new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity Ux", new SerializerWavelet<TGRID, StreamerVelocity<0> >(&PARS), &GRID, &PARS, NULL, StreamerVelocity<0>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_Uy",      new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity Uy", new SerializerWavelet<TGRID, StreamerVelocity<1> >(&PARS), &GRID, &PARS, NULL, StreamerVelocity<1>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_Uz",      new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity Uz", new SerializerWavelet<TGRID, StreamerVelocity<2> >(&PARS), &GRID, &PARS, NULL, StreamerVelocity<2>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_IUI",     new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity magnitude", new SerializerWavelet<TGRID, StreamerVelocityMagnitude>(&PARS), &GRID, &PARS, NULL, StreamerVelocityMagnitude::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_divU",    new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Velocity divergence", new SerializerWavelet<TGRID, StreamerDivU>(&PARS), &GRID, &PARS, NULL, StreamerDivU::EXT, TYPE, HEAVY, StreamerDivU::CLASS, new OperatorType<TGRID,OdivU>(&PROC<LAB,OdivU,TGRID>)))); \
    this->_register(Item(#KEY"_KdivU",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Kdiv(U)", new SerializerWavelet<TGRID, StreamerKDivU>(&PARS), &GRID, &PARS, NULL, StreamerKDivU::EXT, TYPE, HEAVY, StreamerKDivU::CLASS, new OperatorType<TGRID,OdivU>(&PROC<LAB,OdivU,TGRID>)))); \
    this->_register(Item(#KEY"_E",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Total energy", new SerializerWavelet<TGRID, StreamerEnergy>(&PARS), &GRID, &PARS, NULL, StreamerEnergy::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_a2",      new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Alpha 2", new SerializerWavelet<TGRID, StreamerAlpha2>(&PARS), &GRID, &PARS, NULL, StreamerAlpha2::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_K",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": K", new SerializerWavelet<TGRID, StreamerK>(&PARS), &GRID, &PARS, NULL, StreamerK::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_p",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Pressure", new SerializerWavelet<TGRID, StreamerPressure>(&PARS), &GRID, &PARS, NULL, StreamerPressure::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_M",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Mach", new SerializerWavelet<TGRID, StreamerMach>(&PARS), &GRID, &PARS, NULL, StreamerMach::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_c",       new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Speed of sound", new SerializerWavelet<TGRID, StreamerSpeedOfSound>(&PARS), &GRID, &PARS, NULL, StreamerSpeedOfSound::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_Omegax",  new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Vorticity Omegax", new SerializerWavelet<TGRID, StreamerVorticity<0> >(&PARS), &GRID, &PARS, NULL, StreamerVorticity<0>::EXT, TYPE, HEAVY, StreamerVorticity<0>::CLASS, new OperatorType<TGRID,OVort>(&PROC<LAB,OVort,TGRID>)))); \
    this->_register(Item(#KEY"_Omegay",  new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Vorticity Omegay", new SerializerWavelet<TGRID, StreamerVorticity<1> >(&PARS), &GRID, &PARS, NULL, StreamerVorticity<1>::EXT, TYPE, HEAVY, StreamerVorticity<1>::CLASS, new OperatorType<TGRID,OVort>(&PROC<LAB,OVort,TGRID>)))); \
    this->_register(Item(#KEY"_Omegaz",  new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Vorticity Omegaz", new SerializerWavelet<TGRID, StreamerVorticity<2> >(&PARS), &GRID, &PARS, NULL, StreamerVorticity<2>::EXT, TYPE, HEAVY, StreamerVorticity<2>::CLASS, new OperatorType<TGRID,OVort>(&PROC<LAB,OVort,TGRID>)))); \
    this->_register(Item(#KEY"_IOmegaI", new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Vorticity magnitude", new SerializerWavelet<TGRID, StreamerVorticityMagnitude>(&PARS), &GRID, &PARS, NULL, StreamerVorticityMagnitude::EXT, TYPE, HEAVY, StreamerVorticityMagnitude::CLASS, new OperatorType<TGRID,OVort>(&PROC<LAB,OVort,TGRID>)))); \
    this->_register(Item(#KEY"_Qcrit",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Q-criterion", new SerializerWavelet<TGRID, StreamerQcriterion>(&PARS), &GRID, &PARS, NULL, StreamerQcriterion::EXT, TYPE, HEAVY, StreamerQcriterion::CLASS, new OperatorType<TGRID,OQcrit>(&PROC<LAB,OQcrit,TGRID>)))); \
    this->_register(Item(#KEY"_comp0",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 0", new SerializerWavelet<TGRID, StreamerAoScomponent<0> >(&PARS), &GRID, &PARS, NULL, StreamerAoScomponent<0>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp1",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 1", new SerializerWavelet<TGRID, StreamerAoScomponent<1> >(&PARS), &GRID, &PARS, NULL, StreamerAoScomponent<1>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp2",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 2", new SerializerWavelet<TGRID, StreamerAoScomponent<2> >(&PARS), &GRID, &PARS, NULL, StreamerAoScomponent<2>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp3",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 3", new SerializerWavelet<TGRID, StreamerAoScomponent<3> >(&PARS), &GRID, &PARS, NULL, StreamerAoScomponent<3>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp4",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 4", new SerializerWavelet<TGRID, StreamerAoScomponent<4> >(&PARS), &GRID, &PARS, NULL, StreamerAoScomponent<4>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp5",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 5", new SerializerWavelet<TGRID, StreamerAoScomponent<5> >(&PARS), &GRID, &PARS, NULL, StreamerAoScomponent<5>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp6",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 6", new SerializerWavelet<TGRID, StreamerAoScomponent<6> >(&PARS), &GRID, &PARS, NULL, StreamerAoScomponent<6>::EXT, TYPE, HEAVY))); \
    this->_register(Item(#KEY"_comp7",   new ProcessingElement<TENTITY,void*,TGRID>(#INFO": Cell component 7", new SerializerWavelet<TGRID, StreamerAoScomponent<7> >(&PARS), &GRID, &PARS, NULL, StreamerAoScomponent<7>::EXT, TYPE, HEAVY))); \


template <typename TGrid>
class SerializerWaveletBase
{
public:
    virtual void write(const TGrid& grid, const int step_id, const Real t, const std::string basename, const std::string path) = 0;
};

template <typename TGrid, typename TStreamer, template <typename,typename> class TCompressor=SerializerIO_WaveletCompression_MPI_SimpleBlocking>
class SerializerWavelet : public SerializerWaveletBase<TGrid>
{
public:
    SerializerWavelet(ArgumentParser* const p)
    {
        m_serializer.verbose();
        // TODO: (fabianw@mavt.ethz.ch; Wed 19 Oct 2016 12:09:47 PM CEST) We
        // need to transition to wtype=3 as default.  Need to test difference
        // between QPX/SSE of type=1 and non-vectorized type=3.
        m_serializer.set_wtype_write((*p)("wtype").asInt(1));

        // streamer specific settings
        const Real vpeps = (*p)("vpeps").asDouble(1.0e-3);
        if (TStreamer::NAME == "Alpha2")
            m_serializer.set_threshold((*p)("vpeps_a2").asDouble(vpeps));
        else if (TStreamer::NAME == "Pressure")
            m_serializer.set_threshold((*p)("vpeps_p").asDouble(1000.0*vpeps));
        else if (TStreamer::NAME == "Density")
            m_serializer.set_threshold((*p)("vpeps_rho").asDouble(vpeps));
        else if (TStreamer::NAME == "Energy")
            m_serializer.set_threshold((*p)("vpeps_E").asDouble(vpeps));
        else if (TStreamer::NAME == "Qcriterion")
            m_serializer.set_threshold((*p)("vpeps_qcrit").asDouble(vpeps));
        else if (TStreamer::NAME == "Velocity Ux")
            m_serializer.set_threshold((*p)("vpeps_Ux").asDouble(vpeps));
        else if (TStreamer::NAME == "Velocity Uy")
            m_serializer.set_threshold((*p)("vpeps_Uy").asDouble(vpeps));
        else if (TStreamer::NAME == "Velocity Uz")
            m_serializer.set_threshold((*p)("vpeps_Uz").asDouble(vpeps));
        else if (TStreamer::NAME == "Velocity magnitude")
            m_serializer.set_threshold((*p)("vpeps_IUI").asDouble(vpeps));
        else if (TStreamer::NAME == "Vorticity Ox")
            m_serializer.set_threshold((*p)("vpeps_Ox").asDouble(vpeps));
        else if (TStreamer::NAME == "Vorticity Oy")
            m_serializer.set_threshold((*p)("vpeps_Oy").asDouble(vpeps));
        else if (TStreamer::NAME == "Vorticity Oz")
            m_serializer.set_threshold((*p)("vpeps_Oz").asDouble(vpeps));
        else if (TStreamer::NAME == "Vorticity magnitude")
            m_serializer.set_threshold((*p)("vpeps_IOI").asDouble(vpeps));
        else if (TStreamer::NAME == "K")
            m_serializer.set_threshold((*p)("vpeps_K").asDouble(vpeps));
        else if (TStreamer::NAME == "divU")
            m_serializer.set_threshold((*p)("vpeps_divU").asDouble(vpeps));
        else if (TStreamer::NAME == "Mach")
            m_serializer.set_threshold((*p)("vpeps_M").asDouble(vpeps));
        else if (TStreamer::NAME == "Speed of sound")
            m_serializer.set_threshold((*p)("vpeps_c").asDouble(vpeps));
        else if (TStreamer::NAME == "KdivU")
            m_serializer.set_threshold((*p)("vpeps_KdivU").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 0")
            m_serializer.set_threshold((*p)("vpeps_comp0").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 1")
            m_serializer.set_threshold((*p)("vpeps_comp1").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 2")
            m_serializer.set_threshold((*p)("vpeps_comp2").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 3")
            m_serializer.set_threshold((*p)("vpeps_comp3").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 4")
            m_serializer.set_threshold((*p)("vpeps_comp4").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 5")
            m_serializer.set_threshold((*p)("vpeps_comp5").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 6")
            m_serializer.set_threshold((*p)("vpeps_comp6").asDouble(vpeps));
        else if (TStreamer::NAME == "Component 7")
            m_serializer.set_threshold((*p)("vpeps_comp7").asDouble(vpeps));
        else
            m_serializer.set_threshold(vpeps);
    }

    virtual void write(const TGrid& grid, const int step_id, const Real t, const std::string basename, const std::string path)
    {
        m_serializer.Write(grid, basename, path);
    }

private:
    TCompressor<TGrid,TStreamer> m_serializer;
};


typedef BlockLabMPI<Lab> LabMPI;
typedef GridMPI<Grid_t> GridMPI_t;
typedef SliceTypesMPI::Slice<GridMPI_t> SliceGridMPI_t;
typedef SubdomainTypesMPI::Subdomain<GridMPI_t> SubdomainGridMPI_t;

typedef void (*Dumper_h5_full_grid_mpi)(const GridMPI_t&, const int, const Real, const std::string&, const std::string&, const bool);
typedef void (*Dumper_h5_slice_grid_mpi)(const SliceGridMPI_t&, const int, const Real, const std::string&, const std::string&, const bool);
typedef void (*Dumper_h5_subdomain_grid_mpi)(const SubdomainGridMPI_t&, const int, const Real, const std::string&, const std::string&, const bool);

// specialization
template <>
void ProcessingElement<GridMPI_t, Dumper_h5_full_grid_mpi, GridMPI_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_outputFunctor != NULL);
    m_outputFunctor(*m_grid, step_id, t, basename+m_streamerTag, path, true);
}

template <>
void ProcessingElement<EntityProcessor<SliceGridMPI_t>, Dumper_h5_slice_grid_mpi, GridMPI_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_outputFunctor != NULL);
    m_entity->setFunctor(m_outputFunctor);
    m_entity->process_all(step_id, t, basename+m_streamerTag, path);
}

template <>
void ProcessingElement<EntityProcessor<SubdomainGridMPI_t>, Dumper_h5_subdomain_grid_mpi, GridMPI_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_outputFunctor != NULL);
    m_entity->setFunctor(m_outputFunctor);
    m_entity->process_all(step_id, t, basename+m_streamerTag, path);
}

template <>
void ProcessingElement<SerializerWaveletBase<GridMPI_t>, void*, GridMPI_t>::operator()(const int step_id, const Real t, const std::string basename, const std::string& path, const bool verbose)
{
    assert(m_entity != NULL);
    m_entity->write(*m_grid, step_id, t, basename+m_streamerTag, path);
}

template <>
void ProcessingElement<SerializerWaveletBase<GridMPI_t>, void*, GridMPI_t>::dispose()
{
    if (m_entity)
        delete m_entity;
}


template <typename TGrid, template <typename> class TSlice=SliceTypesMPI::Slice, template <typename> class TSubdomain=SubdomainTypesMPI::Subdomain>
class OutputProcessingMPI : public OutputProcessing<TGrid,TSlice,TSubdomain>
{
public:
    OutputProcessingMPI(ArgumentParser& p, TGrid& grid, const bool verbose=true) :
        OutputProcessing<TGrid,TSlice,TSubdomain>(p,grid,verbose)
    {
        MPI_Comm comm = grid.getCartComm();
        MPI_Comm_rank(comm, &m_myrank);
    }

    void show_processor_entries() const override
    {
        std::ostringstream prefix;
        prefix << "[rank=" << m_myrank << "] ";
        this->m_sliceProcessor.showEntities(prefix.str());
        this->m_subdomainProcessor.showEntities(prefix.str());
    }


protected:
    typedef typename OutputProcessing<TGrid,TSlice,TSubdomain>::TFunc_h5_grid  TFunc_h5_gridMPI;
    typedef typename OutputProcessing<TGrid,TSlice,TSubdomain>::TFunc_h5_slice TFunc_h5_sliceMPI;
    typedef typename OutputProcessing<TGrid,TSlice,TSubdomain>::TFunc_h5_subdomain TFunc_h5_subdomainMPI;
    typedef typename OutputProcessing<TGrid,TSlice,TSubdomain>::TProcessorSlice TProcessorSlice;
    typedef typename OutputProcessing<TGrid,TSlice,TSubdomain>::TProcessorSubdomain TProcessorSubdomain;
    typedef typename OutputProcessing<TGrid,TSlice,TSubdomain>::USlice USlice;
    typedef typename OutputProcessing<TGrid,TSlice,TSubdomain>::USubdomain USubdomain;
    typedef SerializerWaveletBase<TGrid> TProcessorVP;

    virtual void _register_h5_grid(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5, HDF5_FullMPI, PElement::H5, true, TGrid, TFunc_h5_gridMPI, TGrid, grid, grid, (this->m_parser), DumpHDF5_MPI, TGrid, process, LabMPI)
    }

    virtual void _register_h5_slice(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5s, HDF5_SliceMPI, PElement::H5SLICE, false, TProcessorSlice, TFunc_h5_sliceMPI, TGrid, (this->m_sliceProcessor), grid, (this->m_parser), DumpSliceHDF5MPI, USlice, process, LabMPI)
    }

    virtual void _register_h5_subdomain(TGrid& grid)
    {
        __REGISTER_H5_ENTITY__(h5sub, HDF5_SubdomainMPI, PElement::H5SUBDOMAIN, true, TProcessorSubdomain, TFunc_h5_subdomainMPI, TGrid, (this->m_subdomainProcessor), grid, (this->m_parser), DumpSubdomainHDF5MPI, USubdomain, process, LabMPI)
    }

    virtual void _register_vp(TGrid& grid)
    {
        __REGISTER_VP_ENTITY__(vp, VP_Wavelet, PElement::VP, true, TProcessorVP, TGrid, grid, (this->m_parser), process, LabMPI)
    }

    virtual void _register_all(TGrid& grid)
    {
        _register_h5_grid(grid);
        _register_h5_slice(grid);
        _register_h5_subdomain(grid);
        _register_vp(grid);
    }

private:
    int m_myrank;
};

#endif /* OUTPUTPROCESSINGMPI_H_TQSCO9A3 */
