/*
 *  Test_Cloud.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 2/24/13.
 *  Completly revised by Ursula Rasthofer in 2016
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_CLOUD_H_XEYR5OG6
#define TEST_CLOUD_H_XEYR5OG6

#include <cassert>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "Test_SteadyState.h"
#include "Shapes.h"

struct CloudICData
{
    double R0;
    double extent[3];
    double epsilon;
    double iface_thickness;
    double prox;
    double buff;
    bool integral;
    bool laplacep;
    std::string laplacep_file;
};


template <typename TGrid, typename TStepper, template <typename> class TSlice=SliceTypes::Slice, template <typename> class TSubdomain=SubdomainTypes::Subdomain>
class Test_Cloud: public virtual Test_SteadyState<TGrid,TStepper,TSlice,TSubdomain>
{
public:
    Test_Cloud(ArgumentParser& P) :
        Test_SteadyState<TGrid,TStepper,TSlice,TSubdomain>(P)
    { }
    virtual ~Test_Cloud() { }


protected:
    typedef typename TGrid::BlockType TBlock;
    typedef typename TBlock::ElementType TElement;

    CloudICData m_clData;

    // case specific
    virtual void _post_save()
    {
        if (this->isroot)
            std::cout << "ODE BC serialization currently not supported on node layer" << std::endl;
    }
    virtual void _post_restart()
    {
        if (this->isroot)
            std::cout << "ODE BC serialization currently not supported on node layer" << std::endl;
    }

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////             TEST CLOUD          ///////////////\n");
        printf("////////////////////////////////////////////////////////////\n");
        typedef typename TGrid::BlockType B;
        std::cout << "Domain size:   [" << this->BPDX*B::sizeX;
        std::cout << " x " << this->BPDY*B::sizeY;
        std::cout << " x " <<  this->BPDZ*B::sizeZ << "]" << std::endl;

        std::cout << "Domain extent: [" << Simulation_Environment::extents[0];
        std::cout << " x " << Simulation_Environment::extents[1];
        std::cout << " x " <<  Simulation_Environment::extents[2] << "]" << std::endl;
    }

    virtual void _setup_parameter()
    {
        Test_SteadyState<TGrid,TStepper,TSlice,TSubdomain>::_setup_parameter();
        this->_parameter_cloud();
    }

    virtual void _ic();

    virtual Seed<shape> _make_shapes()
    {
        std::string b_file = this->parser("-cloud-dat").asString("cloud.dat");
        Seed<shape> myseed;
        myseed.make_shapes(b_file, this->isroot);
        return myseed;
    }

    virtual void _load_pressure_file()
    {
        ReadHDF5<StreamerDummy,DumpReal>(*this->grid, this->m_clData.laplacep_file);
    }

    void _parameter_cloud();
};

///////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////
// small helpers
static inline double _get_energy(const double pressure, const double alpha2)
{
    // REMARK: Simulation_Environment::GAMMA1, ... are stored as
    // Real, but if the are not very fancy values, this is ok
// XXX: [fabianw@mavt.ethz.ch; Tue Aug 28 2018 11:49:38 AM (+0200)] WARNING:
// this assumes kinetic energy is zero; Ke = 0.
    const double g1 = Simulation_Environment::GAMMA1;
    const double g2 = Simulation_Environment::GAMMA2;
    const double pc1 = Simulation_Environment::PC1;
    const double pc2 = Simulation_Environment::PC2;
    const double g1m1Inv = 1.0/(g1-1.0);
    const double g2m1Inv = 1.0/(g2-1.0);

    const double alpha1 = 1.0-alpha2;
    const double gmix_m1Inv = alpha1*g1m1Inv + alpha2*g2m1Inv;
    const double pcmix = alpha1*g1*pc1*g1m1Inv + alpha2*g2*pc2*g2m1Inv;
    const double energy = gmix_m1Inv*pressure + pcmix;

    return energy;
}

// Variations of alpha2 and pressure initialization functions
///////////////////////////////////////////////////////////////////////////////
static inline double _alpha_std(const double p[3], const double alpha2, const CloudICData& data)
{
    return alpha2;
}

static inline double _alpha_tiwari(const double p[3], const double alpha2, const CloudICData& data)
{
#if defined(_BCLABCLOUDSYMABSORB_) || defined(_BCLABCLOUDSYM1DCHARNREF_)
    const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
#else
    const double r = std::sqrt((p[0] - 0.5*data.extent[0])*(p[0] - 0.5*data.extent[0]) + (p[1] - 0.5*data.extent[1])*(p[1] - 0.5*data.extent[1]) + (p[2] - 0.5*data.extent[2])*(p[2] - 0.5*data.extent[2]));
#endif
    return 0.5*(1.0 - std::tanh(0.5 * (r - data.R0)/(data.iface_thickness)));
}

static inline double _pressure_tiwari(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    // initial conditions according to Tiwari et al. 2013
    // careful: this option overwrites all fields
#if defined(_BCLABCLOUDSYMABSORB_) || defined(_BCLABCLOUDSYM1DCHARNREF_)
    const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
#else
    const double r = std::sqrt((p[0] - 0.5*data.extent[0])*(p[0] - 0.5*data.extent[0]) + (p[1] - 0.5*data.extent[1])*(p[1] - 0.5*data.extent[1]) + (p[2] - 0.5*data.extent[2])*(p[2] - 0.5*data.extent[2]));
#endif

    if (r <= data.R0)
        return CloudData::p2;
    else
        return CloudData::p1 - (data.R0/r) * (CloudData::p1 -  CloudData::p2);
}

static inline double _pressure_jump(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    // should not be used: pressure field with jump at the interface
    // causes strong oscillations in the pressure field
    return CloudData::p1*(1.0-alpha2) + CloudData::p2*alpha2;
}

template <int _COS>
static inline double _pressure_eqCloud(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    // entire cloud at constant pressure
    // increase of pressure towards the domain boundary
#if defined(_BCLABCLOUDSYMABSORB_) || defined(_BCLABCLOUDSYM1DCHARNREF_)
    const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
#else
    const double r = std::sqrt((p[0] - 0.5*data.extent[0])*(p[0] - 0.5*data.extent[0]) + (p[1] - 0.5*data.extent[1])*(p[1] - 0.5*data.extent[1]) + (p[2] - 0.5*data.extent[2])*(p[2] - 0.5*data.extent[2]));
#endif
    const double dist = r - data.R0;

    if (dist <= data.buff)
        return CloudData::p2;
    else
    {
        if (_COS)
        {
            const double alpha = M_PI*std::min(1., std::max(0., (data.prox - (dist-data.buff))/data.prox));
            return CloudData::p2 + (0.5 + 0.5 * std::cos(alpha)) * (CloudData::p1 -  CloudData::p2);
        }
        else
            return CloudData::p2 + std::tanh(2.0*dist/data.prox) * (CloudData::p1 -  CloudData::p2);
    }
}

template <int _COS>
static inline double _pressure_tiwariCloud(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    // following ideas presented in Tiwari et al 2015
    const double dist = distance(v_shapes, p);

    if (dist <= data.buff)
        return CloudData::p2;
    else
    {
        if(_COS)
        {
            const double alpha = M_PI*std::min(1., std::max(0., (data.prox - (dist-data.buff))/data.prox));
            return CloudData::p2 + (0.5 + 0.5 * std::cos(alpha)) * (CloudData::p1 -  CloudData::p2);
        }
        else
            return CloudData::p2 + std::tanh(2.0*dist/data.prox) * (CloudData::p1 -  CloudData::p2);
    }
}

// explicit instantiation
template double _pressure_eqCloud<0>(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data); // tanh
template double _pressure_eqCloud<1>(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data); // cos
template double _pressure_tiwariCloud<0>(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data); // tanh
template double _pressure_tiwariCloud<1>(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data); // cos
///////////////////////////////////////////////////////////////////////////////

// global function handles for pressure and alpha
double (*pressure_handle)(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data);
double (*alpha_handle)(const double p[3], const double alpha2, const CloudICData& data);
///////////////////////////////////////////////////////////////////////////////


// ic for cloud: primitive variables
// pressure goes to energy and velocity to momentum
// (p is position in space)
template<typename T>
inline T get_ic_primitives(const double p[3], const std::vector< shape >& v_shapes, const double alpha2, const CloudICData& data)
{
    T out;
    out.clear();

    out.alpha2 = alpha_handle(p, alpha2, data);

    const double pres = pressure_handle(p, v_shapes, alpha2, data);
    out.energy = pres;
#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
    // store the pressure in the dummy as it is needed for characteristic-based boundary conditions
    // HINT: if dummy is need for something else we have to recompute the pressure from energy in _ic()
    out.dummy = pres;
#endif

    return out;
}

template<typename T>
inline T get_ic_conserved(const double p[3], const double h[3], const std::vector< shape >& v_shapes, const CloudICData& data)
{
    T out;
    out.clear();

    if (data.integral)
    {
      T samples[3][3][3];
      T zintegrals[3][3];
      T yzintegrals[3];

      const double x0[3] = {p[0]-0.5*h[0], p[1]-0.5*h[1], p[2]-0.5*h[2]};

      for(int iz=0; iz<3; iz++)
          for(int iy=0; iy<3; iy++)
              for(int ix=0; ix<3; ix++)
              {
                  double mypos[3] = {x0[0]+ix*0.5*h[0], x0[1]+iy*0.5*h[1], x0[2]+iz*0.5*h[2]};
                  for (int d = 0; d < 3; d++)
                    if (Simulation_Environment::BC_PERIODIC[d]) {
                      if (mypos [d] < 0.0) mypos [d] += static_cast<double>(Simulation_Environment::extents [d]);
                      if (mypos [d] > Simulation_Environment::extents[d]) mypos [d] -= static_cast<double>(Simulation_Environment::extents [d]);
                  }
                  if (v_shapes.size()>0)
                    samples[iz][iy][ix] = get_ic_primitives<T>(mypos, v_shapes, eval(v_shapes, mypos), data);
                  else
                    samples[iz][iy][ix] = get_ic_primitives<T>(mypos, v_shapes, 0, data);
              }

      for(int iy=0; iy<3; iy++)
          for(int ix=0; ix<3; ix++)
              zintegrals[iy][ix] = (1/6.) * samples[0][iy][ix]+(2./3) * samples[1][iy][ix]+(1/6.)* samples[2][iy][ix];

      for(int ix=0; ix<3; ix++)
          yzintegrals[ix] = (1./6)*zintegrals[0][ix] + (2./3)*zintegrals[1][ix]+(1./6)*zintegrals[2][ix];

      out = (1./6) * yzintegrals[0]+(2./3) * yzintegrals[1]+(1./6)* yzintegrals[2];
    }
    else
    {
      if (v_shapes.size()>0)
        out = get_ic_primitives<T>(p, v_shapes, eval(v_shapes, p), data);
      else
        out = get_ic_primitives<T>(p, v_shapes, 0, data);
    }

    // now we compute the conserved quantities
    out.alpha1rho1 = (1.0-out.alpha2)*CloudData::rho1;
    out.alpha2rho2 = out.alpha2*CloudData::rho2;
    // we assume zero velocities as also done below for the energy
//    const double rho_mix = out.alpha1rho1 + out.alpha2rho2;
//    out.ru = rho_mix*out.ru;
//    out.rv = rho_mix*out.rv;
//    out.rw = rho_mix*out.rw;

    if (!data.laplacep)
        out.energy = _get_energy(out.energy, out.alpha2);

   return out;
}

///////////////////////////////////////////////////////////////////////////////
// CLASS IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////
template <typename TGrid, typename TStepper, template <typename> class TSlice, template <typename> class TSubdomain>
void Test_Cloud<TGrid,TStepper,TSlice,TSubdomain>::_ic()
{
    if (this->isroot)
        std::cout << "Cloud Initial condition..." << std::endl;

    const double G1 = Simulation_Environment::GAMMA1-1;
    const double G2 = Simulation_Environment::GAMMA2-1;
    const double PC1 = Simulation_Environment::PC1;
    const double PC2 = Simulation_Environment::PC2;
    const double F1 = G1*PC1;
    const double F2 = G2*PC2;

    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();
#ifndef _NONUNIFORM_BLOCK_
    const double h[3] = {vInfo[0].h_gridpoint, vInfo[0].h_gridpoint, vInfo[0].h_gridpoint};
    m_clData.iface_thickness = m_clData.epsilon * h[0];
#else
    {
        // compute the smallest geometric mean for m_clData.iface_thickness.
        // This may not be correct for all cases, the assumption here is that
        // the bubbles are resolved in the _finest_ resolved region of the
        // non-uniform mesh.  The value computed here is only used for the
        // Tiwari-type initial conditions.  The iterations below are over the
        // whole computational domain, each rank will compute the same value.
        const unsigned int NX = this->grid->getMeshMap(0).ncells();
        const unsigned int NY = this->grid->getMeshMap(1).ncells();
        const unsigned int NZ = this->grid->getMeshMap(2).ncells();
        double h_min = HUGE_VAL;
        for (unsigned int iz = 0; iz < NZ; ++iz)
        {
            const double hz = this->grid->getMeshMap(2).cell_width(iz);

            for (unsigned int iy = 0; iy < NY; ++iy)
            {
                const double hy = this->grid->getMeshMap(1).cell_width(iy);

                for (unsigned int ix = 0; ix < NX; ++ix)
                {
                    const double hx = this->grid->getMeshMap(0).cell_width(ix);
                    const double h = std::pow(hx*hy*hz, 1.0/3.0);
                    h_min = (h<h_min) ? h : h_min;
                }
            }
        }

        m_clData.iface_thickness = m_clData.epsilon * h_min;
    }
#endif /* _NONUNIFORM_BLOCK_ */

    Seed<shape> myseed = this->_make_shapes();

    // rasthofer May 2016
    // set-up boundary ghost cells
#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
    BGC::BPDX = this->BPDX;
    BGC::BPDY = this->BPDY;
    BGC::BPDZ = this->BPDZ;

    const int xpesize = this->parser("-xpesize").asInt(1);
    const int ypesize = this->parser("-ypesize").asInt(1);
    const int zpesize = this->parser("-zpesize").asInt(1);

    const int NBX = this->BPDX * xpesize;
    const int NBY = this->BPDY * ypesize;
    const int NBZ = this->BPDZ * zpesize;

    BGC::pamb = CloudData::p1;
    BGC::L = this->parser("-L-boundary").asDouble(1.0);
    BGC::lambda = this->parser("-lambda").asDouble(0.75);

    // loop all block and check whether they are close to the boundary
    for(int i=0; i<(int)vInfo.size(); i++)
    {
        BlockInfo info = vInfo[i];

        BGC::BoundaryBlock dummy;
        if (info.index[0] == 0)
            BGC::bgblock_list_dir0_side0.push_back(dummy);
        if (info.index[0] == (NBX-1))
            BGC::bgblock_list_dir0_side1.push_back(dummy);
        if (info.index[1] == 0)
            BGC::bgblock_list_dir1_side0.push_back(dummy);
        if (info.index[1] == (NBY-1))
            BGC::bgblock_list_dir1_side1.push_back(dummy);
        if (info.index[2] == 0)
            BGC::bgblock_list_dir2_side0.push_back(dummy);
        if (info.index[2] == (NBZ-1))
            BGC::bgblock_list_dir2_side1.push_back(dummy);
    }

    // some debug output
    /*
       {
       cout << "\n-----------------------------------------------------------\n";
       cout << "Boundary blocks:\n" ;
       cout << "xdir lhs: " << BGC::bgblock_list_dir0_side0.size() << " blocks\n";
       cout << "xdir rhs: " << BGC::bgblock_list_dir0_side1.size() << " blocks\n";
       cout << "ydir lhs: " << BGC::bgblock_list_dir1_side0.size() << " blocks\n";
       cout << "ydir rhs: " << BGC::bgblock_list_dir1_side1.size() << " blocks\n";
       cout << "zdir lhs: " << BGC::bgblock_list_dir2_side0.size() << " blocks\n";
       cout << "zdir rhs: " << BGC::bgblock_list_dir2_side1.size() << " blocks\n";
       cout << "-----------------------------------------------------------" << endl;
       }
       */
#endif

    // if a Laplace solution is used for initial pressure, load it from
    // external FEM solution
    if (m_clData.laplacep)
        this->_load_pressure_file();

    // setting maxlayer-block=0 is fine for all ic, except for tiwari-cloud
    // in this case, each rank as to see its closest bubbles
    // recommended value for this case is 4, but you should carefully check
    // your ic before subitting a production run
    const int maxlayer = this->parser("maxlayer-block").asInt(0);

#pragma omp parallel
    {
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif

        {
#pragma omp for
            for(int i=0; i<(int)vInfo.size(); i++)
            {
                BlockInfo info = vInfo[i];
                double mystart[3] = {info.origin[0], info.origin[1], info.origin[2]};
                double myextent[3] = { info.block_extent[0], info.block_extent[1], info.block_extent[2] } ;

                TBlock& b = *(TBlock*)info.ptrBlock;

                int s[3] = {0,0,0};
                int e[3] = {TBlock::sizeX,TBlock::sizeY,TBlock::sizeZ};

#if defined(_CHARACTERISTIC_1D_BOUNDARY_)
                int delta_ix = 0;
                int delta_iy = 0;
                int delta_iz = 0;
                int bnx = 3;
                int bny = 3;

                if (info.index[0] == 0)
                    s[0] = -3;
                if (info.index[0] == (NBX-1))
                    e[0] = TBlock::sizeX+3;
                if (info.index[1] == 0)
                    s[1] = -3;
                if (info.index[1] == (NBY-1))
                    e[1] = TBlock::sizeY+3;
                if (info.index[2] == 0)
                    s[2] = -3;
                if (info.index[2] == (NBZ-1))
                    e[2] = TBlock::sizeZ+3;
#endif

#ifdef _KEEPALL_
                std::vector<shape> myshapes = myseed.get_shapes();
#else
                // April 20, 2016: rasthofer removed due to extended block bounding box below
                //vector<shape> myshapes = myseed.get_shapes().size() ? myseed.retain_shapes(info.origin, myextent).get_shapes() : myseed.get_shapes();
                //printf("processing %d out of %d. for this block i have: %d bubbles\n", i, (int)vInfo.size(), myshapes.size());

                std::vector<shape> myshapes;

                if (myseed.get_shapes().size() != 0)
                {
                    int bbox =1;
                    int layerplus = 0;
                    while(true)
                    {
                        double ms[3], me[3];
                        for (size_t rr=0; rr<3; rr++)
                        {
                            ms[rr] = mystart[rr] - bbox*info.block_extent[rr];
                            me[rr] = (1+2*bbox)*myextent[rr];
                        }

                        myshapes = myseed.retain_shapes(ms,me).get_shapes();

                        if (myshapes.size() == 0)
                            bbox++;
                        else
                        {
                            if (layerplus<maxlayer)
                            {
                                bbox++;
                                layerplus++;
                            }
                            else
                                break;
                        }
                    }
                }
                else
                    myshapes = myseed.get_shapes();
#endif

                {
                    for(int iz=s[2]; iz<e[2]; iz++)
                        for(int iy=s[1]; iy<e[1]; iy++)
                            for(int ix=s[0]; ix<e[0]; ix++)
                            {
                                double p[3];
                                info.pos(p, ix, iy, iz);
#ifdef _NONUNIFORM_BLOCK_
                                double h[3];
                                info.spacing(h, ix, iy, iz);
#endif /* _NONUNIFORM_BLOCK_ */

                                TElement myvalues;

                                myvalues  = get_ic_conserved<TElement>(p, h, myshapes, m_clData);

#if !defined(_CHARACTERISTIC_1D_BOUNDARY_)
                                if (m_clData.laplacep)
                                {
                                    b(ix,iy,iz).alpha1rho1 = myvalues.alpha1rho1;
                                    b(ix,iy,iz).alpha2rho2 = myvalues.alpha2rho2;
                                    b(ix,iy,iz).ru = myvalues.ru;
                                    b(ix,iy,iz).rv = myvalues.rv;
                                    b(ix,iy,iz).rw = myvalues.rw;
                                    b(ix,iy,iz).energy = _get_energy(b(ix,iy,iz).dummy, myvalues.alpha2);
                                    b(ix,iy,iz).alpha2 = myvalues.alpha2;
                                }
                                else
                                    b(ix,iy,iz) = myvalues;
#else
                                if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX) // regular cell in block
                                {
                                    if (m_clData.laplacep)
                                    {
                                        b(ix,iy,iz).alpha1rho1 = myvalues.alpha1rho1;
                                        b(ix,iy,iz).alpha2rho2 = myvalues.alpha2rho2;
                                        b(ix,iy,iz).ru = myvalues.ru;
                                        b(ix,iy,iz).rv = myvalues.rv;
                                        b(ix,iy,iz).rw = myvalues.rw;
                                        b(ix,iy,iz).energy = _get_energy(b(ix,iy,iz).dummy, myvalues.alpha2);
                                        b(ix,iy,iz).alpha2 = myvalues.alpha2;
                                    }
                                    else
                                        b(ix,iy,iz) = myvalues;
                                }
                                else
                                {
                                    BGC::BoundaryFluidElement mybgc;
                                    mybgc.rho = myvalues.alpha1rho1+myvalues.alpha2rho2;
                                    // IMPORTANT REMARK: if we need the dummy for something else
                                    // we can recompute the pressure from energy
                                    mybgc.pressure = myvalues.dummy;
                                    const double myrhoinv = 1.0/mybgc.rho;
                                    assert(!isnan(mybgc.rho));
                                    assert(!isnan(mybgc.pressure));
                                    assert(!isnan(myrhoinv));

                                    if (iz<0 && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY);
                                        delta_ix = 0;
                                        delta_iy = 0;
                                        delta_iz = 3;
                                        bnx = TBlock::sizeX;
                                        bny = TBlock::sizeY;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir2_side0[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir2_side0[bgb_index].fe[bgc_index].u=myvalues.rw*myrhoinv;
                                    }
                                    else if  (iz>=TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[1]%BGC::BPDY);
                                        delta_ix = 0;
                                        delta_iy = 0;
                                        delta_iz = -TBlock::sizeZ;
                                        bnx = TBlock::sizeX;
                                        bny = TBlock::sizeY;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir2_side1[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir2_side1[bgb_index].fe[bgc_index].u=myvalues.rw*myrhoinv;
                                    }
                                    else if (iz>=0 && iz<TBlock::sizeZ && iy<0 && ix>=0 && ix<TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ);
                                        delta_iy = 3;
                                        delta_ix = 0;
                                        delta_iz = 0;
                                        bnx = TBlock::sizeX;
                                        bny = 3;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir1_side0[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir1_side0[bgb_index].fe[bgc_index].u=myvalues.rv*myrhoinv;
                                    }
                                    else if (iz>=0 && iz<TBlock::sizeZ && iy>=TBlock::sizeY && ix>=0 && ix<TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[0]%BGC::BPDX+BGC::BPDX*(info.index[2]%BGC::BPDZ);
                                        delta_iy = -TBlock::sizeY;
                                        delta_ix = 0;
                                        delta_iz = 0;
                                        bnx = TBlock::sizeX;
                                        bny = 3;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir1_side1[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir1_side1[bgb_index].fe[bgc_index].u=myvalues.rv*myrhoinv;
                                    }
                                    else if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix<0)
                                    {
                                        const int bgb_index = info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ);
                                        delta_ix = 3;
                                        delta_iy = 0;
                                        delta_iz = 0;
                                        bnx = 3;
                                        bny = TBlock::sizeY;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir0_side0[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir0_side0[bgb_index].fe[bgc_index].u=myvalues.ru*myrhoinv;
                                    }
                                    else if (iz>=0 && iz<TBlock::sizeZ && iy>=0 && iy<TBlock::sizeY && ix>=TBlock::sizeX)
                                    {
                                        const int bgb_index = info.index[1]%BGC::BPDY+BGC::BPDY*(info.index[2]%BGC::BPDZ);
                                        delta_ix = -TBlock::sizeX;
                                        delta_iy = 0;
                                        delta_iz = 0;
                                        bny = TBlock::sizeY;
                                        bnx = 3;
                                        const int bgc_index = (ix+delta_ix) + bnx * (iy+delta_iy) + bnx*bny * (iz+delta_iz);

                                        BGC::bgblock_list_dir0_side1[bgb_index].fe[bgc_index]=mybgc;
                                        BGC::bgblock_list_dir0_side1[bgb_index].fe[bgc_index].u=myvalues.ru*myrhoinv;
                                    }
                                }
#endif
                            }
                }
            }
        }
    }

    if (this->isroot)
        cout << "done." << endl;
}

template <typename TGrid, typename TStepper, template <typename> class TSlice, template <typename> class TSubdomain>
void Test_Cloud<TGrid,TStepper,TSlice,TSubdomain>::_parameter_cloud()
{
    // this data is defined in TestLabs.h
    CloudData::rho1 = this->parser("rho1").asDouble(1000.0);
    CloudData::rho2 = this->parser("rho2").asDouble(1.0);
    CloudData::p1   = this->parser("p1").asDouble(100.0);
    CloudData::p2   = this->parser("p2").asDouble(0.0234);

    // pack relevant cloud data passed to pressure and alpha handles
    m_clData.R0      = this->parser("-R0").asDouble(1.0);
    m_clData.extent[0]  = Simulation_Environment::extents[0];
    m_clData.extent[1]  = Simulation_Environment::extents[1];
    m_clData.extent[2]  = Simulation_Environment::extents[2];
    m_clData.epsilon = this->parser("-eps-sharp").asDouble(0.75);
    if (m_clData.epsilon <= 0.0)
        m_clData.epsilon = 0.75;
// XXX: [fabianw@mavt.ethz.ch; Wed Oct 24 2018 02:57:56 PM (+0200)] Note:
// m_clData.iface_thickness depends on the geometry of the grid and is set in
// the _ic() method.  It is required for some initial condondition variations
// implemented above
    m_clData.iface_thickness = 1.0;
    m_clData.buff     = this->parser("-buffer").asDouble(0.0);
    m_clData.laplacep = false;
    m_clData.integral = this->parser("-ic-integral").asBool(true);

    // assign pressure and alpha handles
    const std::string ic_p = this->parser("-ic-pressure-cloud").asString("eqcloud");
    const bool cosine = this->parser("-cosine").asBool(true);
    pressure_handle = NULL;
    if (ic_p == "jump")
        pressure_handle = &_pressure_jump;
    else if (ic_p == "tiwari")
        pressure_handle = &_pressure_tiwari;
    else if (ic_p == "eqcloud")
    {
        if (cosine) pressure_handle = &_pressure_eqCloud<1>; // cos
        else        pressure_handle = &_pressure_eqCloud<0>; // tanh
        this->parser.set_strict_mode();
        m_clData.prox = this->parser("-proximity").asDouble();
        this->parser.unset_strict_mode();
    }
    else if (ic_p == "tiwari-cloud")
    {
        if (cosine) pressure_handle = &_pressure_tiwariCloud<1>; // cos
        else        pressure_handle = &_pressure_tiwariCloud<0>; // tanh
        this->parser.set_strict_mode();
        m_clData.prox = this->parser("-proximity").asDouble();
        this->parser.unset_strict_mode();
    }
    else if (ic_p == "laplacian")
    {
        m_clData.laplacep = true;
        pressure_handle = &_pressure_jump;
        this->parser.set_strict_mode();
        m_clData.laplacep_file = this->parser("-ic-pressure-file").asString();
        this->parser.unset_strict_mode();
    }
    assert(pressure_handle != NULL);

    alpha_handle = NULL;
    if (ic_p != "tiwari")
        alpha_handle = &_alpha_std;
    else
        alpha_handle = &_alpha_tiwari;
    assert(alpha_handle != NULL);
}

#endif /* TEST_CLOUD_H_XEYR5OG6 */
