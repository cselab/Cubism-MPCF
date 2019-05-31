/*
 *  Test_EllipsoidDrop.h
 *  MPCFnode
 *
 *  Created by Fabian Wermelinger on 8/31/17.
 *  Copyright 2017 ETH Zurich. All rights reserved.
 *
 */
#ifndef TEST_ELLIPSOIDDROP_H_GAFJNIX5
#define TEST_ELLIPSOIDDROP_H_GAFJNIX5

#include <omp.h>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <limits>

#include "Test_SteadyState.h"


template <typename TGrid, typename TStepper, template <typename> class TSlice=SliceTypes::Slice, template <typename> class TSubdomain=SubdomainTypes::Subdomain>
class Test_EllipsoidDrop: public virtual Test_SteadyState<TGrid,TStepper,TSlice,TSubdomain>
{
public:
    Test_EllipsoidDrop(ArgumentParser& parser) :
        Test_SteadyState<TGrid,TStepper,TSlice,TSubdomain>(parser) { }
    virtual ~Test_EllipsoidDrop() { }


protected:
    typedef typename TGrid::BlockType TBlock;
    typedef typename TBlock::ElementType TElement;

    Real m_R0;
    Real m_alpha;
    Real m_ia, m_ib, m_ic;

    virtual void _print_case_header()
    {
        printf("////////////////////////////////////////////////////////////\n");
        printf("////////////   TEST ELLIPSOID DROP/BUBBLE    ///////////////\n");
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
        this->_parameter_ellipsoid_drop();
    }

    virtual void _init()
    {
        Test_SteadyState<TGrid,TStepper,TSlice,TSubdomain>::_init();
    }

    virtual void _ic();

    virtual void _analysis()
    {
        if (this->ANALYSISPERIOD != 0 && this->step_id%this->ANALYSISPERIOD == 0 && !this->bRESTART)
        {
            this->profiler.push_start("ANALYSIS");
            this->_ellipsoid_drop_simple_analysis();
            this->profiler.pop_stop();
        }
    }

    void _parameter_ellipsoid_drop();
    void _ellipsoid_drop_simple_analysis();
    TElement _get_IC(const double p[3]);
    TElement _integral(const double p[3], const double h);
};


///////////////////////////////////////////////////////////////////////////////
// Class Implementation
///////////////////////////////////////////////////////////////////////////////
template <typename TGrid, typename TStepper, template <typename> class TSlice, template <typename> class TSubdomain>
void Test_EllipsoidDrop<TGrid,TStepper,TSlice,TSubdomain>::_ic()
{
    std::vector<BlockInfo> vInfo = this->grid->getBlocksInfo();
#ifdef _NONUNIFORM_BLOCK_
    const double h = this->m_nonuniform->minimum_cell_width();
#else
    const double h = vInfo[0].h_gridpoint;
#endif /* _NONUNIFORM_BLOCK_ */

#pragma omp parallel
    {
#ifdef _USE_NUMA_
        const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
        const int mynode = omp_get_thread_num() / cores_per_node;
        numa_run_on_node(mynode);
#endif

#pragma omp for
        for(int i=0; i<(int)vInfo.size(); i++)
        {
            BlockInfo info = vInfo[i];
            TBlock& b = *(TBlock*)info.ptrBlock;
            for(int iz=0; iz<TBlock::sizeZ; iz++)
                for(int iy=0; iy<TBlock::sizeY; iy++)
                    for(int ix=0; ix<TBlock::sizeX; ix++)
                    {
                        double pos[3];
                        info.pos(pos, ix, iy, iz);
                        b(ix,iy,iz) = _integral(pos, 0.5*h);
                    }
        }
    }
}


template <typename TGrid, typename TStepper, template <typename> class TSlice, template <typename> class TSubdomain>
void Test_EllipsoidDrop<TGrid,TStepper,TSlice,TSubdomain>::_parameter_ellipsoid_drop()
{
    // // mixture epsilon
    // m_epsilon = this->parser("-epsilon").asDouble(0.0);

    // initial bubble radius
    m_R0 = this->parser("R0").asDouble(1.0);

    // initial stretching of m_R0 along x-axis
    m_alpha = this->parser("alpha").asDouble(1.0);

    // ellipsoid axes (assumes a = b), inverse values 1/a^2, 1/b^2, 1/c^2:
    m_ia = std::pow(1.0 / (m_alpha * m_R0), 2);
    m_ib = m_ia;
    m_ic = std::pow((m_alpha*m_alpha) / m_R0, 2);
}


template <typename TGrid, typename TStepper, template <typename> class TSlice, template <typename> class TSubdomain>
void Test_EllipsoidDrop<TGrid,TStepper,TSlice,TSubdomain>::_ellipsoid_drop_simple_analysis()
{
}


template <typename TGrid, typename TStepper, template <typename> class TSlice, template <typename> class TSubdomain>
typename Test_EllipsoidDrop<TGrid,TStepper,TSlice,TSubdomain>::TElement Test_EllipsoidDrop<TGrid,TStepper,TSlice,TSubdomain>::_get_IC(const double p[3])
{
    const Real x = p[0]-0.5*Simulation_Environment::extents[0];
    const Real y = p[1]-0.5*Simulation_Environment::extents[1];
    const Real z = p[2]-0.5*Simulation_Environment::extents[2];
    const Real d = (x*x*m_ia + y*y*m_ib + z*z*m_ic - 1.0) * m_R0;
    const double eps = Simulation_Environment::EPSILON;
    const double alpha = M_PI*std::min(1., std::max(0., (double)(d + eps)/(2 * eps)));
    const double a2 = 0.5 + 0.5 * std::cos(alpha); // 1 = inside ellipsoid
    const double a1 = 1.0 - a2;

    // pressure (assume linear variation across interface)
    const double pressure = Simulation_Environment::P1 + a2*Simulation_Environment::SIGMA*2.0/m_R0; // second term Laplace pressure

    // material mixture
    const double G1 = 1.0/(Simulation_Environment::GAMMA1-1.0);
    const double G2 = 1.0/(Simulation_Environment::GAMMA2-1.0);
    const double F1 = Simulation_Environment::GAMMA1*Simulation_Environment::PC1;
    const double F2 = Simulation_Environment::GAMMA2*Simulation_Environment::PC2;
    const double Gmix = a1*G1 + a2*G2;
    const double Pmix = a1*F1*G1 + a2*F2*G2;
    const double Emix = pressure*Gmix + Pmix; // Ekin = 0

    TElement el;
    el.alpha1rho1 = a1*Simulation_Environment::RHO1;
    el.alpha2rho2 = a2*Simulation_Environment::RHO2;
    el.ru = 0.0;
    el.rv = 0.0;
    el.rw = 0.0;
    el.energy = Emix;
    el.alpha2 = a2;
    return el;
}


template <typename TGrid, typename TStepper, template <typename> class TSlice, template <typename> class TSubdomain>
typename Test_EllipsoidDrop<TGrid,TStepper,TSlice,TSubdomain>::TElement Test_EllipsoidDrop<TGrid,TStepper,TSlice,TSubdomain>::_integral(const double p[3], const double h) // h should be half the grid spacing h
{
    TElement samples[3][3][3];
    TElement zintegrals[3][3];
    TElement yzintegrals[3];

    const double x0[3] = {p[0] - h, p[1] - h, p[2] - h};

    for(int iz=0; iz<3; iz++)
        for(int iy=0; iy<3; iy++)
            for(int ix=0; ix<3; ix++)
            {
                const double mypos[3] = {x0[0]+ix*h, x0[1]+iy*h, x0[2]+iz*h};
                samples[iz][iy][ix] = _get_IC(mypos);
            }

    for(int iy=0; iy<3; iy++)
        for(int ix=0; ix<3; ix++)
            zintegrals[iy][ix] = (1/6.) * samples[0][iy][ix]+(2./3) * samples[1][iy][ix]+(1/6.)* samples[2][iy][ix];

    for(int ix=0; ix<3; ix++)
        yzintegrals[ix] = (1./6)*zintegrals[0][ix] + (2./3)*zintegrals[1][ix]+(1./6)*zintegrals[2][ix];

    return (1./6) * yzintegrals[0]+(2./3) * yzintegrals[1]+(1./6)* yzintegrals[2];
}

#endif /* TEST_ELLIPSOIDDROP_H_GAFJNIX5 */
