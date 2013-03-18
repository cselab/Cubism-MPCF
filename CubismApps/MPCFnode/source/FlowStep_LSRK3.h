/*
 *  FlowStep_LSRK3.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <StencilInfo.h>

#ifdef _USE_NUMA_
#include <numa.h>
#include <omp.h>
#endif

#include "Types.h"

namespace LSRK3data 
{
	extern Real gamma1;
	extern Real gamma2 ;
	extern Real smoothlength ;
	extern int verbosity;
	extern float PEAKPERF_CORE, PEAKBAND;
	extern int NCORES, TLP;
    extern Real pc1;
    extern Real pc2 ;
	extern string dispatcher;
	extern int step_id;
	extern int ReportFreq;
	extern Real nu1, nu2;
    
	template < typename Kernel , typename Lab>
	struct FlowStep
	{
		StencilInfo stencil;
		
		Real a, dtinvh;
		
		int stencil_start[3];
		int stencil_end[3];
		
        FlowStep(Real a, Real dtinvh): a(a), dtinvh(dtinvh), stencil(-3,-3,-3,4,4,4, false, 7, 0,1,2,3,4,5,6)
		{
			stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;		
			stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
		}

		FlowStep(const FlowStep& c): a(c.a), dtinvh(c.dtinvh), stencil(-3,-3,-3,4,4,4, false, 7, 0,1,2,3,4,5,6)
		{
			stencil_start[0] = stencil_start[1] = stencil_start[2] = -3;		
			stencil_end[0] = stencil_end[1] = stencil_end[2] = 4;
		}
		
		inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock& o) const
		{
			Kernel kernel(a, dtinvh);
			
			const Real * const srcfirst = &lab(-3,-3,-3).rho;
			const int labSizeRow = lab.template getActualSize<0>();
			const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
			Real * const destfirst = &o.tmp[0][0][0][0];
			
			kernel.compute(srcfirst, FluidBlock::gptfloats, labSizeRow, labSizeSlice, 
						   destfirst, FluidBlock::gptfloats, FluidBlock::sizeX, FluidBlock::sizeX*FluidBlock::sizeY);
		}
	};
	
    template < typename Kernel, typename Lab >
    struct Diffusion
    {
        StencilInfo stencil;
        Real nu1, nu2, dtinvh;
        int stencil_start[3];
        int stencil_end[3];
        
        Diffusion(const Real _dtinvh, const Real _nu1=0, const Real _nu2=0): dtinvh(_dtinvh), nu1(_nu1), nu2(_nu2), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,5)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }
        
        Diffusion(const Diffusion& c): dtinvh(c.dtinvh), nu1(c.nu1), nu2(c.nu2), stencil(-1,-1,-1,2,2,2, true, 5, 0,1,2,3,5)
        {
            stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
            stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
        }
        
        inline void operator()(Lab& lab, const BlockInfo& info, FluidBlock& o) const
        {
            Kernel kernel(1, nu1, nu2, max((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), min((Real)1/(LSRK3data::gamma1-1), (Real)1/(LSRK3data::gamma2-1)), info.h_gridpoint, LSRK3data::smoothlength, dtinvh);
            
            const Real * const srcfirst = &lab(-1,-1,-1).rho;
            const int labSizeRow = lab.template getActualSize<0>();
            const int labSizeSlice = labSizeRow*lab.template getActualSize<1>();
            Real * const destfirst =  &o.tmp[0][0][0][0];
            kernel.compute(srcfirst, FluidBlock::gptfloats, labSizeRow, labSizeSlice,
                           destfirst, FluidBlock::gptfloats, FluidBlock::sizeX, FluidBlock::sizeX*FluidBlock::sizeY);
        }
    };
    
	template < typename Kernel >
	struct Update
	{
		Real b;
		BlockInfo * ary;
		
	public:
		
		Update(float b, BlockInfo * ary): b(b), ary(ary) { }
		Update(const Update& c): b(c.b), ary(c.ary) { } 
	    
		void omp(const int N)
		{
#pragma omp parallel
			{
#ifdef _USE_NUMA_
				const int cores_per_node = numa_num_configured_cpus() / numa_num_configured_nodes();
                const int mynode = omp_get_thread_num() / cores_per_node;
                numa_run_on_node(mynode);
#endif
                Kernel kernel(b);
                
#pragma omp for schedule(runtime)
                for(int r=0; r<N; ++r)
                {
                    FluidBlock & block = *(FluidBlock *)ary[r].ptrBlock;
                    kernel.compute(&block.tmp[0][0][0][0], &block.data[0][0][0].rho, block.gptfloats);
                }
			}
		}
	};
}

class FlowStep_LSRK3
{
private:
    FluidGrid& grid; 
    
protected:
    const Real gamma1, gamma2;
    Real smoothlength, h;
    Real current_time;
    Real PEAKPERF_CORE, PEAKBAND;
    
    string blockdispatcher;
    
    const int verbosity;
    Real pc1, pc2;
    
    bool bAwk;
    
    Real _computeSOS(bool bAwk=false);
    
    ArgumentParser parser;
    Profiler * profiler;
    
public:
    Real CFL;
    
    FlowStep_LSRK3(FluidGrid& grid, const Real CFL, const Real gamma1, const Real gamma2, ArgumentParser& parser, const int verbosity=1, Profiler* profiler=NULL, const Real pc1=0, const Real pc2=0, const bool bAwk=false):
    grid(grid), CFL(CFL), gamma1(gamma1), gamma2(gamma2), parser(parser), verbosity(verbosity), profiler(profiler), pc1(pc1), pc2(pc2), bAwk(bAwk) 
    {
        parser.unset_strict_mode();
        
        PEAKPERF_CORE = parser("-pp").asDouble(27.2*2);
        PEAKBAND = parser("-pb").asDouble(19);
        blockdispatcher = parser("-dispatcher").asString("");
        
        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        h = vInfo[0].h_gridpoint;
        smoothlength = (Real)(parser("-mollfactor").asInt())*sqrt(3.)*h;
        Simulation_Environment::EPSILON = smoothlength;
    }
    
    Real operator()(const Real max_dt);
    
    void set_current_time(const Real _current_time) {current_time=_current_time;}
    
    void set_constants();
    
    void set_CFL(const Real _CFL) {CFL = _CFL;}
};
