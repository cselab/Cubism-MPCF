//
//  PoissonSolver_MPI.h
//  MPCFcluster
//
//  Created by Babak Hejazialhosseini on 2/10/12.
//  Copyright 2011 ETH Zurich. All rights reserved.
//
#pragma once

#include <GridMPI.h>
#include <BlockInfo.h>
#include <Types.h>

#ifdef _USE_FFTW_
#include "PoissonSolverPeriodic.h"
#endif
#include <BlockLabMPI.h>
#include <BlockProcessingMPI.h>

typedef GridMPI< FluidGrid > GridType;
typedef BlockLabMPI< BlockLab<FluidBlock, tbb::scalable_allocator > > LabMPI_PossionSolver;

using namespace LSRK3MPIdata;

struct CurlOp
{
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];
    
    CurlOp(): stencil(-1,-1,-1,2,2,2, false, 4, 0,1,2,3)
    {
        stencil_start[0] = -1;
        stencil_start[1] = -1;
        stencil_start[2] = -1;
        
        stencil_end[0] = +2;
        stencil_end[1] = +2;
        stencil_end[2] = +2;
    }
    
    CurlOp(const CurlOp& c): stencil(-1,-1,-1,2,2,2, false, 4, 0,1,2,3)
    {
        stencil_start[0] = -1;
        stencil_start[1] = -1;
        stencil_start[2] = -1;
        
        stencil_end[0] = +2;
        stencil_end[1] = +2;
        stencil_end[2] = +2;
    }
    
	template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& b) const
    {
        typedef BlockType B;
        typedef typename BlockType::ElementType E;
        const Real hInvHalf = (Real)(0.5/info.h_gridpoint);
        
		for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    b.tmp[iz][iy][ix][1]= hInvHalf*lab(ix,iy,iz).rho*(lab(ix,iy+1,iz).w-lab(ix,iy-1,iz).w - (lab(ix,iy,iz+1).v-lab(ix,iy,iz-1).v));
                    b.tmp[iz][iy][ix][2]= hInvHalf*lab(ix,iy,iz).rho*(lab(ix,iy,iz+1).u-lab(ix,iy,iz-1).u - (lab(ix+1,iy,iz).w-lab(ix-1,iy,iz).w));
                    b.tmp[iz][iy][ix][3]= hInvHalf*lab(ix,iy,iz).rho*(lab(ix+1,iy,iz).v-lab(ix-1,iy,iz).v - (lab(ix,iy+1,iz).u-lab(ix,iy-1,iz).u));
                    b.tmp[iz][iy][ix][4]= 0.5*(pow(b(ix,iy,iz).u,2) + pow(b(ix,iy,iz).v,2) + pow(b(ix,iy,iz).w,2))/b(ix,iy,iz).rho;
                    b.tmp[iz][iy][ix][5]= b(ix,iy,iz).levelset;
                }
    }
};  

template <typename BlockType>
struct tmpToSV
{
    BlockInfo * ary;
    
public:  
    tmpToSV(BlockInfo * ary): ary(ary) { }
    
    void operator()(blocked_range<int> range) const
    {
        for(int r=range.begin(); r<range.end(); ++r)
        {
            BlockType & b = *(BlockType *)ary[r].ptrBlock;
            
            for(int iz=0; iz<BlockType::sizeZ; iz++)
                for(int iy=0; iy<BlockType::sizeY; iy++)
                    for(int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        b(ix,iy,iz).u = b.tmp[iz][iy][ix][1];
                        b(ix,iy,iz).v = b.tmp[iz][iy][ix][2];
                        b(ix,iy,iz).w = b.tmp[iz][iy][ix][3];
                        b(ix,iy,iz).energy = b.tmp[iz][iy][ix][4];
                        b(ix,iy,iz).levelset = b.tmp[iz][iy][ix][5];
                    }
        }
    }
};

struct Velocity_ContractionOp
{
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];
    
    Velocity_ContractionOp(): stencil(-1,-1,-1,2,2,2, false, 4, 0,1,2,3)
    {
        stencil_start[0] = -1;
        stencil_start[1] = -1;
        stencil_start[2] = -1;
        
        stencil_end[0] = +2;
        stencil_end[1] = +2;
        stencil_end[2] = +2;
    }
    
    Velocity_ContractionOp(const Velocity_ContractionOp& c): stencil(-1,-1,-1,2,2,2, false, 4, 0,1,2,3)
	{
        stencil_start[0] = -1;
        stencil_start[1] = -1;
        stencil_start[2] = -1;
        
        stencil_end[0] = +2;
        stencil_end[1] = +2;
        stencil_end[2] = +2;
	}
    
    template<typename LabType, typename BlockType>
    inline void operator()(LabType& lab, const BlockInfo& info, BlockType& b) const
    {
        typedef BlockType B;
        typedef typename BlockType::ElementType E;
        const Real hInv2 = pow((Real)(0.5/info.h_gridpoint),2);

		for(int iz=0; iz<B::sizeZ; iz++)
            for(int iy=0; iy<B::sizeY; iy++)
                for(int ix=0; ix<B::sizeX; ix++)
                {
                    const Real dudx = lab(ix+1,iy,iz).u/lab(ix+1,iy,iz).rho - lab(ix-1,iy,iz).u/lab(ix-1,iy,iz).rho;
                    const Real dudy = lab(ix,iy+1,iz).u/lab(ix,iy+1,iz).rho - lab(ix,iy-1,iz).u/lab(ix,iy-1,iz).rho;
                    const Real dudz = lab(ix,iy,iz+1).u/lab(ix,iy,iz+1).rho - lab(ix,iy,iz-1).u/lab(ix,iy,iz-1).rho;
                    const Real dvdx = lab(ix+1,iy,iz).v/lab(ix+1,iy,iz).rho - lab(ix-1,iy,iz).v/lab(ix-1,iy,iz).rho;
                    const Real dvdy = lab(ix,iy+1,iz).v/lab(ix,iy+1,iz).rho - lab(ix,iy-1,iz).v/lab(ix,iy-1,iz).rho;
                    const Real dvdz = lab(ix,iy,iz+1).v/lab(ix,iy,iz+1).rho - lab(ix,iy,iz-1).v/lab(ix,iy,iz-1).rho;
                    const Real dwdx = lab(ix+1,iy,iz).w/lab(ix+1,iy,iz).rho - lab(ix-1,iy,iz).w/lab(ix-1,iy,iz).rho;
                    const Real dwdy = lab(ix,iy+1,iz).w/lab(ix,iy+1,iz).rho - lab(ix,iy-1,iz).w/lab(ix,iy-1,iz).rho;
                    const Real dwdz = lab(ix,iy,iz+1).w/lab(ix,iy,iz+1).rho - lab(ix,iy,iz-1).w/lab(ix,iy,iz-1).rho;
                    
                    b.tmp[iz][iy][ix][0] = (dudx*dudx+dvdy*dvdy+dwdz*dwdz + 2*(dudy*dvdx+dudz*dwdx+dvdz*dwdy))*hInv2;
                }
    }
};

template <typename TGrid>
struct CurlPsi
{
    TGrid & grid;
    
    CurlPsi(TGrid & _grid): grid(_grid) {}
    
    void execute()
    {            
        CurlOp curl;        
        SynchronizerMPI& synch = ((TGrid&)grid).sync(curl);        
        while (!synch.done())
        {                
            vector<BlockInfo> avail = synch.avail(2);
            BlockProcessingMPI::process< LabMPI_PossionSolver >(avail, curl, (TGrid&)grid);
        }
    }
};

template <typename TGrid>
struct VelocityContraction
{
    TGrid & grid;
    
    VelocityContraction(TGrid & _grid): grid(_grid) {}
    
    void execute()
    {        
        Velocity_ContractionOp vel_cont;
        SynchronizerMPI& synch = ((TGrid&)grid).sync(vel_cont);
        while (!synch.done())
        {                
            vector<BlockInfo> avail = synch.avail(2);
            BlockProcessingMPI::process< LabMPI_PossionSolver >(avail, vel_cont, (TGrid&)grid);
        }
    }
};

template <typename TGrid>
class PoissonSolver_MPI
{
    TGrid& grid;
    Real gamma;
    
public: 
    
    PoissonSolver_MPI(TGrid& _grid, const Real _gamma): grid(_grid), gamma(_gamma) { }
    
    void solve_poisson()
    {             
#ifdef _USE_FFTW_
		CurlPsi< TGrid > curlpsi(grid);
        VelocityContraction< TGrid > vel_cont(grid);
        vector<BlockInfo> vInfo = grid.getResidentBlocksInfo();
		
		PeriodicPoisson3D_FFTW_MPI<TGrid,FluidBlock> poisson0(&grid, 0, 1, &vInfo.front());
        parallel_for(blocked_range<int>(0, vInfo.size()), poisson0, auto_partitioner());
        poisson0.Fourier2Space();
        MPI_Barrier(MPI_COMM_WORLD);
		
        poisson0.remove_average();
        poisson0.LaplacianInFourier_4th();
        parallel_for(blocked_range<int>(0, vInfo.size()), poisson0, auto_partitioner());
        poisson0.Space2Fourier();
		MPI_Barrier(MPI_COMM_WORLD);
		poisson0.Free();
        MPI_Barrier(MPI_COMM_WORLD);
		
        PeriodicPoisson3D_FFTW_MPI<TGrid,FluidBlock> poisson1(&grid, 1, 1, &vInfo.front());
        parallel_for(blocked_range<int>(0, vInfo.size()), poisson1, auto_partitioner());
        poisson1.Fourier2Space();
        MPI_Barrier(MPI_COMM_WORLD);
        
        poisson1.remove_average();
        poisson1.LaplacianInFourier_4th();
        parallel_for(blocked_range<int>(0, vInfo.size()), poisson1, auto_partitioner());
        poisson1.Space2Fourier();
        MPI_Barrier(MPI_COMM_WORLD);
		poisson1.Free();
        MPI_Barrier(MPI_COMM_WORLD);
		
        PeriodicPoisson3D_FFTW_MPI<TGrid,FluidBlock> poisson2(&grid, 2, 1, &vInfo.front());
        parallel_for(blocked_range<int>(0, vInfo.size()), poisson2, auto_partitioner());
        poisson2.Fourier2Space();
        MPI_Barrier(MPI_COMM_WORLD);
        
        poisson2.remove_average();
        poisson2.LaplacianInFourier_4th();	
        parallel_for(blocked_range<int>(0, vInfo.size()), poisson2, auto_partitioner());
        poisson2.Space2Fourier();
        MPI_Barrier(MPI_COMM_WORLD);
		poisson2.Free();
        MPI_Barrier(MPI_COMM_WORLD);
		
		curlpsi.execute();
        
        MPI_Barrier(MPI_COMM_WORLD);
        tmpToSV<FluidBlock> tmp2sv(&vInfo.front());
        MPI_Barrier(MPI_COMM_WORLD);
        parallel_for(blocked_range<int>(0, vInfo.size()), tmp2sv, auto_partitioner());
        
        MPI_Barrier(MPI_COMM_WORLD);
        vel_cont.execute();
        MPI_Barrier(MPI_COMM_WORLD);
		
		PeriodicPoisson3D_FFTW_MPI<TGrid,FluidBlock> poisson_pressure(&grid, 0, 0, &vInfo.front(), gamma);
        parallel_for(blocked_range<int>(0, vInfo.size()), poisson_pressure, auto_partitioner()); 
        poisson_pressure.Fourier2Space();
        MPI_Barrier(MPI_COMM_WORLD);
        
        poisson_pressure.remove_average();
        poisson_pressure.LaplacianInFourier_4th();
        parallel_for(blocked_range<int>(0, vInfo.size()), poisson_pressure, auto_partitioner()); 
        poisson_pressure.Space2Fourier();
        MPI_Barrier(MPI_COMM_WORLD);
		poisson_pressure.Free();
#else
#warning USE OF FFTW WAS DISABLED AT COMPILE TIME
#endif
    }
};
