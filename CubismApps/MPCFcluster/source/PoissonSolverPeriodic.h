/*
 *  PoissonSolverPeriodic.h
 *  Cubism
 *
 *  Created by Babak Hejazialhosseini on 6/24/09.
 *  Copyright 2009 ETHZ. All rights reserved.
 *
 */
#pragma once

#include <assert.h>
#include "mpi.h"

#include "BlockInfo.h"

#define MPI_FP MPI_DOUBLE
#include "dfftw_mpi.h"
#include "dfftw.h"

template<typename GridType, typename BlockType>
struct PeriodicPoisson3D_FFTW_MPI
{
private:
	static const unsigned int dim = 3;
	size_t nx;
	size_t ny;	
	size_t nz;
	size_t nbpercpux;
	size_t nbpercpuy;	
	size_t nbpercpuz;
	int local_nx;
	int local_x_start;
	int local_ny_after_transpose;
	int local_y_start_after_transpose;
	int total_local_size;
	int status;
	fftw_complex *local_rhs; 
	fftw_complex *local_work;
	fftwnd_mpi_plan plan_fwd;
	fftwnd_mpi_plan plan_bwd;
	GridType *grid;
	int iComp;
    bool bVel;
    Real gamma;
    BlockInfo * ary;
    
public:
	PeriodicPoisson3D_FFTW_MPI(GridType *_grid, const int _iComp, const bool _bVel, BlockInfo * _ary, const Real _gamma=0): ary(_ary), gamma(_gamma)
	{
		// Initialize some variables
		grid = _grid;
		status = 0;
        
		nx = grid->getBlocksPerDimension(0)*_BLOCKSIZE_;
		ny = grid->getBlocksPerDimension(1)*_BLOCKSIZE_;
		nz = grid->getBlocksPerDimension(2)*_BLOCKSIZE_;
		
		nbpercpux = grid->getResidentBlocksPerDimension(0);
		nbpercpuy = grid->getResidentBlocksPerDimension(1);
		nbpercpuz = grid->getResidentBlocksPerDimension(2);		
		
		int size[dim] = {nx,ny,nz};
        
		// Create plans
		plan_fwd = fftwnd_mpi_create_plan(MPI_COMM_WORLD, dim, size, FFTW_FORWARD, FFTW_ESTIMATE);
		plan_bwd = fftwnd_mpi_create_plan(MPI_COMM_WORLD, dim, size, FFTW_BACKWARD, FFTW_ESTIMATE);
		
		// Get sizes
		fftwnd_mpi_local_sizes(plan_fwd, &local_nx, &local_x_start, &local_ny_after_transpose, &local_y_start_after_transpose, &total_local_size);
        
		//Allocate fftw vectors
		local_rhs  = (fftw_complex*) malloc(sizeof(fftw_complex) * total_local_size);
		local_work = (fftw_complex*) malloc(sizeof(fftw_complex) * total_local_size);
		
		iComp = _iComp;
        bVel = _bVel;
	}
	
	~PeriodicPoisson3D_FFTW_MPI(){}
	
	void Free(void)
	{
		free(local_rhs);
		free(local_work);
		fftwnd_mpi_destroy_plan(plan_fwd);
		fftwnd_mpi_destroy_plan(plan_bwd);
	}
    
    void operator()(blocked_range<int> range) const
    {
        for(int r=range.begin(); r<range.end(); ++r)
        {
            FluidBlock & o = *(FluidBlock *)ary[r].ptrBlock; 
            BlockInfo info = ary[r];
            
            const int idx[dim] = {info.index[0],info.index[1],info.index[2]};
            
            if (bVel)
                switch (status) {
                    case 0:
                    {
						for(int iz=0; iz<BlockType::sizeZ; iz++)
                            for(int iy=0; iy<BlockType::sizeY; iy++)
                                for(int ix=0; ix<BlockType::sizeX; ix++)
                                {
                                    const size_t IDX = (size_t)(idx[2]*BlockType::sizeZ+iz) +
                                    (size_t)(idx[1]*BlockType::sizeY+iy)*BlockType::sizeZ*nbpercpuz +
                                    (size_t)(idx[0]*BlockType::sizeX+ix)*BlockType::sizeY*BlockType::sizeZ*nbpercpuz*nbpercpuy;
                                    
                                    c_re(local_rhs[IDX]) = - ((iComp==0)? o(ix,iy,iz).u: ( (iComp==1)? o(ix,iy,iz).v: o(ix,iy,iz).w))/o(ix,iy,iz).rho;
                                    c_im(local_rhs[IDX]) = 0.0;
                                }
                    }
                        break;
                    case 1:
                    {
                        const double factor = 1.0/((double)ny*(double)nz*(double)nx);
                        for(int iz=0; iz<BlockType::sizeZ; iz++)
                            for(int iy=0; iy<BlockType::sizeY; iy++)
                                for(int ix=0; ix<BlockType::sizeX; ix++)
                                {
                                    const int IDX = (size_t)(idx[2]*BlockType::sizeZ+iz) +
                                    (size_t)(idx[1]*BlockType::sizeY+iy)*BlockType::sizeZ*nbpercpuz +
                                    (size_t)(idx[0]*BlockType::sizeX+ix)*BlockType::sizeY*BlockType::sizeZ*nbpercpuz*nbpercpuy;
                                    
                                    if (iComp==0)
                                    {
                                        assert(!isnan(c_re(local_rhs[IDX])*factor));
                                        o(ix,iy,iz).u = c_re(local_rhs[IDX])*factor;
                                    }
                                    else if (iComp==1)
                                    {
                                        assert(!isnan(c_re(local_rhs[IDX])*factor));
                                        o(ix,iy,iz).v = c_re(local_rhs[IDX])*factor;
                                    }
                                    else if (iComp==2)
                                    {
                                        assert(!isnan(c_re(local_rhs[IDX])*factor));
                                        o(ix,iy,iz).w = c_re(local_rhs[IDX])*factor;
                                    }
                                    else
                                        abort();
                                }
                    }
                        break;
                    default:
                        break;
                }
            else
                switch (status) {
                    case 0:
                    {
                        for(int iz=0; iz<BlockType::sizeZ; iz++)
                            for(int iy=0; iy<BlockType::sizeY; iy++)
                                for(int ix=0; ix<BlockType::sizeX; ix++)
                                {
                                    const size_t IDX = (size_t)(idx[2]*BlockType::sizeZ+iz) +
                                    (size_t)(idx[1]*BlockType::sizeY+iy)*BlockType::sizeZ*nbpercpuz +
                                    (size_t)(idx[0]*BlockType::sizeX+ix)*BlockType::sizeY*BlockType::sizeZ*nbpercpuz*nbpercpuy;
                                    
                                    c_re(local_rhs[IDX]) =  -(gamma-1)/gamma * o.tmp[iz][iy][ix][0];
                                    c_im(local_rhs[IDX]) = 0.0;
                                }
                    }
                        break;
                    case 1:
                    {
                        const double factor = 1.0/((double)ny*(double)nz*(double)nx);
                        for(int iz=0; iz<BlockType::sizeZ; iz++)
                            for(int iy=0; iy<BlockType::sizeY; iy++)
                                for(int ix=0; ix<BlockType::sizeX; ix++)
                                {                                   
                                    const size_t IDX = (size_t)(idx[2]*BlockType::sizeZ+iz) +
                                    (size_t)(idx[1]*BlockType::sizeY+iy)*BlockType::sizeZ*nbpercpuz +
                                    (size_t)(idx[0]*BlockType::sizeX+ix)*BlockType::sizeY*BlockType::sizeZ*nbpercpuz*nbpercpuy;
                                    
                                    // for radius=0.666
                                    // 1.0 -> M=0.5
                                    // .30 -> M=1.1
                                    // .20 -> M=1.5
                                    // .15 -> M=2.0
                                    // .12 -> M=2.8
                                    // .10 -> M=4.3
                                    
                                    // for radius=0.333
                                    // .25 -> M=0.52
                                    // .20 -> M=0.59
                                    // .10 -> M=0.85
                                    // .04 -> M=1.42
                                    // .025 -> M=2
                                    
                                    assert(!isnan(c_re(local_rhs[IDX])*factor));

                                    
                                    const double polya = 0.2;
                                    const double r_pressure_rho = c_re(local_rhs[IDX])*factor+polya;
                                    const double pressure = pow(r_pressure_rho/pow(polya,1./gamma),gamma/(gamma-1.));
                                    const double density = pressure/r_pressure_rho;
                                    
                                    o(ix,iy,iz).rho = density;
                                    o(ix,iy,iz).u *= density;
                                    o(ix,iy,iz).v *= density;
                                    o(ix,iy,iz).w *= density;
                                    
                                    o(ix,iy,iz).energy = pressure/(gamma-1) + 0.5*(o(ix,iy,iz).u*o(ix,iy,iz).u+o(ix,iy,iz).v*o(ix,iy,iz).v+o(ix,iy,iz).w*o(ix,iy,iz).w)/o(ix,iy,iz).rho;
                                }
                    }
                        break;
                    default:
                        break;
                }
        }
    }
    
    void remove_average()
    {       
        double average_re, average_im, average_re_Global, average_im_Global;
        
        for (size_t z=0; z<nz; z++)
            for (size_t y=0; y<ny; y++)
                for (size_t x=0; x<(size_t)local_nx; x++)					
                {
                    const size_t index = (size_t)z + (size_t)(y*nz) + (size_t)(x*nz*ny);
                    average_re+=c_re(local_rhs[index]);
                    average_im+=c_im(local_rhs[index]);
                }
        
        MPI_Allreduce(&average_re, &average_re_Global, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&average_im, &average_im_Global, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        
        const bool isroot = MPI::COMM_WORLD.Get_rank()==0;
        if (isroot) cout << "The average is " << average_re_Global/((double)nx*(double)ny*(double)nz) << endl;
        
        for (size_t z=0; z<nz; z++)
            for (size_t y=0; y<ny; y++)
                for (size_t x=0; x<(size_t)local_nx; x++)					
                {
                    const size_t index = (size_t)z + (size_t)(y*nz) + (size_t)(x*nz*ny);
                    c_re(local_rhs[index]) -= average_re_Global/((double)nx*(double)ny*(double)nz);
                    c_im(local_rhs[index]) -= average_im_Global/((double)nx*(double)ny*(double)nz);
                }
    }
    
    void LaplacianInFourier_2nd()
    {       
        fftwnd_mpi(plan_fwd, 1, local_rhs, local_work, FFTW_NORMAL_ORDER);
        
        const double factor = pow(1.0/(double)max(nx,max(ny,nz)),2.0);
        
        for (size_t z=0; z<nz; z++)
            for (size_t y=0; y<ny; y++)
                for (size_t x=0; x<local_nx; x++)					
                {
                    const size_t index = (size_t)z + (size_t)(y*nz) + (size_t)(x*nz*ny);
                    
                    const double denom = cos(2.0*M_PI*(double)(x+local_x_start)/(double)nx)+cos(2.0*M_PI*(double)y/(double)ny)+cos(2.0*M_PI*(double)z/(double)nz)-3.0;
                    const double inv_denom = (denom==0)? 0.0:1.0/denom;
                    c_re(local_rhs[index]) = 0.5 * c_re(local_rhs[index]) * inv_denom * factor;
                    c_im(local_rhs[index]) = 0.5 * c_im(local_rhs[index]) * inv_denom * factor;
                }
        
        fftwnd_mpi(plan_bwd, 1, local_rhs, local_work, FFTW_NORMAL_ORDER);		
    }
    
    void LaplacianInFourier_4th()
    {        
        fftwnd_mpi(plan_fwd, 1, local_rhs, local_work, FFTW_NORMAL_ORDER);
        
        const bool isroot = MPI::COMM_WORLD.Get_rank()==0;
        if (isroot) cout << "The zeroth element of the FF-transformed array is " << c_re(local_rhs[0]) << endl;
        
        const double factor = pow(1./(double)max(nx,max(ny,nz)),2.);
        
        for (size_t z=0; z<nz; z++)
            for (size_t y=0; y<ny; y++)
                for (size_t x=0; x<(size_t)local_nx; x++)			
                {
                    const size_t index = (size_t)z + (size_t)(y*nz) + (size_t)(x*nz*ny);
                    
                    const double denom = 32.*(cos(2.*M_PI*(double)(x+local_x_start)/(double)nx)+cos(2.*M_PI*(double)y/(double)ny)+cos(2.0*M_PI*(double)z/(double)nz))
                    -2.*(cos(4.*M_PI*(double)(x+local_x_start)/(double)nx)+cos(4.*M_PI*(double)y/(double)ny)+cos(4.*M_PI*(double)z/(double)nz))		
                    -90.;
                    const double inv_denom = (denom==0)? 0.:1./denom;
                    c_re(local_rhs[index]) = 12. * c_re(local_rhs[index]) * inv_denom * factor;
                    c_im(local_rhs[index]) = 12. * c_im(local_rhs[index]) * inv_denom * factor;
                }
        
        fftwnd_mpi(plan_bwd, 1, local_rhs, local_work, FFTW_NORMAL_ORDER);
    }
    
    void Fourier2Space (void)
    {
        assert(status==0);
        status = 1;
    }
    
    void Space2Fourier (void)
    {
        assert(status==1);
        status=0;
    }
};
