/*
 *  SurfaceTension_Test.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 2/20/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <typeinfo>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

#include "TestTypes.h"
#include "Timer.h"

class SurfaceTension_Test
{	
	Real _heaviside(const Real phi, const Real EPSILON)
	{
		const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/EPSILON)));
		const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
		const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
		return (x<0 ? val_xneg : val_xpos);    
	}
	
	inline Real _getGamma(Real gamma1, Real gamma2, Real phi, const Real EPSILON)
	{		
		const Real HS = _heaviside(phi, EPSILON);
		assert(HS>=0);
		assert(HS<=1);
		assert(gamma1>0);
		assert(gamma2>0);
		assert(gamma1*HS + gamma2*(((Real)1)-HS)>=1);
		return (gamma1*HS + gamma2*(((Real)1)-HS));
	}
	
	void _initialize(Block& lab)
	{
		for(int iz = 0; iz<_BLOCKSIZE_; iz++)
			for(int iy = 0; iy<_BLOCKSIZE_; iy++)
				for(int ix = 0; ix<_BLOCKSIZE_; ix++)
				{
					const double val = 1;
					
					lab(ix, iy, iz).dsdt.r = val;
					lab(ix, iy, iz).dsdt.u = val;
					lab(ix, iy, iz).dsdt.v = val;
					lab(ix, iy, iz).dsdt.w = val;
					lab(ix, iy, iz).dsdt.levelset = val;
					lab(ix, iy, iz).dsdt.s = val;
				}		
	}
	
	void _initialize_lab(TestLab_S2& lab, const double h, const double eps, const double radius, const double gamma1, const double gamma2)
	{
		for(int iz = -1; iz<_BLOCKSIZE_+1; iz++)
			for(int iy = -1; iy<_BLOCKSIZE_+1; iy++)
				for(int ix = -1; ix<_BLOCKSIZE_+1; ix++)
				{
					const double x = (ix+0.5)*h;
					const double y = (iy+0.5)*h;
					const double z = (iz+0.5)*h;
					
					const double r = sqrt(pow((x-0.5)/0.2,2)+pow((y-0.5)/0.12,2)+pow((z-0.5)/0.15,2));
					const double bubble = _heaviside(r-radius, eps);
					
					lab(ix, iy, iz).s.r = 10*bubble+0.01*(1.-bubble);
					lab(ix, iy, iz).s.u = x*x;
					lab(ix, iy, iz).s.v = (1-y)*(1-y);
					lab(ix, iy, iz).s.w = z*z;
					lab(ix, iy, iz).s.levelset = 1./(_getGamma(gamma1, gamma2, r-radius, eps)-1);
					lab(ix, iy, iz).s.s = 10*lab(ix, iy, iz).s.levelset;
				}
	}
	
	void _gold(TestLab_S2& lab, Block& block, const Real a, const Real _dtinvh, const Real G1, const Real G2, const Real _h, const Real _sigma=1);
	
	void _compare(Block& _a, Block& _b, double accuracy, string kernelname);
	
	template<typename FS> double _benchmark(FS fs, const int NBLOCKS, const int NTIMES)
	{
		TestLab_S2 * lab = new TestLab_S2[NBLOCKS];
		Block * block = new Block[NBLOCKS];
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize_lab(lab[i], 1./_BLOCKSIZE_, 1, 0.1, 1/fs.G1+1, 1/fs.G2+1);
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(block[i]);
		
		double tCOMPUTE = 0;
		
		//measure performance
		{
			Timer timer;
						
			//run FS
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+2;
			const int slicesrcs =  (_BLOCKSIZE_+2)*(_BLOCKSIZE_+2);
			
			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;
			
			const Real * const srcfirst = &(lab[0])(-1,-1, -1).s.r;
			Real * const dstfirst = &(block[0])(0,0,0).dsdt.r;
			
			if (NBLOCKS == 1)
				fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
			
			timer.start();
			for(int i=0; i<NTIMES; i++)
			{
				const Real * const srcfirst = &(lab[i%NBLOCKS])(-1,-1, -1).s.r;
				Real * const dstfirst = &(block[i%NBLOCKS])(0,0,0).dsdt.r;
				
				fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
			}
			tCOMPUTE = timer.stop();
		}
		
		delete [] lab;
		delete [] block;
		
		return tCOMPUTE;
	}
	
	double _benchmarkGold(const int NBLOCKS, const int NTIMES, const double a, const double dtinvh, const double sigma, const double h, const double radius, const double gamma1, const double gamma2)
	{
		TestLab_S2 * lab = new TestLab_S2[NBLOCKS];
		Block * blockgold = new Block[NBLOCKS];
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize_lab(lab[i], 1./_BLOCKSIZE_, 1., radius, gamma1, gamma2);
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(blockgold[i]);
		
		double tGOLD = 0;
		
		//measure performance
		{
			Timer timer;
			
			//run gold
			if (NBLOCKS == 1)
				_gold(lab[0], blockgold[0], a, dtinvh, sigma, h, gamma1, gamma2);
			
			timer.start();
			for(int i=0; i<NTIMES; i++)
				_gold(lab[i%NBLOCKS], blockgold[i%NBLOCKS], a, dtinvh, sigma, h, gamma1, gamma2);
			tGOLD = timer.stop();
		}
		
		delete [] lab;
		delete [] blockgold;
		
		return tGOLD;
	}
	
public:
	
	template<typename FS> void accuracy(FS& fs, double accuracy=1e-4, bool=false)
	{
		TestLab_S2 * lab = new TestLab_S2;
		_initialize_lab(*lab, 1./_BLOCKSIZE_, 1., 0.1, 1/fs.G1+1, 1/fs.G2+1);
				
		Block * blockgold = new Block;
		_initialize(*blockgold);
		
		_gold(*lab, *blockgold, fs.a, fs.dtinvh, fs.G1, fs.G2, fs.h, fs.sigma);//34.164);//
		
		Block * block = new Block;
		_initialize(*block);

		{
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+2;
			const int slicesrcs =  (_BLOCKSIZE_+2)*(_BLOCKSIZE_+2);
			
			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;
			
			const Real * const srcfirst = &(*lab)(-1,-1, -1).s.r;
			Real * const dstfirst = &(*block)(0,0,0).dsdt.r;
			
			fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs,
					   dstfirst, dstfloats, rowdsts, slicedsts);
		}
		
		delete lab;
		
		_compare(*blockgold, *block, accuracy, typeid(fs).name());
		
		delete blockgold;
		delete block;		
	}
	
	template<typename FS> void performance(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{
		const Real gamma1 = 1/fs.G1 + 1;
		const Real gamma2 = 1/fs.G2 + 1;
		const Real dtinvh = fs.dtinvh;
		const Real a = fs.a;
		const Real sigma = fs.sigma;
		
		float tGOLD = 0, tCOMPUTE = 0;
		int COUNT = 0;
		
#pragma omp parallel 
		{
#pragma omp critical
			{
				COUNT++;
			}
		}

		//measure performance
		{
			Timer timer;
			
			//run gold
#pragma omp parallel
			{
				const double t =  _benchmarkGold(NBLOCKS, NTIMES, a, dtinvh, sigma, 1./_BLOCKSIZE_, 0.1, gamma1, gamma2);
				
#pragma omp critical
				{
					tGOLD += t;
				}
			}
			
			tGOLD /= COUNT;
			
			//run FS
#pragma omp parallel
			{
				const double t =  _benchmark(fs, NBLOCKS, NTIMES);
				
#pragma omp critical
				{
					tCOMPUTE += t;
				}
			}
			
			tCOMPUTE /= COUNT;
		}
		
		string implname = typeid(fs).name();

		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE, false);
		printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
	}
	
	template<typename FS> void profile(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{	
		const double t =  _benchmark(fs, NBLOCKS, NTIMES);
		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, t, false);
	}
};
