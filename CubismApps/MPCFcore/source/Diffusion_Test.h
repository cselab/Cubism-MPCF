//
//  Diffusion_Test.h
//  MPCFcore
//
//  Created by Babak Hejazialhosseini on 9/16/11.
//  Copyright 2011 ETH Zurich. All rights reserved.
//

//This code is NOT optimized and is a second-order accurate
//shear-stress tensor naive implementation.

#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <string>

#include <BlockInfo.h>
#include "TestTypes.h"
#include "Timer.h"

class Diffusion_Test
{
	void _initialize(TestLab_S2& lab, const double h)
	{
		/* Alternatively, we load the test initial condition from a file:
		 
		 FILE * f = fopen("test0", "r");
		 
		 int c = 0;
		 while (!feof(f)) {
		 float r,u,v,w,s,l;
		 fscanf(f, "%e %e %e %e %e %e\n", &r, &u, &v, &w, &s, &l);
		 
		 const int ix = (c % 18) -1;
		 const int iy = (c/18 % 18) - 1;
		 const int iz = (c/18/18) -1;
		 
		 lab(ix, iy, iz).s.r = r;
		 lab(ix, iy, iz).s.u = u;
		 lab(ix, iy, iz).s.v = v;
		 lab(ix, iy, iz).s.w = w;
		 lab(ix, iy, iz).s.levelset = s;
		 lab(ix, iy, iz).s.s = l;
		 
		 c++;
		 }
		 
		 fclose(f);
		 */
		
		for(int iz = -1; iz<_BLOCKSIZE_+1; iz++)
			for(int iy = -1; iy<_BLOCKSIZE_+1; iy++)
				for(int ix = -1; ix<_BLOCKSIZE_+1; ix++)
				{
					const double x = (ix+0.5)*h;
					const double y = (iy+0.5)*h;
					const double z = (iz+0.5)*h;
					
					const double r = sqrt(pow((x-0.5)/0.2,2)+pow((y-0.5)/0.12,2)+pow((z-0.5)/0.15,2));
					
					lab(ix, iy, iz).s.r = 1 + (x/16 +y/8 + z/4)*0.05;
					lab(ix, iy, iz).s.u = x*x-z;//x*x-z;
					lab(ix, iy, iz).s.v = (1-y)*(1-y)*z;
					lab(ix, iy, iz).s.w = z*z/2;
					lab(ix, iy, iz).s.levelset = 213;
					lab(ix, iy, iz).s.s = x-0.4;
				}
	}
	
	void _initialize(Block& block)
	{
		for(int iz = 0; iz<_BLOCKSIZE_; iz++)
			for(int iy = 0; iy<_BLOCKSIZE_; iy++)
				for(int ix = 0; ix<_BLOCKSIZE_; ix++)
				{
					const double val = 1;
					
					block(ix, iy, iz).dsdt.r = val;
					block(ix, iy, iz).dsdt.u = val;
					block(ix, iy, iz).dsdt.v = val;
					block(ix, iy, iz).dsdt.w = val;
					block(ix, iy, iz).dsdt.levelset = val;
					block(ix, iy, iz).dsdt.s = val;
				}
	}
	
	double _benchmarkGold(const int NBLOCKS, const int NTIMES, const Real nu1, const Real nu2, const Real G1, const Real G2, const Real a, const Real dtinvh, const Real h)
	{
		TestLab_S2 * lab = new TestLab_S2[NBLOCKS];
		Block * blockgold = new Block[NBLOCKS];
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(lab[i], h);
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(blockgold[i]);
		
		double tGOLD = 0;
		
		//measure performance
		{
			Timer timer;
			
			//run gold
			if (NBLOCKS == 1)
				_gold(lab[0], blockgold[0], nu1, nu2, G1, G2, a, dtinvh, h);
			
			timer.start();
			for(int i=0; i<NTIMES; i++)
				_gold(lab[i%NBLOCKS], blockgold[i%NBLOCKS], nu1, nu2, G1, G2, a, dtinvh, h);
			tGOLD = timer.stop();
		}
		
		delete [] lab;
		delete [] blockgold;
		
		return tGOLD;
	}
	
	void _gold(TestLab_S2& lab, Block& block, const Real _nu1, const Real _nu2, const Real _g1, const Real _g2, const Real a, const Real _dtinvh, const Real _h);
	
	void _compare(Block& _a, Block& _b, double accuracy, std::string kernelname);
	
	template<typename FS> double _benchmark(FS fs, const int NBLOCKS, const int NTIMES)
	{
		TestLab_S2 * lab = new TestLab_S2[NBLOCKS];
		Block * block = new Block[NBLOCKS];
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(lab[i], fs.h);
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(block[i]);
		
		double tCOMPUTE = 0;
		
		//measure performance
		{
			Timer timer;
			
			//timer.start();
			
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
	
public:
	
	template<typename FS> void accuracy(FS& fs, double accuracy=1e-4, bool=false)
	{
		TestLab_S2 * lab = new TestLab_S2;
		const double EPSILON = 1e-6;
		_initialize(*lab, fs.h);
		
		Block * blockgold = new Block;
		_initialize(*blockgold);
		
		_gold(*lab, *blockgold, fs.nu1, fs.nu2, fs.G1, fs.G2, fs.a, fs.dtinvh, fs.h);
		
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
		const Real nu1 = fs.nu1;
		const Real nu2 = fs.nu2;
		const Real G1 = fs.G1;
		const Real G2 = fs.G2;
		const Real dtinvh = fs.dtinvh;
		const Real a = fs.a;
		const Real sigma = fs.sigma;
		const Real h = fs.h;
		
		float tGOLD = 0, tCOMPUTE = 0;
		int COUNT = 0;
		
#pragma omp parallel 
		{
#pragma omp critical
			{
				COUNT++;
			}
		}
		printf("NTHREADS is %d\n", COUNT);
		
		
		//measure performance
		{			
			//run gold
#pragma omp parallel
			{
				const double t =  _benchmarkGold(NBLOCKS, NTIMES, nu1, nu2, G1, G2, a, dtinvh, h);
				
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
		
		//report performance
		std::string implname = typeid(fs).name();
		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE, false);
		printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
	}
	
	template<typename FS> void profile(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{		
		int COUNT = 0;
		
#pragma omp parallel 
		{
#pragma omp critical
			{
				COUNT++;
			}
		}
		
		printf("NTHREADS is %d\n", COUNT);
		double tCOMPUTE = 0;
		
#pragma omp parallel
		{
			const double t =  _benchmark(fs, NBLOCKS, NTIMES);
			
#pragma omp critical
			{
				tCOMPUTE += t;
			}
		}
		
		tCOMPUTE /= COUNT;
		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE, false);
	}
};
