/*
 *  Convection_Test.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include <Timer.h>

#include "common.h"
#include "TestTypes.h"

using namespace std;

class Test_Convection
{
	void _toPrimitive(GP& p, Real out[7]) 
	{
		out[0] = p.s.r;
		out[1] = p.s.u/out[0];
		out[2] = p.s.v/out[0];
		out[3] = p.s.w/out[0];
        out[4] = (p.s.s - (p.s.u*p.s.u+p.s.v*p.s.v+p.s.w*p.s.w)*(((Real)0.5)/p.s.r) - p.s.P)/p.s.G;
		out[5] = p.s.G;
        out[6] = p.s.P;

		assert(p.s.r>0);
		assert(!isnan(p.s.r));
		assert(!isnan(p.s.u));
		assert(!isnan(p.s.v));
		assert(!isnan(p.s.w));
		assert(!isnan(p.s.s));
		assert(!isnan(p.s.G));
        assert(!isnan(p.s.P));

  		for(int i=0; i<7; i++)
			assert(!isnan(out[i]));
	}
	
	virtual void _initialize(TestLab& lab, Block& block);
	void _gold(TestLab& lab, Block& block, const Real h);
	void _print(Block& block);

	Real dtinvh;
	
public:
	
	template<typename FS> void accuracy(FS& fs, double accuracy=1e-4, bool bAwk=false)
	{
		dtinvh = fs.dtinvh;
		
		TestLab * lab = new TestLab;
		Block * blockgold = new Block;
		Block * block = new Block;
		
		_initialize(*lab, *blockgold);
		_initialize(*lab, *block);
		_gold(*lab, *blockgold, 1);
		
		//run FS
		{
			const Real * const srcfirst = &(*lab)(-3,-3, -3).s.r;
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+6;
			const int slicesrcs =  (_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);
			
			Real * const dstfirst = &(*block)(0,0,0).dsdt.r;
			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;
			
			fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
		}
		
		printAccuracyTitle();
		string kernelname = typeid(fs).name();
		
		block->compare(*blockgold, accuracy, kernelname);
		printEndLine();
		
		
		delete block;
		delete blockgold;
		delete lab;
	}
	
protected:
	template<typename FS> double _benchmark(FS fs, const int NBLOCKS, const int NTIMES)
	{
		TestLab * lab = new TestLab[NBLOCKS];
		Block * block = new Block[NBLOCKS];
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(lab[i], block[i]);
		
		double tCOMPUTE = 0;
		
		//measure performance
		{
			Timer timer;
			
			timer.start();
			
			//run FS
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+6;
			const int slicesrcs =  (_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);
			
			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;
			
			const Real * const srcfirst = &(lab[0])(-3,-3, -3).s.r;
			Real * const dstfirst = &(block[0])(0,0,0).dsdt.r;
			
			if (NBLOCKS == 1)
				fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
			
			timer.start();
			for(int i=0; i<NTIMES; i++)
				fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
			tCOMPUTE = timer.stop();
		}
		
		delete [] lab;
		delete [] block;
		
		return tCOMPUTE;
	}
	
public:
	
	template<typename FS> void performance(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{
		dtinvh = fs.dtinvh;
		
		int COUNT = 0;
		
#pragma omp parallel 
		{
#pragma omp critical
			{
				COUNT++;
			}
		}
		
		double tGOLD = 0, tCOMPUTE = 0;
		
#pragma omp parallel 
		{
			TestLab * lab = new TestLab[NBLOCKS];
			Block * blockgold = new Block[NBLOCKS];
			
			for(int i=0; i<NBLOCKS; i++)
				_initialize(lab[i], blockgold[i]);
			
			Timer timer;
			
			//run gold
			if (NBLOCKS == 1)
				_gold(lab[0], blockgold[0], 1);
			
#pragma omp barrier			
			
			timer.start();
			for(int i=0; i<NTIMES; i++)
				_gold(lab[i%NBLOCKS], blockgold[i%NBLOCKS], 1);
			const double tgold = timer.stop();
			
#pragma omp barrier
			
			const double t = _benchmark(fs, NBLOCKS, NTIMES);
			
#pragma omp barrier
			
#pragma omp critical			
			{
				tGOLD += tgold;
				tCOMPUTE += t;
			}
		}
		
		tCOMPUTE /= COUNT;
		tGOLD /= COUNT;
		
		
		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE);
		
		printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
	}
	
	template<typename FS> void profile(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{
		dtinvh = fs.dtinvh;
		
		int COUNT = 0;
		
#pragma omp parallel 
		{
#pragma omp critical
			{
				COUNT++;
			}
		}
		
		double tCOMPUTE = 0;
		
#pragma omp parallel 
		{
			TestLab * lab = new TestLab[NBLOCKS];
			Block * block = new Block[NBLOCKS];
			
			for(int i=0; i<NBLOCKS; i++)
				_initialize(lab[i], block[i]);
			
			//float tCOMPUTE = 0;
			
			const double t = _benchmark(fs, NBLOCKS, NTIMES);
			
#pragma omp critical			
			{
				tCOMPUTE += t;
			}
			
			delete block;
			delete lab;	
		}
		
		tCOMPUTE /= COUNT;
		
		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE);
	}
};
