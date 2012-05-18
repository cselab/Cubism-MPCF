/*
 *  FlowStep_Test.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include <Timer.h>
#include "TestTypes.h"

#include "Convection_CPP.h"
//#ifdef _SSE_
//#include "FlowStep_SSE_diego.h"
//#endif
//#if _ALIGNBYTES_ % 32 == 0
//#include "FlowStep_AVX_diego.h"
//#endif

using namespace std;

class FlowStep_Test
{
protected:
	Real gamma1, gamma2, smoothlength, dtinvh;
	
	static Real _heavyside(Real phi, Real h)
	{
		const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/h)));
		const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
		const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
		return x<0 ? val_xneg : val_xpos;
	}
	
	inline Real _getGamma(Real phi)
	{		
		const Real HS = _heavyside(phi, smoothlength);
		assert(HS>=0);
		assert(HS<=1);
		assert(gamma1>0);
		assert(gamma2>0);
		assert(gamma1*HS + gamma2*(((Real)1)-HS)>=1);
		return (gamma1*HS + gamma2*(((Real)1)-HS));
	}
	
	void _toPrimitive(GP& p, Real out[6]) 
	{
		out[0] = p.s.r;
		out[1] = p.s.u/out[0];
		out[2] = p.s.v/out[0];
		out[3] = p.s.w/out[0];
		
		out[4] = (p.s.s - (p.s.u*p.s.u+p.s.v*p.s.v+p.s.w*p.s.w)*(((Real)0.5)/p.s.r))*(_getGamma(p.s.levelset)-(Real)1);
		out[5] = p.s.levelset;
		
		assert((p.s.s > (pow(p.s.u,2)+pow(p.s.v,2)+pow(p.s.w,2))*(0.5/p.s.r)));
		assert((_getGamma(p.s.levelset)-1.) > 0);
		assert(p.s.r>0);
		assert(!isnan(p.s.r));
		assert(!isnan(p.s.u));
		assert(!isnan(p.s.v));
		assert(!isnan(p.s.w));
		assert(!isnan(p.s.s));
		assert(!isnan(p.s.levelset));
		
		for(int i=0; i<6; i++)
			assert(!isnan(out[i]));
		
		assert(out[4]>0);
	}
	
	virtual void _initialize(TestLab& lab, Block& block);
	void _gold(TestLab& lab, Block& block, const Real h);
	void _print(Block& block);
	void _compare(Block& a, Block& b, double accuracy=1e-4, bool bAwk=false, string kernelname="");
	
public:
	
	template<typename FS> void accuracy(FS& fs, double accuracy=1e-4, bool bAwk=false)
	{
		gamma1 = fs.gamma1;
		gamma2 = fs.gamma2;
		smoothlength = fs.smoothlength;
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
		/*if (bAwk)
		{
			if (typeid(fs) == typeid(FlowStep_CPP))
				kernelname = "FlowStep_CPP";
#ifdef _SSE_
			else if (typeid(fs) == typeid(FlowStep_SSE_diego))
				kernelname = "FlowStep_SSE_diego";
#endif
#if _ALIGNBYTES_ % 32 == 0 && defined(_AVX_)
			else if (typeid(fs) == typeid(FlowStep_AVX_diego))
				kernelname = "FlowStep_AVX_diego";
#endif
			else
			{
				cout << "Kernel name not inserted for awk\n";
				abort();
			}
		}*/
		
		_compare(*block, *blockgold, accuracy, bAwk, kernelname);
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
		gamma1 = fs.gamma1;
		gamma2 = fs.gamma2;
		smoothlength = fs.smoothlength;
		dtinvh = fs.dtinvh;
		
		TestLab * lab = new TestLab[NBLOCKS];
		Block * blockgold = new Block[NBLOCKS];
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(lab[i], blockgold[i]);
		
		float tGOLD = 0, tCOMPUTE = 0;
		int COUNT = 0;
		
#pragma omp parallel 
		{
#pragma omp critical
			{
				COUNT++;
			}
		}
//		printf("NTHREADS is %d\n", COUNT);

		
		//measure performance
		{
			Timer timer;
			
			//run gold
			if (NBLOCKS == 1)
				_gold(lab[0], blockgold[0], 1);
			
			timer.start();
			for(int i=0; i<NTIMES; i++)
				_gold(lab[i%NBLOCKS], blockgold[i%NBLOCKS], 1);
			tGOLD = timer.stop();
			
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
/*
		{
			if (typeid(fs) == typeid(FlowStep_CPP))
				implname = "FlowStep_CPP";
#ifdef _SSE_
			else if (typeid(fs) == typeid(FlowStep_SSE_diego))
				implname = "FlowStep_SSE_diego";
#endif
#if _ALIGNBYTES_ % 32 == 0 && defined(_AVX_)
			else if (typeid(fs) == typeid(FlowStep_AVX_diego))
				implname = "FlowStep_AVX_diego";
#endif
			else
			{
				cout << "Kernel name not inserted for awk\n";
				abort();
			}
		}*/
		
		fs.printflops(PEAKPERF, PEAKBAND, 1, NTIMES, 1, tCOMPUTE);
		//cout << setprecision(4) << "\tMEMORY FOOTPRINT: "<< (lab->kB() +  Block::kB() + fs.kB(false))/1024<< " MB" << endl;
		
		//how to solve this??????
		//printf("\tGOLD IS %.2f GFLOP/s\n", FlowStep_CPP::getGFLOP()/tGOLD);	
		printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
		
		delete [] blockgold;
		delete [] lab;	
	}
	
	template<typename FS> void profile(FS& fs, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{
		gamma1 = fs.gamma1();
		gamma2 = fs.gamma2();
		smoothlength = fs.smoothlength();
		dtinvh = fs.dtinvh();
		
		TestLab * lab = new TestLab[NBLOCKS];
		Block * block = new Block[NBLOCKS];
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(lab[i], block[i]);
		
		float tCOMPUTE = 0;
		
		//measure performance
		{
			Timer timer;
			
			//run gold
			timer.start();
			
			//run FS
			const int srcfloats  = sizeof(GP)/sizeof(Real);
			const int rowsrcs = _BLOCKSIZE_+6;
			const int slicesrcs =  (_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);
			
			const int dstfloats =  sizeof(GP)/sizeof(Real);
			const int rowdsts =  _BLOCKSIZE_;
			const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;
			
			timer.start();
			
			for(int i=0; i<NTIMES; i++)
			{
				const Real * const srcfirst = &(lab[i%NBLOCKS])(-3,-3, -3).s.r;
				Real * const dstfirst = &(block[i%NBLOCKS])(0,0,0).dsdt.r;
				
				fs.compute(srcfirst, srcfloats, rowsrcs, slicesrcs, dstfirst, dstfloats, rowdsts, slicedsts);
			}
			
			tCOMPUTE = timer.stop();
		}
		
		//FLOP estimation
		const double GFLOPCONVERT = 31*(_BLOCKSIZE_+6)*(_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);
		const double GFLOPFLUX = 911*_BLOCKSIZE_*(_BLOCKSIZE_+1);
		const double GFLOPRHS = (12+12+6)*_BLOCKSIZE_*_BLOCKSIZE_;
		const double GFLOPCOPYBACK = 9*_BLOCKSIZE_*_BLOCKSIZE_;
		const double GFLOP = NTIMES*1.*(GFLOPCONVERT + GFLOPFLUX + _BLOCKSIZE_*(GFLOPFLUX*3 + GFLOPRHS + GFLOPCOPYBACK))/1e9;
		
		//FLOP/s estimation
		const double AIConvert   = .5;
		const double AIWENO      = 2.5;
		const double AIExtraTerm = .125;
		const double AICharVel   = 1.3;
		const double AIHLLERho   = .571;
		const double AIHLLEVel   = .25;
		const double AIHLLEPVel  = .584;
		const double AIHLLEE     = 1.283;
		const double AIPRHS      = .125;
		const double AICopyBack  = .161;
		const double EPERFCONVERT  = min(AIConvert*  PEAKBAND, PEAKPERF);
		const double EPERWENO      = min(AIWENO*     PEAKBAND, PEAKPERF);
		const double EPERXTRATERM  = min(AIExtraTerm*PEAKBAND, PEAKPERF);
		const double EPERCHARVEL   = min(AICharVel*  PEAKBAND, PEAKPERF);
		const double EPERHLLERHO   = min(AIHLLERho*  PEAKBAND, PEAKPERF);
		const double EPERHLLEVEL   = min(AIHLLEVel*  PEAKBAND, PEAKPERF);
		const double EPERHLLEPVEL  = min(AIHLLEPVel* PEAKBAND, PEAKPERF);
		const double EPERHLLEE     = min(AIHLLEE*    PEAKBAND, PEAKPERF);
		const double EPERFPRHS     = min(AIPRHS*     PEAKBAND, PEAKPERF); 
		const double EPERFCOPYBACK = min(AICopyBack* PEAKBAND, PEAKPERF);
		
		//execution time estimation
		const double TCONVERT= GFLOPCONVERT/EPERFCONVERT;
		const double TFLUX = _BLOCKSIZE_*(_BLOCKSIZE_+1)*(60*12./EPERWENO + 8./EPERXTRATERM + 52./EPERCHARVEL + 16*2./EPERHLLERHO + 21./EPERHLLEPVEL + 2*14./EPERHLLEVEL + 77./EPERHLLEE);//GFLOPFLUX/EPERFFLUX;
		const double TRHS = GFLOPRHS/EPERFPRHS;
		const double TCOPYBACK = GFLOPCOPYBACK/EPERFCOPYBACK;
		const double TEXPECTED = NTIMES*1.*(TCONVERT + TFLUX + _BLOCKSIZE_*(TFLUX*3 + TRHS + TCOPYBACK));
		
		printPerformanceTitle();
		cout << setprecision(4) << "\tMEMORY FOOTPRINT: "<< (lab->kB() +  block->kB() + fs.kB(false))/1024<< " MB" << endl;
		printf("\tASSUMING PP: %.2f GFLOP/s  PB: %.2f GB/s\n", PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tRIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tTHIS ONE IS %.2f GFLOP/s\n", GFLOP/tCOMPUTE);	
		printf("\tTIME PER BLOCK: %.2f ms (expected %.2f ms)\n",  1e3*tCOMPUTE/NTIMES, 1e3*TEXPECTED/NTIMES);
		printf("\tEFFICIENCY: %.2f%%, HW-UTILIZATION: %.2f%%\n", 100.*TEXPECTED/tCOMPUTE, 100*(GFLOP/tCOMPUTE*1e9)/PEAKPERF);
		if (bAwk)
		{
			cout << "no awk output there!\n";
			abort();
		}
		printEndLine();
		
		delete block;
		delete lab;	
	}
};

class FlowStep_Test_Babak : public FlowStep_Test
{
	
	static Real _heavyside(Real phi, Real h)
	{
		const Real epsilon = h;
		const Real relDist = phi / epsilon;
		if(relDist>1.0) return 0.0;
		if(relDist<-1.0) return 1.0;
		if(relDist<0.0) 
		{
			const Real mapx = 1.0+relDist;
			return -mapx*mapx*0.5+1.0;
		} 
		else
		{
			const Real mapx = relDist-1.0;
			return 0.5*mapx*mapx;
		}
	}
	
	inline Real _getGamma(Real phi)
	{		
		const Real HS = _heavyside(phi, smoothlength);
		assert(HS>=0);
		assert(HS<=1);
		assert(gamma1>0);
		assert(gamma2>0);
		assert(gamma1*HS + gamma2*(1.0-HS)>1);
		return (gamma1*HS + gamma2*(1.0-HS));
	}
	
	void _initialize(TestLab& lab, Block& block);
};
