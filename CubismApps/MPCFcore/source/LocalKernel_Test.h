/*
 *  MaxSpeedOfSound_Test.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 6/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#include "MaxSpeedOfSound.h"
#include "Update.h"
#include "output.h"

class LocalKernel_Test
{
protected:
	Real gamma1, gamma2, smoothlength;
	
	virtual void _initialize(Block& block)
	{
		srand48(61651);
		
		for(int iz = 0; iz<_BLOCKSIZE_; iz++)
			for(int iy = 0; iy<_BLOCKSIZE_; iy++)
				for(int ix = 0; ix<_BLOCKSIZE_; ix++)
				{
					//memset(&block(ix, iy, iz), 0, sizeof(GP));
					
					block(ix, iy, iz).s.r = drand48()+iz;
					block(ix, iy, iz).s.u = drand48()+ix;
					block(ix, iy, iz).s.v = drand48()+iy;
					block(ix, iy, iz).s.w = drand48()+(ix*iy/(double) _BLOCKSIZE_);
					block(ix, iy, iz).s.s = drand48()+1+(iy*iz)/(double) _BLOCKSIZE_;
					block(ix, iy, iz).s.levelset = -1 +  (iz+ix+iy)/(double) _BLOCKSIZE_;
					
					block(ix, iy, iz).dsdt.r = drand48()+iz;
					block(ix, iy, iz).dsdt.u = drand48()+ix;
					block(ix, iy, iz).dsdt.v = drand48()+iy;
					block(ix, iy, iz).dsdt.w = drand48()+(ix*iy/(double) _BLOCKSIZE_);
					block(ix, iy, iz).dsdt.s = drand48()+1+(iy*iz)/(double) _BLOCKSIZE_;
					block(ix, iy, iz).dsdt.levelset = -1 +  (iz+ix+iy)/(double) _BLOCKSIZE_;
				}
	}
	
	virtual void _initialize_from_file(Block& block)
	{
		FILE * f = fopen("test_rank0", "r");
		
		int c = 0;
		while (!feof(f)) {
			float r,u,v,w,s,l;
			const int values = fscanf(f, "%e %e %e %e %e %e\n", &r, &u, &v, &w, &s, &l);
			assert(values==6);

			const int ix = (c % 16);
			const int iy = (c/16 % 16) ;
			const int iz = (c/16/16);
			
			block(ix, iy, iz).s.r = r;
			block(ix, iy, iz).s.u = u;
			block(ix, iy, iz).s.v = v;
			block(ix, iy, iz).s.w = w;
			block(ix, iy, iz).s.levelset = s;
			block(ix, iy, iz).s.s = l;
			
			c++;
		}
		
		fclose(f);
	}
	
	
	void _compare(Block& _a, Block& _b, double accuracy=1e-4, bool bAwk=false, string kernelname="")
	{
		double maxe[6]={0,0,0,0,0,0};
		double sume[6]={0,0,0,0,0,0};
		
		for(int iz = 0; iz<_BLOCKSIZE_; iz++)
			for(int iy = 0; iy<_BLOCKSIZE_; iy++)
				for(int ix = 0; ix<_BLOCKSIZE_; ix++)
				{
					StateVector a = _a(ix, iy, iz).dsdt;
					StateVector b = _b(ix, iy, iz).dsdt;
					const double s[6]  = {
						b.r ,
						b.u ,
						b.v ,
						b.w ,
						b.s ,
						b.levelset
					};
					
					const double e[6]  = {
						b.r - a.r,
						b.u - a.u,
						b.v - a.v,
						b.w - a.w,
						b.s - a.s,
						b.levelset - a.levelset
					};
					
					for(int i=0; i<6; i++)
						if (fabs(e[i])/fabs(s[i])>accuracy && fabs(e[i])>accuracy) printf("significant error at %d %d %d -> e=%e (rel is %e)\n", ix, iy, iz, e[i], e[i]/s[i]);
					
					for(int i=0; i<6; i++)
						maxe[i] = max(fabs(e[i]), maxe[i]);
					
					for(int i=0; i<6; i++)
						sume[i] += fabs(e[i]);
				}
		
		printf("\tLinf discrepancy:\t");
		for(int i=0; i<6; i++)
			printf("%.2e ", maxe[i]);
		
		cout << endl;
		printf("\tL1 (dh=1):       \t");
		for(int i=0; i<6; i++)
			printf("%.2e ", sume[i]);
		cout<<endl;
		
		if (bAwk) awkAcc(kernelname,
						 maxe[0], maxe[1], maxe[2], maxe[3], maxe[4], maxe[5],
						 sume[0], sume[1], sume[2], sume[3], sume[4], sume[5]);
	}
	
public:
	
	template<typename TSOS>
	void accuracy(TSOS& kernel, MaxSpeedOfSound_CPP& refkernel, double accuracy=1e-4, bool bAwk=false)
	{
		Block * blockgold = new Block;
		Block * block = new Block;
		
		_initialize(*blockgold);
		_initialize(*block);
		
		const Real v1 = refkernel.compute(&(*blockgold)(0,0,0).s.r, sizeof(GP)/sizeof(Real));
		
		//run Kernel
		const Real v2 = kernel.compute(&(*blockgold)(0,0,0).s.r, sizeof(GP)/sizeof(Real));
					
		printAccuracyTitle();
		printf("ERROR: %e (relative: %e)\n", v1-v2, (v1-v2)/max(fabs(v1), max(fabs(v2), (Real)accuracy)));
		printEndLine();
		
		delete block;
		delete blockgold;
	}
	
	template<typename TUPDATE>
	void accuracy(TUPDATE& kernel, Update_CPP& refkernel, double accuracy=1e-4, bool bAwk=false)
	{
		Block * blockgold = new Block;
		Block * block = new Block;
		
		_initialize(*blockgold);
		_initialize(*block);
		
		refkernel.compute(&(*blockgold)(0,0,0).dsdt.r, &(*blockgold)(0,0,0).s.r, sizeof(GP)/sizeof(Real));
		
		//run Kernel
		kernel.compute(&(*block)(0,0,0).dsdt.r, &(*block)(0,0,0).s.r, sizeof(GP)/sizeof(Real));
		
		printAccuracyTitle();
		_compare(*block, *blockgold, accuracy, bAwk, "Update_SSE");
		printEndLine();
		
		delete block;
		delete blockgold;
	}
	
	template<typename TKernel> double _benchmarkSOS(TKernel kernel, const int NBLOCKS, const int NTIMES)
	{
		Block * block = new Block[NBLOCKS];
	
		for(int i=0; i<NBLOCKS; i++)
			_initialize(block[i]);
		
		Timer timer;
		timer.start();
		
		if (NBLOCKS == 1)
			kernel.compute(&block[0](0,0,0).s.r, sizeof(GP)/sizeof(Real));

		for(int i=0; i<NTIMES; i++)
			kernel.compute(&block[i%NBLOCKS](0,0,0).s.r, sizeof(GP)/sizeof(Real));
		
		const double tCOMPUTE = timer.stop();
			
		delete [] block;
		
		return tCOMPUTE;
	}
	
	template<typename TKernel> double _benchmarkUP(TKernel kernel, const int NBLOCKS, const int NTIMES)
	{
		Block * block = new Block[NBLOCKS];
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(block[i]);
		
		Timer timer;
		timer.start();
		
		if (NBLOCKS == 1)
			kernel.compute(&block[0](0,0,0).dsdt.r, &block[0](0,0,0).s.r, sizeof(GP)/sizeof(Real));		
		
		for(int i=0; i<NTIMES; i++)
			kernel.compute(&block[i%NBLOCKS](0,0,0).dsdt.r, &block[i%NBLOCKS](0,0,0).s.r, sizeof(GP)/sizeof(Real));		
		
		const double tCOMPUTE = timer.stop();
		
		delete [] block;
		
		return tCOMPUTE;
	}
	
	
	template<typename TSOS>
	void performance(TSOS& kernel, MaxSpeedOfSound_CPP& refkernel, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{
		Block * blockgold = new Block[NBLOCKS];
		Block * block = new Block[NBLOCKS];
		
		const double OI = 120./24.;
		const double GFLOP = _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*120.*NTIMES/1.e9;
		const double EPERF = min(OI*PEAKBAND,PEAKPERF);
		const double TEXPECTED = 1.e9*GFLOP/EPERF;
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(blockgold[i]);
	
		float tGOLD = 0, tCOMPUTE = 0;
		int COUNT = 0;
		
#pragma omp parallel 
		{
#pragma omp critical
			{
				COUNT++;
			}
		}
		//printf("NTHREADS is %d\n", COUNT);
		
		//measure performance
		{
			Timer timer;
			
			if (NBLOCKS == 1)
				refkernel.compute(&blockgold[0](0,0,0).s.r, sizeof(GP)/sizeof(Real));

			//run gold
			timer.start();
			for(int i=0; i<NTIMES; i++)
				refkernel.compute(&blockgold[i%NBLOCKS](0,0,0).s.r, sizeof(GP)/sizeof(Real));
			tGOLD = timer.stop();
			
#pragma omp parallel
			{
				const double t = _benchmarkSOS(kernel, NBLOCKS, NTIMES);
#pragma omp critical
				{
					tCOMPUTE += t;
				}				
			}
			tCOMPUTE /= COUNT;
		}
		
		TSOS::printflops(PEAKPERF, PEAKBAND, 1, 1, NTIMES, tCOMPUTE, false);
		cout << setprecision(4) << "\tMEMORY FOOTPRINT: "<< (block->kB())/1024 << " MB" << endl;
		printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
		if (bAwk) awkShortPerf("MaxSOS_SSE", PEAKPERF, PEAKBAND, GFLOP, tGOLD, tCOMPUTE, EPERF, TEXPECTED, NTIMES, block->kB(), OI);
		printEndLine();
		
		delete [] block;
		delete [] blockgold;
	}
	
	template<typename TUPDATE>
	void performance(TUPDATE& kernel, Update_CPP& refkernel, const double PEAKPERF = 2.66*8/(sizeof(Real)/4)*1e9, const double PEAKBAND = 4.5*1e9, const int NBLOCKS=8*8*8, const int NTIMES=100, bool bAwk=false)
	{
		Block * blockgold = new Block[NBLOCKS];
		Block * block = new Block[NBLOCKS];
		
		const double GFLOP = _BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*2.*6.*NTIMES/1.e9;
		const double OI = 2./8.;
		const double EPERF = min(OI*PEAKBAND,PEAKPERF);
		const double TEXPECTED = 1.e9*GFLOP/EPERF;
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(blockgold[i]);
		
		for(int i=0; i<NBLOCKS; i++)
			_initialize(block[i]);
		
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
			if (NBLOCKS==1)
				refkernel.compute(&blockgold[0](0,0,0).dsdt.r, &blockgold[0](0,0,0).s.r, sizeof(GP)/sizeof(Real));
			
			timer.start();
			for(int i=0; i<NTIMES; i++)
				refkernel.compute(&blockgold[i%NBLOCKS](0,0,0).dsdt.r, &blockgold[i%NBLOCKS](0,0,0).s.r, sizeof(GP)/sizeof(Real));
			tGOLD = timer.stop();
			
#pragma omp parallel
			{
				const double t = _benchmarkUP(kernel, NBLOCKS, NTIMES);
#pragma omp critical
				{
					tCOMPUTE += t;
				}				
			}
			tCOMPUTE /= COUNT;
		}
		
		TUPDATE::printflops(PEAKPERF, PEAKBAND, 1, 1, NTIMES, tCOMPUTE, false);
		cout << setprecision(4) << "\tMEMORY FOOTPRINT: "<< (block->kB())/1024<< " MB" << endl;

		printf("\tGAIN-OVER-GOLD: %.2fX\n", tGOLD/tCOMPUTE);
		if (bAwk) awkShortPerf("Update_SSE", PEAKPERF, PEAKBAND, GFLOP, tGOLD, tCOMPUTE, EPERF, TEXPECTED, NTIMES, block->kB(), OI);
		printEndLine();
		
		delete [] block;
		delete [] blockgold;
	}
};
