/*
 *  Flowstep_CPPfloat.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/6/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#define _ALIGNBYTES 16

#include <cassert>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <typeinfo>
#include "output.h"

using namespace std;

#ifndef WENOEPS
#define WENOEPS 1.0e-6
#endif

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

template < int _SX, int _EX, int _SY, int _EY, typename TReal=Real > 
struct SOA2D
{
	static const int _CPERALIGNBYTES = _ALIGNBYTES_/sizeof(TReal);
	
	static const int SX = _CPERALIGNBYTES*((_SX - (_CPERALIGNBYTES-1))/_CPERALIGNBYTES);
	static const int EX = _CPERALIGNBYTES*((_EX + (_CPERALIGNBYTES-1))/_CPERALIGNBYTES);
	static const int NX = _EX - _SX;
	static const int NY = _EY - _SY;
	
	static const int PITCH = EX - SX;
	
	TReal __attribute__((aligned(_ALIGNBYTES_))) data[NY][PITCH];
	
	SOA2D()
	{
		assert(((size_t)(&data[0][0]) % _ALIGNBYTES_) == 0);
	}
	
	inline TReal operator()(const int ix, const int iy) const
	{
		assert(ix >= _SX); assert(ix < _EX);
		assert(iy >= _SY); assert(iy < _EY);
		
		return data[iy-_SY][ix-SX];
	}
	
	inline const TReal * ptr(const int ix, const int iy) const
	{
		assert(ix >= SX); assert(ix < EX);
		assert(iy >= _SY); assert(iy < _EY);
		
		return &data[iy-_SY][ix-SX];
	}
	
	inline TReal& ref(const int ix, const int iy)
	{
		assert(ix >= SX); assert(ix < EX);
		assert(iy >= _SY); assert(iy < _EY);
		
		return data[iy-_SY][ix-SX];
	}
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(SOA2D)/1024.;
	}
	
	inline SOA2D& operator= (const SOA2D& c)
	{
		for(int y=0; y<NY; ++y)
			for(int x=0; x<PITCH; ++x)
				data[y][x] = c.data[y][x];
		
		return *this;
	}
};

template< int _SX, int _EX, int _SY, int _EY, int _NSLICES, typename TReal=Real>
struct RingSOA2D
{
	int currslice;
	
	SOA2D<_SX, _EX, _SY, _EY, TReal> slices[_NSLICES];
	
	RingSOA2D(): currslice(0){ }
	
	inline const SOA2D<_SX, _EX, _SY, _EY, TReal>& operator()(const int relativeid=0) const
	{
		return slices[(relativeid + currslice + _NSLICES) % _NSLICES];
	}
	
	inline SOA2D<_SX, _EX, _SY, _EY, TReal>& ref(const int relativeid=0)
	{
		return slices[(relativeid + currslice + _NSLICES) % _NSLICES];
	}
	
	void next(){ currslice = (currslice + 1) % _NSLICES; } 
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(RingSOA2D)/1024.;
	}
};

typedef RingSOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3, 6> RingInputSOA;
typedef SOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3> InputSOA;
typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_> TempSOA;
typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_, 2> RingTempSOA;
typedef RingSOA2D<0,_BLOCKSIZE_+1, 0, _BLOCKSIZE_, 3> RingTempSOA3;
typedef SOA2D<0, _BLOCKSIZE_, 0, _BLOCKSIZE_> OutputSOA;

template< class SOAtype > void _check(const SOAtype& a, const SOAtype& b)
{
	for(int iy=0; iy<SOAtype::NY; iy++)
        for(int ix=0; ix<SOAtype::NX; ix++)
		{
			const double va = a(ix, iy);
			const double vb = b(ix, iy);
			
			const double err = va-vb;
			const double tol = 1e-5;
			const double relerr = err/max(max(tol, va), vb);
			
			if (relerr >= tol && err > tol)
				printf("ix:%d iy:%d a:%.15f b:%.15f\t\ttol:%e err:%e rel.err:%e \n", ix, iy, a(ix, iy), b(ix, iy), tol, err, relerr);
			
			assert(relerr < tol || err < tol);
		}
}

class FlowStep_CPP
{	
protected:
	
	//soa input
	RingInputSOA ringrho, ringu, ringv, ringw, ringp, ringls;
	RingTempSOA wenorho, wenou, wenov, wenop, charvel;
	RingTempSOA3 wenow, wenols;
	RingTempSOA fluxrho, fluxu, fluxv, fluxw, fluxp, fluxls;
	OutputSOA sumls, divu;
	OutputSOA rhsrho, rhsu, rhsv, rhsw, rhss, rhsls;
	
	Real m_gamma1, m_gamma2, m_smoothlength;
	Real m_a, m_dtinvh;
    Real m_pc1, m_pc2;
	
	//16 FLOP/call, 4 Real/call, rsingle=1.125
	Real _getgamma(const Real phi);
	Real _getPC(const Real phi);
    
	//63 FLOP/pt, 12 Real/pt, rsingle=0.5
	void _convert(const Real * const gptfirst, const int gptfloats, const int rowgpts);
	
	//82 FLOP/pt, 6 Real/pt, rsingle=2.5
	void _xweno_minus(const InputSOA& in, TempSOA& out);
	void _xweno_pluss(const InputSOA& in, TempSOA& out);
	void _yweno_minus(const InputSOA& in, TempSOA& out);
	void _yweno_pluss(const InputSOA& in, TempSOA& out);
	void _zweno_minus(const int relid, const RingInputSOA& in, TempSOA& out);
	void _zweno_pluss(const int relid, const RingInputSOA& in, TempSOA& out);
	
	//2 FLOP/pt, 6 Real/pt, rsingle=0.084
	void _xextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp);
	
	//4 FLOP/pt, 8 Real/pt, rsingle=0.125
	void _yextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp);
	void _zextraterm(const TempSOA& um, const TempSOA& up, const TempSOA& lm, const TempSOA& lp);
	
	//80 FLOP/pt, 10 Real/pt, rsingle=1.3
	void _char_vel(const TempSOA& rminus, const TempSOA& rplus, 
				   const TempSOA& vminus, const TempSOA& vplus,
				   const TempSOA& pminus, const TempSOA& pplus,
				   const TempSOA& lminus, const TempSOA& lplus, 
				   TempSOA& out_minus, TempSOA& out_plus);
	
	//17 FLOP/pt, 7 Real/pt, rsingle=0.571
	void _hlle_rho(const TempSOA& rm, const TempSOA& rp,
				   const TempSOA& vm, const TempSOA& vp,
				   const TempSOA& am, const TempSOA& ap,
				   TempSOA& out);
	
	//19 FLOP/pt, 9 Real/pt, rsingle=0.25
	void _hlle_vel(const TempSOA& rminus, const TempSOA& rplus,
				   const TempSOA& vminus, const TempSOA& vplus,
				   const TempSOA& vdminus, const TempSOA& vdplus,
				   const TempSOA& aminus, const TempSOA& aplus,
				   TempSOA& out);
	
	//21 FLOP/pt, 9 Real/pt, rsingle=0.584
	void _hlle_pvel(const TempSOA& rminus, const TempSOA& rplus,
					const TempSOA& vminus, const TempSOA& vplus,
					const TempSOA& pminus, const TempSOA& pplus,
					const TempSOA& aminus, const TempSOA& aplus,
					TempSOA& out);
	//73 FLOP/pt, 15 Real/pt, rsingle=1.283
	void _hlle_e(const TempSOA& rminus, const TempSOA& rplus,
				 const TempSOA& vdminus, const TempSOA& vdplus,
				 const TempSOA& v1minus, const TempSOA& v1plus,
				 const TempSOA& v2minus, const TempSOA& v2plus,
				 const TempSOA& pminus, const TempSOA& pplus,
				 const TempSOA& lsminus, const TempSOA& lsplus, 
				 const TempSOA& aminus, const TempSOA& aplus,
				 TempSOA& out);
	
	//1 FLOP/pt, 3 Real/pt, rsingle=0.08
	void _xrhsadd(const TempSOA& flux, OutputSOA& rhs);
	
	//2 FLOP/pt, 4 Real/pt, rsingle=0.125
	void _yrhsadd(const TempSOA& flux, OutputSOA& rhs);
	void _zrhsadd(const TempSOA& fback, const TempSOA& fforward, OutputSOA& rhs);
	
	//21 FLOP/pt, 14 Real/pt, rsingle=0.161
	void _copyback(Real * const gptfirst, const int gptfloats, const int rowgpts);
	
	//911 FLOP/pt, 156 Real/pt, rsingle=1.45 (unfair average!)
	void _xflux(const int relsliceid);
	void _yflux(const int relsliceid);
	void _zflux(const int relsliceid);
	
	//6 FLOP/pt, 18 Real/pt, rsingle=0.083
	void _xrhs();
	
	//12 FLOP/pt, 24 Real/pt, rsingle=0.125
	void _yrhs();
	void _zrhs();
	
	void _next() 
	{ 
		ringrho.next(); ringu.next(); ringv.next(); ringw.next(); ringp.next(); ringls.next();
	}
	
	void _flux_next()
	{
		fluxrho.next(); fluxu.next(); fluxv.next(); fluxw.next(); fluxp.next(); fluxls.next();
		
		wenow.next(); wenols.next();
	}
	
public:
	
	FlowStep_CPP(Real a=0, Real dtinvh=1, Real gamma1=2.5, Real gamma2=2.1, Real smoothlength=1, Real pc1=0, Real pc2=0):
	m_a(a), m_dtinvh(dtinvh), m_gamma1(gamma1), m_gamma2(gamma2), m_smoothlength(smoothlength), m_pc1(pc1), m_pc2(pc2)
	{}
	
	void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
				 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts);
	
	float kB(const bool bVerbos=true) const;
	
	Real gamma1() const { return m_gamma1; }
	Real gamma2() const { return m_gamma2; }
	Real smoothlength() const { return m_smoothlength; }
	Real dtinvh() const { return m_dtinvh; }	
	
	static double getGFLOP(const int NBLOCKS=1)
	{
		const double FLOPCONVERT = 63*(_BLOCKSIZE_+6)*(_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);
		const double FLOPFLUX = 1232*_BLOCKSIZE_*(_BLOCKSIZE_+1);
		const double FLOPRHS = (12+12+6)*_BLOCKSIZE_*_BLOCKSIZE_;
		const double FLOPCOPYBACK = 21*_BLOCKSIZE_*_BLOCKSIZE_;
		return NBLOCKS*1.*(FLOPCONVERT + FLOPFLUX + _BLOCKSIZE_*(FLOPFLUX*3 + FLOPRHS + FLOPCOPYBACK))/1e9;
	}
	
	static double getTEXPECTED(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES=1, const int NBLOCKS=1)
	{
		const float PEAKPERF = PEAKPERF_CORE*NCORES;
		const float OITotal     = 1e9*(getGFLOP(NBLOCKS)/NBLOCKS)/((_BLOCKSIZE_+6)*(_BLOCKSIZE_+6)*(_BLOCKSIZE_+6)*sizeof(Real)*6 + 2*_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*sizeof(Real)*6);
		const double EPERFTotal    = min(OITotal*    PEAKBAND, PEAKPERF);

		return getGFLOP(NBLOCKS)*1e9/EPERFTotal;
	}
	
	static void printflops(const float PEAKPERF_CORE, const float PEAKBAND, const int NCORES, const int NT, const int NBLOCKS, const float MEASUREDTIME, const bool bAwk=false, const string implname="FS")
	{
		const float PEAKPERF = PEAKPERF_CORE*NCORES;
		
		//FLOP estimation
		const double FLOPCONVERT = 63*(_BLOCKSIZE_+6)*(_BLOCKSIZE_+6)*(_BLOCKSIZE_+6);
		const double FLOPRHS = (12+12+6)*_BLOCKSIZE_*_BLOCKSIZE_;
		const double FLOPCOPYBACK = 21*_BLOCKSIZE_*_BLOCKSIZE_;
		const double GFLOP = getGFLOP(NBLOCKS);
		
		//FLOP/s estimation
		// for OI use GFLOP/block, these are Arithmetic Intensities
		const float AIConvert   = .5f;
		const float AIWENO      = 2.5f;
		const float AIExtraTerm = .125f;
		const float AICharVel   = 1.3f;
		const float AIHLLERho   = .571f;
		const float AIHLLEVel   = .25f;
		const float AIHLLEPVel  = .584f;
		const float AIHLLEE     = 1.283f;
		const float AIPRHS      = .125f;
		const float AICopyBack  = .161f;
		const float OITotal     = 1e9*(GFLOP/NBLOCKS)/((_BLOCKSIZE_+6)*(_BLOCKSIZE_+6)*(_BLOCKSIZE_+6)*sizeof(Real)*6 + 2*_BLOCKSIZE_*_BLOCKSIZE_*_BLOCKSIZE_*sizeof(Real)*6);
		//printf("OITOTAL IS %f\n", OITotal); exit(0);

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
		const double EPERFTotal    = min(OITotal*    PEAKBAND, PEAKPERF);
		
	/*	cout << "EPERFCONVERT: " << EPERFCONVERT*1e-9 << endl;
		cout << "EPERWENO: " << EPERWENO*1e-9 << endl;
		cout << "EPERXTRATERM: " << EPERXTRATERM*1e-9 << endl;
		cout << "EPERCHARVEL: " << EPERCHARVEL*1e-9 << endl;
		cout << "EPERHLLERHO: " << EPERHLLERHO*1e-9 << endl;
		cout << "EPERHLLEVEL: " << EPERHLLEVEL*1e-9 << endl;
		cout << "EPERHLLEPVEL: " << EPERHLLEPVEL*1e-9 << endl;
		cout << "EPERHLLEE: " << EPERHLLEE*1e-9 << endl;
		cout << "EPERFPRHS: " << EPERFPRHS*1e-9 << endl;
		cout << "EPERFCOPYBACK: " << EPERFCOPYBACK*1e-9 << endl;
	*/	
		//execution time estimation
		const double TCONVERT= FLOPCONVERT/EPERFCONVERT;
		const double TFLUX = _BLOCKSIZE_*(_BLOCKSIZE_+1)*(60*12./EPERWENO + 8./EPERXTRATERM + 52./EPERCHARVEL + 16*2./EPERHLLERHO + 21./EPERHLLEPVEL + 2*14./EPERHLLEVEL + 77./EPERHLLEE);//GFLOPFLUX/EPERFFLUX;
		const double TRHS = FLOPRHS/EPERFPRHS;
		const double TCOPYBACK = FLOPCOPYBACK/EPERFCOPYBACK;
		const double TEXPECTED = GFLOP*1e9/EPERFTotal;
		const double TEXPECTEDAI = NBLOCKS*1.*(TCONVERT + TFLUX + _BLOCKSIZE_*(TFLUX*3 + TRHS + TCOPYBACK));
		
		//printf("=================  PERFORMANCE   =================\n");
		printPerformanceTitle();
		//cout << setprecision(4) << "MEMORY FOOTPRINT: "<< (lab->kB() +  block->kB() + fs.kB(false))/1024<< " MB" << endl;
		cout << "\tFS TOTAL GFLOPS " << GFLOP <<endl ;
		printf("\tFS ASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tFS RIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tFS THIS ONE IS %.2f GFLOP/s,\t\"PER FS\" %.2f FLOP/B (OI),\t\"PER FS\" %.2f FLOP/B (AI)\n", GFLOP/MEASUREDTIME, OITotal,GFLOP*1e9/TEXPECTEDAI/PEAKBAND);
		printf("\tFS TIME PER BLOCK: %.5f ms (expected %.5f ms [OI],\texpected %.5f ms [AI])\n",  1e3*MEASUREDTIME/NBLOCKS, 1e3*TEXPECTED/NBLOCKS, 1e3*TEXPECTEDAI/NBLOCKS);
		cout << "\tFS Expected Performance is: "<< 1.e-9*EPERFTotal << " GFLOP/s [OI]" << endl;
		cout << "\tFS Expected Performance is: "<< GFLOP/TEXPECTEDAI << " GFLOP/s [AI]" << endl;
		printf("\tFS EFFICIENCY: %.2f%%, HW-UTILIZATION: %.2f%% [OI]\n", 100.*1.e9*GFLOP/EPERFTotal/MEASUREDTIME, 100*(GFLOP/MEASUREDTIME*1e9)/PEAKPERF);
		printf("\tFS EFFICIENCY: %.2f%%, HW-UTILIZATION: %.2f%% [AI]\n", 100.*TEXPECTEDAI/MEASUREDTIME, 100*(GFLOP/MEASUREDTIME*1e9)/PEAKPERF);
		
		if (bAwk)
		{			
			awkMCorePredictions(implname,
								AIConvert, AIWENO, AIExtraTerm, AICharVel, AIHLLERho, AIHLLEVel, AIHLLEPVel, AIHLLEE, AIPRHS, AICopyBack, OITotal,
								EPERFCONVERT, EPERWENO, EPERXTRATERM, EPERCHARVEL, EPERHLLERHO, EPERHLLEVEL, EPERHLLEPVEL, EPERHLLEE, EPERFPRHS, EPERFCOPYBACK, EPERFTotal, GFLOP/TEXPECTED);
			awkMCore(implname, GFLOP, PEAKPERF_CORE, PEAKPERF, PEAKBAND, MEASUREDTIME, 1.e9*GFLOP/EPERFTotal, NBLOCKS, /*NCORES,*/ NT);
		}
		printEndLine();
		
		{
			static int counter = 0;
			
			std::ofstream report("report.txt", counter==0? ios::out : ios::app);
			
			if(counter == 0)
				report<<"STEPID\tGFLOP/s\tEXPECTED GFLOP/s"<<endl;
			
			report<<counter<<"\t"<<GFLOP/MEASUREDTIME<<"\t"<<GFLOP/TEXPECTED << endl;
			
			counter++;
		}
	}
};
