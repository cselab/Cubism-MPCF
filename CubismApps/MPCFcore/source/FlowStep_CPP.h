/*
 *  Flowstep_CPPfloat.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/6/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cassert>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <typeinfo>

#include "common.h"
#include "SOA2D.h"

using namespace std;

class FlowStep_CPP
{	
protected:

	typedef RingSOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3, 6> RingInputSOA;
	typedef SOA2D<-3, _BLOCKSIZE_+3, -3, _BLOCKSIZE_+3> InputSOA;
	typedef SOA2D<0,_BLOCKSIZE_+1, 0,_BLOCKSIZE_> TempSOA;
	typedef RingSOA2D<0, _BLOCKSIZE_+1, 0,_BLOCKSIZE_, 2> RingTempSOA;
	typedef RingSOA2D<0,_BLOCKSIZE_+1, 0, _BLOCKSIZE_, 3> RingTempSOA3;
	
	struct AssumedType { Real r, u, v, w, s, l; }; //Assumed input/output grid point type
	
	//working dataset
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
		
		//execution time estimation
		const double TCONVERT= FLOPCONVERT/EPERFCONVERT;
		const double TFLUX = _BLOCKSIZE_*(_BLOCKSIZE_+1)*(60*12./EPERWENO + 8./EPERXTRATERM + 52./EPERCHARVEL + 16*2./EPERHLLERHO + 21./EPERHLLEPVEL + 2*14./EPERHLLEVEL + 77./EPERHLLEE);//GFLOPFLUX/EPERFFLUX;
		const double TRHS = FLOPRHS/EPERFPRHS;
		const double TCOPYBACK = FLOPCOPYBACK/EPERFCOPYBACK;
		const double TEXPECTED = GFLOP*1e9/EPERFTotal;
		const double TEXPECTEDAI = NBLOCKS*1.*(TCONVERT + TFLUX + _BLOCKSIZE_*(TFLUX*3 + TRHS + TCOPYBACK));
		
		printPerformanceTitle();

		cout << "\tFS TOTAL GFLOPS " << GFLOP <<endl ;
		printf("\tFS ASSUMING PP: %.2f GFLOP/s (PER CORE), %.2f GFLOP/s (OVERALL)\n\tPB: %.2f GB/s (OVERALL)\n", PEAKPERF_CORE*1e-9, PEAKPERF*1e-9, PEAKBAND*1e-9);
		printf("\tFS RIDGE AT %.2f FLOP/B\n", PEAKPERF/PEAKBAND);
		printf("\tFS THIS ONE IS %.2f GFLOP/s,\t\"PER FS\" %.2f FLOP/B (OI),\t\"PER FS\" %.2f FLOP/B (AI)\n", GFLOP/MEASUREDTIME, OITotal,GFLOP*1e9/TEXPECTEDAI/PEAKBAND);
		printf("\tFS TIME PER BLOCK: %.5f ms (expected %.5f ms [OI],\texpected %.5f ms [AI])\n",  1e3*MEASUREDTIME/NBLOCKS, 1e3*TEXPECTED/NBLOCKS, 1e3*TEXPECTEDAI/NBLOCKS);
		cout << "\tFS Expected Performance is: "<< 1.e-9*EPERFTotal << " GFLOP/s [OI]" << endl;
		cout << "\tFS Expected Performance is: "<< GFLOP/TEXPECTEDAI << " GFLOP/s [AI]" << endl;
		printf("\tFS EFFICIENCY: %.2f%%, HW-UTILIZATION: %.2f%% [OI]\n", 100.*1.e9*GFLOP/EPERFTotal/MEASUREDTIME, 100*(GFLOP/MEASUREDTIME*1e9)/PEAKPERF);
		printf("\tFS EFFICIENCY: %.2f%%, HW-UTILIZATION: %.2f%% [AI]\n", 100.*TEXPECTEDAI/MEASUREDTIME, 100*(GFLOP/MEASUREDTIME*1e9)/PEAKPERF);
		
		printEndLine();
		
		{
			static int counter = 0;
			
			std::ofstream report("report.txt", counter==0? ios::out : ios::app);
			
			if(counter == 0)
				report<<"STEPID\tGFLOP/s\tEXPECTED GFLOP/s"<<endl;
			
			report << counter << "\t"<<GFLOP/MEASUREDTIME<<"\t"<<GFLOP/TEXPECTED << endl;
			
			counter++;
		}
	}
};
