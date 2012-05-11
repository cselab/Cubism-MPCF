/*
 *  output.h
 *  MPCFcore
 *
 *  Created by Christian Conti on 1/12/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <iostream>
#include <iomanip>

using namespace std;

const int widthL = 100;
const int widthS = 80;

static void printKernelName(string kernelname)
{
	//cout << endl;
	//cout << setw(widthL) << setfill('=');
	//cout << endl << endl;
	cout << endl;
	for (int i=0; i<100; i++)
		cout << "=";
	cout << endl << endl;
	cout << "KERNEL - " << kernelname.c_str() << endl << endl;
}

static void printEndKernelTest()
{
	cout << endl;
	for (int i=0; i<100; i++)
		cout << "=";
	cout << endl << endl;
}

static void printAccuracyTitle()
{
	string acc = " ACCURACY ";
	int l = (80-acc.size())/2;
	cout << "\t";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << acc.c_str();
	cout << "";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << endl;
}

static void printPerformanceTitle()
{
	string perf = " PERFORMANCE ";
	int l = (80-perf.size())/2;
	cout << "\t";
	for (int i=0; i<l; i++)
		cout << "-";
	cout << perf.c_str();
	cout << "";
	for (int i=0; i<80-l-(int)perf.size(); i++)
		cout << "-";
	cout << endl;
}

static void printEndLine()
{
	cout << "\t";
	for (int i=0; i<80; i++)
		cout << "-";
	cout << endl;
}

static void awkAcc(string kernelname,
				   double linf0, double linf1, double linf2, double linf3, double linf4, double linf5,
				   double l10, double l11, double l12, double l13, double l14, double l15)
{
	// should also output name
	cout << "awkAccLinf\t" <<  kernelname.c_str();
	cout << "\t" << setprecision(4) << linf0;
	cout << "\t" << setprecision(4) << linf1;
	cout << "\t" << setprecision(4) << linf2;
	cout << "\t" << setprecision(4) << linf3;
	cout << "\t" << setprecision(4) << linf4;
	cout << "\t" << setprecision(4) << linf5;
	cout << endl;
	
	cout << "awkAccL1\t" <<  kernelname.c_str();
	cout << "\t" << setprecision(4) << l10;
	cout << "\t" << setprecision(4) << l11;
	cout << "\t" << setprecision(4) << l12;
	cout << "\t" << setprecision(4) << l13;
	cout << "\t" << setprecision(4) << l14;
	cout << "\t" << setprecision(4) << l15;
	cout << endl;
}

static void awkFSPredictions(string kernelname,
							 const double OIConvert,
							 const double OIWENO,
							 const double OIExtraTerm,
							 const double OICharVel,
							 const double OIHLLERho,
							 const double OIHLLEVel,
							 const double OIHLLEPVel,
							 const double OIHLLEE,
							 const double OIPRHS,
							 const double OICopyBack,
							 const double EPERFCONVERT,
							 const double EPERFWENO,
							 const double EPERFXTRATERM,
							 const double EPERFCHARVEL,
							 const double EPERFHLLERHO,
							 const double EPERFHLLEVEL,
							 const double EPERFHLLEPVEL,
							 const double EPERFHLLEE,
							 const double EPERFPRHS,
							 const double EPERFCOPYBACK)
{
	cout << "awkPredictions\t" << kernelname.c_str();
	cout << "\t" << setprecision(4) << OIConvert;
	cout << "\t" << setprecision(4) << OIWENO;
	cout << "\t" << setprecision(4) << OIExtraTerm;
	cout << "\t" << setprecision(4) << OICharVel;
	cout << "\t" << setprecision(4) << OIHLLERho;
	cout << "\t" << setprecision(4) << OIHLLEVel;
	cout << "\t" << setprecision(4) << OIHLLEPVel;
	cout << "\t" << setprecision(4) << OIHLLEE;
	cout << "\t" << setprecision(4) << OIPRHS;
	cout << "\t" << setprecision(4) << OICopyBack;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFCONVERT;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFWENO;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFXTRATERM;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFCHARVEL;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFHLLERHO;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFHLLEVEL;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFHLLEPVEL;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFHLLEE;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFPRHS;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFCOPYBACK;
	cout << endl;
}

static void awkMCore(string kernelname,
					 const double TOTGFLOP,
					 const double PEAKPERF_CORE,
					 const double PEAKPERF,
					 const double PEAKBAND,
					 const double MEASUREDTIME,
					 const double TEXPECTED,
					 const int NBLOCKS,
					 /*const int NCORES,*/
					 const int NT)
{
	cout << "awkMCore\t" << kernelname.c_str();
	cout << "\t" << NT;
	cout << "\t" << setprecision(4) << TOTGFLOP;
	cout << "\t" << setprecision(4) << PEAKPERF_CORE*1e-9;
	cout << "\t" << setprecision(4) << PEAKPERF*1e-9;
	cout << "\t" << setprecision(4) << PEAKBAND*1e-9;
	cout << "\t" << setprecision(4) << PEAKPERF/PEAKBAND;
	cout << "\t" << setprecision(4) << TOTGFLOP/MEASUREDTIME;
	cout << "\t" << setprecision(4) << TOTGFLOP*1e9/TEXPECTED/PEAKBAND;
	cout << "\t" << setprecision(4) << 1e3*MEASUREDTIME/NBLOCKS;
	cout << "\t" << setprecision(4) << 1e3*TEXPECTED/NBLOCKS;
	cout << "\t" << setprecision(4) << TOTGFLOP/TEXPECTED;
	cout << "\t" << setprecision(4) << 100.*TEXPECTED/MEASUREDTIME;
	cout << "\t" << setprecision(4) << 100*(TOTGFLOP/MEASUREDTIME*1e9)/PEAKPERF;
	cout << "\t" << setprecision(4) << MEASUREDTIME;
	cout << endl;
}

static void awkMCorePredictions(string kernelname,
								const double OI,
								const double EPERF)
{
	cout << "awkMCorePredictions\t" << kernelname.c_str();
	cout << "\t" << setprecision(4) << OI;
	cout << "\t" << setprecision(4) << 1.e-9*EPERF;
	cout << endl;
}

static void awkMCorePredictions(string kernelname,
								const double OIConvert,
								const double OIWENO,
								const double OIExtraTerm,
								const double OICharVel,
								const double OIHLLERho,
								const double OIHLLEVel,
								const double OIHLLEPVel,
								const double OIHLLEE,
								const double OIPRHS,
								const double OICopyBack,
								const double OITotal,
								const double EPERFCONVERT,
								const double EPERFWENO,
								const double EPERFXTRATERM,
								const double EPERFCHARVEL,
								const double EPERFHLLERHO,
								const double EPERFHLLEVEL,
								const double EPERFHLLEPVEL,
								const double EPERFHLLEE,
								const double EPERFPRHS,
								const double EPERFCOPYBACK,
								const double EPERFTotal,
								const double AIEPERFTotal)
{
	cout << "awkMCorePredictions\t" << kernelname.c_str();
	cout << "\t" << setprecision(4) << OIConvert;
	cout << "\t" << setprecision(4) << OIWENO;
	cout << "\t" << setprecision(4) << OIExtraTerm;
	cout << "\t" << setprecision(4) << OICharVel;
	cout << "\t" << setprecision(4) << OIHLLERho;
	cout << "\t" << setprecision(4) << OIHLLEVel;
	cout << "\t" << setprecision(4) << OIHLLEPVel;
	cout << "\t" << setprecision(4) << OIHLLEE;
	cout << "\t" << setprecision(4) << OIPRHS;
	cout << "\t" << setprecision(4) << OICopyBack;
	cout << "\t" << setprecision(4) << OITotal;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFCONVERT;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFWENO;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFXTRATERM;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFCHARVEL;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFHLLERHO;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFHLLEVEL;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFHLLEPVEL;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFHLLEE;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFPRHS;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFCOPYBACK;
	cout << "\t" << setprecision(4) << 1.e-9*EPERFTotal;
	cout << "\t" << setprecision(4) << 1.e-9*AIEPERFTotal;
	cout << endl;
}

//#ifdef _SSE_
static void awkShortPerf(string kernelname,
						 double PEAKPERF,
						 double PEAKBAND,
						 double GFLOP,
						 double tGOLD,
						 double tCOMPUTE,
						 double EPERF,
						 double tEXPECTED,
						 int NTIMES,
						 double kB,
						 double OI)
{
	// should also output name
	cout << "awkShortPerf\t" <<  kernelname.c_str();
	cout << "\t" << setprecision(4) << kB/1024.;
	cout << "\t" << setprecision(4) << PEAKPERF*1e-9;
	cout << "\t" << setprecision(4) << PEAKBAND*1e-9;
	cout << "\t" << setprecision(4) << PEAKPERF/PEAKBAND;
	cout << "\t" << setprecision(4) << GFLOP/tGOLD;
	cout << "\t" << setprecision(4) << GFLOP/tCOMPUTE;
	cout << "\t" << setprecision(4) << EPERF/1.e9;
	cout << "\t" << setprecision(4) << tGOLD/tCOMPUTE;
	cout << "\t" << setprecision(4) << tCOMPUTE*1e3;
	cout << "\t" << setprecision(4) << tEXPECTED*1e3;
	cout << "\t" << setprecision(4) << 1e3*tCOMPUTE/NTIMES;
	cout << "\t" << setprecision(4) << 1e3*tEXPECTED/NTIMES;
	cout << "\t" << setprecision(4) << 100*tEXPECTED/tCOMPUTE;
	cout << "\t" << setprecision(4) << 100*(GFLOP/tCOMPUTE*1e9)/PEAKPERF;
	cout << "\t" << setprecision(4) << OI;
	cout << endl;
}

static void awkPerf(string kernelname,
					double PEAKPERF,
					double PEAKBAND,
					double tGOLD,
					double tCOMPUTE,
					double GFLOP,
					double TEXPECTED,
					int NTIMES,
					double labkB,
					double blockkB,
					double fskB)
{
	// should also output name
	cout << "awkPerf\t" <<  kernelname.c_str();
	cout << "\t" << setprecision(4) << (labkB +  blockkB + fskB)/1024;
	cout << "\t" << setprecision(4) << PEAKPERF*1e-9;
	cout << "\t" << setprecision(4) << PEAKBAND*1e-9;
	cout << "\t" << setprecision(4) << PEAKPERF/PEAKBAND;
	cout << "\t" << setprecision(4) << GFLOP/tGOLD;
	cout << "\t" << setprecision(4) << GFLOP/tCOMPUTE;
	cout << "\t" << setprecision(4) << tGOLD/tCOMPUTE;
	cout << "\t" << setprecision(4) << 1e3*tCOMPUTE/NTIMES;
	cout << "\t" << setprecision(4) << 1e3*TEXPECTED/NTIMES;
	cout << "\t" << setprecision(4) << 100.*TEXPECTED/tCOMPUTE;
	cout << "\t" << setprecision(4) << 100.*(GFLOP/tCOMPUTE*1e9)/PEAKPERF;
	cout << endl;
}
/*#else
static void awkShortPerf(double PEAKPERF,
						 double PEAKBAND,
                         double GFLOP,
						 double tGOLD,
						 double EPERF,
						 double tEXPECTED,
                         int NTIMES,
						 double kB)
{
	// should also output name
	cout << "awkShortPerf";
	cout << "\t" << setprecision(4) << kB/1024.;
	cout << "\t" << setprecision(4) << PEAKPERF*1e-9;
	cout << "\t" << setprecision(4) << PEAKBAND*1e-9;
	cout << "\t" << setprecision(4) << PEAKPERF/PEAKBAND;
	cout << "\t" << setprecision(4) << GFLOP/tGOLD;
	cout << "\t" << setprecision(4) << EPERF/1.e9;
	cout << "\t" << setprecision(4) << tGOLD*1e3;
	cout << "\t" << setprecision(4) << tEXPECTED*1e3;
	cout << "\t" << setprecision(4) << 1e3*tGOLD/NTIMES;
	cout << "\t" << setprecision(4) << 1e3*tEXPECTED/NTIMES;
	cout << "\t" << setprecision(4) << 100*tEXPECTED/tGOLD;
	cout << "\t" << setprecision(4) << 100*(GFLOP/tGOLD*1e9)/PEAKPERF;
	cout << endl;
}

void awkPerf(double PEAKPERF,
			 double PEAKBAND,
			 double tGOLD,
			 double tCOMPUTE,
			 double GFLOP,
			 double TEXPECTED,
			 int NTIMES,
			 double labkB,
			 double blockkB,
			 double fskB)
{
	// should also output name
	cout << "awkPerf";
	cout << "\t" << setprecision(4) << (labkB +  blockkB + fskB)/1024;
	cout << "\t" << setprecision(4) << PEAKPERF*1e-9;
	cout << "\t" << setprecision(4) << PEAKBAND*1e-9;
	cout << "\t" << setprecision(4) << PEAKPERF/PEAKBAND;
	cout << "\t" << setprecision(4) << GFLOP/tGOLD;
	cout << "\t" << setprecision(4) << GFLOP/tCOMPUTE;
	cout << "\t" << setprecision(4) << tGOLD/tCOMPUTE;
	cout << "\t" << setprecision(4) << 1e3*tCOMPUTE/NTIMES;
	cout << "\t" << setprecision(4) << 1e3*TEXPECTED/NTIMES;
	cout << "\t" << setprecision(4) << 100.*TEXPECTED/tCOMPUTE;
	cout << "\t" << setprecision(4) << 100.*(GFLOP/tCOMPUTE*1e9)/PEAKPERF;
	cout << endl;
}
#endif*/
