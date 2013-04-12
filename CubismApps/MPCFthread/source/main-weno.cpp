#include <iostream>
#include <limits>
#include <string>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>

#ifndef WENOEPS
#define WENOEPS 1.e-6f
#endif

#ifdef _PROFILE_
#include <mpi.h>
extern "C" void HPM_Start(char *);
extern "C" void HPM_Stop(char *);
#endif



#ifndef _PRECDIV_
static const bool precdiv = false;
#else
static const bool precdiv = true;
#endif

#ifdef __xlC__
extern
#ifdef __cplusplus
"builtin"
#endif
void __alignx (int n, const void *addr);
#endif

#include "common.h"
#include "Weno_CPP.h"
#include "Weno_QPX.h"
#ifdef _WENO3_
#include "Weno_QPX_3rdOrder.h"
#endif
Real * myalloc(const int NENTRIES,  const bool verbose )
{
	const bool initialize = true;
	
	enum { alignment_bytes = 32 } ;
	
	Real * tmp = NULL;
	
	const int result = posix_memalign((void **)&tmp, alignment_bytes, sizeof(Real) * NENTRIES);
	assert(result == 0);
	
	if (initialize)
	{
		for(int i=0; i<NENTRIES; ++i)
			tmp[i] = drand48();
		
		if (verbose)
		{
			for(int i=0; i<NENTRIES; ++i)
				printf("tmp[%d] = %f\n", i, tmp[i]);
			printf("==============\n");
		}
	}
	
	return tmp;
}

class Timer
{
	struct timeval t_start, t_end;
	struct timezone t_zone; 
	
public:
	
	void start()
	{
		gettimeofday(&t_start,  &t_zone);
	}
	
	double stop()
	{		
		gettimeofday(&t_end,  &t_zone);
		return (t_end.tv_usec  - t_start.tv_usec)*1e-6  + (t_end.tv_sec  - t_start.tv_sec);
	}
};

template<typename T>
void check_error(const double tol, T ref[], T val[], const int N)
{
	static const bool verbose = false;
	
	for(int i=0; i<N; ++i)
	{
		assert(!std::isnan(ref[i]));
		assert(!std::isnan(val[i]));
		
		const double err = ref[i] - val[i];
		const double relerr = err/std::max((double)std::numeric_limits<T>::epsilon(), (double)std::max(fabs(val[i]), fabs(ref[i])));
		
		if (verbose) printf("+%1.1e,", relerr);
		
		if (fabs(relerr) >= tol && fabs(err) >= tol)
			printf("\n%d: %e %e -> %e %e\n", i, ref[i], val[i], err, relerr);
		
		assert(fabs(relerr) < tol || fabs(err) < tol);
	}
	
	if (verbose) printf("\t");
}

template<typename T>
void check_error(const double tol, T ref, T val)
{
	check_error(tol, &ref, &val, 1);
}

void print_performance(const int NENTRIES, const double t, const std::string title, const bool printheader)
{
	printf("============  %s  ================\n", title.c_str());
#ifdef _WENO3_
	const double FLOP_PER_ITERATION = (19 + 3 * 4) ;
#else
	const double FLOP_PER_ITERATION = 78;
#endif
	const double GFLOP_RATE = FLOP_PER_ITERATION * 1e-9 * NENTRIES / t; 
	const double AI = FLOP_PER_ITERATION / sizeof(Real) / 6.;
	const double PP = 3.2 * 4;
	const double PB = 0.67 * 4;
	const double expected = std::min(PP, PB * AI);
	
	if (printheader)
		printf("AI: %.2f FLOP/B -> given PP = %.2f and PB = %.2f we expect %.2f GFLOP/s\n", AI, PP, PB, expected);
	
	printf("PERFORMANCE: %.2f GFLOP/s TIME: %.3f ms\n", GFLOP_RATE, 1e3 * t);
	printf("We are on %.2f %% of the max achievable performance\n", GFLOP_RATE / expected * 100);
	
	if (printheader) 
		printf("By the way the ridge point is at %.2f FLOP/BYTE\n", PP/PB);
}

void benchmark(int argc, char *  argv[], const int NENTRIES_, const int NTIMES, const bool verbose, std::string benchmark_name)
{
	const int NENTRIES = 4 * (NENTRIES_ / 4);
	
	printf("nentries set to %e\n", (float)NENTRIES);
	
	Timer timer;
	
	Real * const a = myalloc(NENTRIES, verbose);
	Real * const b = myalloc(NENTRIES, verbose);
	Real * const c = myalloc(NENTRIES, verbose);
	Real * const d = myalloc(NENTRIES, verbose);
	Real * const e = myalloc(NENTRIES, verbose);
	Real * const f = myalloc(NENTRIES, verbose);
	Real * const gold = myalloc(NENTRIES, verbose);
	Real * const goldp = myalloc(NENTRIES, verbose);
	Real * const result = myalloc(NENTRIES, verbose);
	Real * const resultp = myalloc(NENTRIES, verbose);
	
	WenoReference<Weno_FunctorMinus> weno_reference;
	WenoReference<Weno_FunctorPlus> weno_referencep;
	weno_reference(a, b, c, d, e, gold, NENTRIES);
	weno_referencep(b, c, d, e, f, goldp, NENTRIES);
	
	Weno_QPX_fused asd;
	//WenoQPX_3rdOrder_Minus asd;
	
	for(int i=0 ;i<NENTRIES; i+= 4)
	{
	//, resultm + i, resultp + i
	//	vec_sta(asd._weno(vec_lda(0L, a + i), vec_lda(0L, b + i), vec_lda(0L, c + i), 
	//			  vec_lda(0L, d + i), vec_lda(0L, e + i)), 0L, result + i);
		asd.weno_minus_plus_fused_opt2(vec_lda(0L, a + i), vec_lda(0L, b + i), vec_lda(0L, c + i), 
									   vec_lda(0L, d + i), vec_lda(0L, e + i), vec_lda(0L, f + i), result + i, resultp + i);
	}
	
	printf("finito. \n");
	
	/* lets just touch the results */
	//weno_reference(a, b, c, d, e, gold, NENTRIES);
	/*
	timer.start();
#pragma omp parallel
	{
		for(int i=0; i<NTIMES; ++i)
			weno_reference(a, b, c, d, e, gold, NENTRIES);
	}
	const double tref = timer.stop();
	
	if (verbose)
		for(int i=0; i<NENTRIES; ++i)
			printf("gold[%d] = %f\n", i, gold[i]);
	
	/* lets just touch the results * /
	weno_reference(a, b, c, d, e, result, NENTRIES);
	
#ifdef _WENO3_
	WenoQPX<WenoQPX_3rdOrder_Minus> wenoqpx;
#else
	WenoQPX<WenoQPX_MinusFunctor> wenoqpx;
#endif
	
	timer.start();
#pragma omp parallel
	{
		for(int i=0; i<NTIMES; ++i)
			wenoqpx(a, b, c, d, e, result, NENTRIES);
	}
	const double t = timer.stop();
	
	/* lets just touch the results * /
	weno_reference(a, b, c, d, e, result, NENTRIES);
	
#ifdef _WENO3_
	WenoQPXUnrolled<WenoQPX_3rdOrder_Minus> wenoqpx_unrolled;
#else
	WenoQPXUnrolled<WenoQPX_MinusFunctor> wenoqpx_unrolled;
#endif
	
	timer.start();
#pragma omp parallel 
	{
		for(int i=0; i<NTIMES; ++i)
			wenoqpx_unrolled(a, b, c, d, e, result, NENTRIES);
	}
	const double t_unrolled = timer.stop();
	
	
	
#ifdef _PROFILE_
	{
		static bool binit = false;
		if (!binit)
		{
			MPI_Init(&argc, &argv);
			binit = true;
		}
		printf("Profiling now ...");
		HPM_Start(const_cast<char *>(("WenoCPP"+ benchmark_name).c_str()));
		for(int i=0; i<NTIMES; ++i)
			weno_reference(a, b, c, d, e, result, NENTRIES);
		HPM_Stop(const_cast<char *>(("WenoCPP"+ benchmark_name).c_str()));
		
		HPM_Start(const_cast<char *>(("WenoQPX"+ benchmark_name).c_str()));
		for(int i=0; i<NTIMES; ++i)
			wenoqpx(a, b, c, d, e, result, NENTRIES);
		HPM_Stop(const_cast<char *>(("WenoQPX"+ benchmark_name).c_str()));
		
		HPM_Start(const_cast<char *>(("WenoQPXUnrolled"+ benchmark_name).c_str()));
		for(int i=0; i<NTIMES; ++i)
			wenoqpx_unrolled(a, b, c, d, e, result, NENTRIES);
		HPM_Stop(const_cast<char *>(("WenoQPXUnrolled"+ benchmark_name).c_str()));
		
		printf("done.\n");
	}
#endif
	
	if (verbose)
		for(int i=0; i<NENTRIES; ++i)
			printf("SOLVER[%d] = %f\n", i, result[i]);
	*/
	const double tol = precdiv? std::numeric_limits<Real>::epsilon() * 250 : 1e-5;
	printf("minus: verifying accuracy with tolerance %.5e...", tol);	
	check_error(tol, gold, result, NENTRIES);
	printf("passed!\n");
	
	printf("plus: verifying accuracy with tolerance %.5e...", tol);	
	check_error(tol, goldp, resultp, NENTRIES);
	printf("passed!\n");
	
	/*print_performance(NENTRIES, tref / NTIMES, "gold", true);
	print_performance(NENTRIES, t / NTIMES, "solver", false);
	print_performance(NENTRIES, t_unrolled / NTIMES, "unrolled solver", false);
	*/
	free(a);
	free(b);
	free(c);
	free(d);
	free(e);
	free(gold);
	free(result);
}

int main (int argc, char *  argv[]) 
{	
    std::cout << "Hello, weno benchmark!\n";
	
	const bool debug = false;
	
	if (debug)
	{
		benchmark(argc, argv, 4, 1, true, "debug");
		return 0;
	}
	
	/* performance on cache hits */
	{
		const double desired_kb =  16 * 4 * 0.5; /* we want to fill 50% of the dcache */
		const int nentries =  16 * (int)(pow(32 + 6, 2) * 4);//floor(desired_kb * 1024. / 7 / sizeof(Real));
		const int ntimes = (int)std::floor(2. / (1e-7 * nentries));
		
		for(int i=0; i<4; ++i)
		{
			printf("*************** PEAK-LIKE BENCHMARK (RUN %d) **************************\n", i);
			benchmark(argc, argv, nentries, ntimes, false, "cache");
		}
	}
	
	/* performance on data streams */
	{
		const double desired_mb =  128 * 4;
		const int nentries =  (int)floor(desired_mb * 1024. * 1024. / 7 / sizeof(Real));
		
		for(int i=0; i<4; ++i)
		{
			printf("*************** STREAM-LIKE BENCHMARK (RUN %d) **************************\n", i);
			benchmark(argc, argv, nentries, 1, false, "stream");
		}
	}
	
#ifdef _PROFILE_
	MPI_Finalize();
#endif
    return 0;
}

