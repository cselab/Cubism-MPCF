/*
 *  check_error.h
 *  ComputeExpansions
 *
 *  Created by Diego Rossinelli on 1/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cmath>
#include <limits>

using namespace std;

template<typename T>
void check_error(const double tol, T ref[], T val[], const int N)
{
	static const bool verbose = false;
	int nan_ctr = 0;
	int dif_ctr = 0;

//	printf("%d\n",  N);
	for(int i=0; i<N; ++i)
	{
		assert(!std::isnan(ref[i]));
		assert(!std::isnan(val[i]));

		if (isnan(ref[i]) || isnan(val[i])) {
			nan_ctr++;
			continue;
		}
		const double err = ref[i] - val[i];
//		const double relerr = err/std::max((double)std::numeric_limits<T>::epsilon(), std::max(fabs(val[i]), fabs(ref[i])));
		const double relerr = err/std::max((double)std::numeric_limits<T>::epsilon(), (double)std::max(fabs(val[i]), fabs(ref[i])));

		if (verbose) printf("+%1.1e,", relerr);

		if (fabs(relerr) >= tol && fabs(err) >= tol) {
			dif_ctr++;
			printf("\n%d: %e %e -> %e %e\n", i, ref[i], val[i], err, relerr);
		}

		assert(fabs(relerr) < tol || fabs(err) < tol);
	}

	if (verbose) printf("\t");

	if (true) { //(verbose)
		if (nan_ctr) printf("\nnan_ctr=%d", nan_ctr);
		if (dif_ctr) printf("\ndif_ctr=%d", dif_ctr);
		if (nan_ctr || dif_ctr) printf("\n");
	}
}


template<typename T>
void check_error(const double tol, T ref, T val)
{
	check_error(tol, &ref, &val, 1);
}

