/*
 *  WaveletTexture3D.cpp
 *  
 *
 *  Created by Diego Rossinelli on 3/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include <cstdio>
#include <cstring>
#include <vector>

using namespace std;

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include "WaveletTexture3D.h"
#include "../../MPCFnode/source/WaveletCompressor.h"

//typedef WaveletCompressorGeneric<_VOXELS_> WaveletCompressorASD;

vector<unsigned char> WaveletTexture3D::compress() 
{ 
	printf("compression here...\n"); return vector<unsigned char>();
	WaveletCompressorGeneric<_VOXELS_>ss asd;
	
	asd.compress(1e-3, true, data);
}

