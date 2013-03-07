/*
 *  WaveletCompressor.h
 *  
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <cstring>

class WaveletCompressor
{
	enum 
	{ 
		BS3 = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_,
		BITSETSIZE = (BS3 + 7) / 8,
		BUFMAXSIZE = BITSETSIZE + sizeof(Real) * BS3
	};

	__attribute__((__aligned__(16))) unsigned char bufcompression[BUFMAXSIZE];
	
	size_t bufsize;
	
public: 
			
	void * data()
	{
		return bufcompression;
	}
	
	size_t compress(const Real threshold, const bool float16, const Real data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_]);
	
	void decompress(const bool float16, size_t bytes, Real data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_]);
};