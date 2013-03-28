/*
 *  WaveletCompressor.cpp
 *  
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#ifdef _SP_COMP_
typedef float Real;
#else
typedef double Real;
#endif

#include <cstdio>
#include <bitset>
#include <cassert>

using namespace std;

#include "WaveletCompressor.h"
#include "WaveletsOnInterval.h"

static const bool lifting_scheme = false;

template<int N>
void serialize_bitset(bitset<N> mybits, unsigned char * const buf, const int nbytes)
{
	assert(nbytes == (N + 7) / 8);
	
	const int nicebits = 8 * (N / 8);
		
	for(int i = 0, B = 0; i < nicebits; i += 8, ++B)
	{
		unsigned char c = 0;
		
		for(int b = 0; b < 8; ++b)
			c |= mybits[i + b] << b;
		
		buf[B] = c;
	}
	
	if (nicebits < N)
	{
		unsigned char c = 0;
		
		for(int b = nicebits; b < N; ++b)
			c |= mybits[b] << (b - nicebits);
		
		buf[nbytes - 1] = c;
	}
}

template<int N>
int deserialize_bitset(bitset<N>& mybits, const unsigned char * const buf, const int nbytes)
{
	assert(nbytes == (N + 7) / 8);
	
	int sum = 0;

	const int nicebits = 8 * (N / 8);
	
	for(int i = 0, B = 0; i < nicebits; i += 8, ++B)
	{
		const unsigned char c = buf[B];
		
		for(int b = 0; b < 8; ++b)
			sum += mybits[i + b] = (c >> b) & 1;		
	}
	
	if (nicebits < N)
	{
		const unsigned char c = buf[nbytes - 1];
		
		for(int b = nicebits; b < N; ++b)
			sum += mybits[b] = (c >> (b - nicebits)) & 1;
	}
	
	return sum;
}

unsigned short _cvt2f16(const float xfloat)
{
	assert(sizeof(unsigned short) == 2);

	const unsigned int x = *(unsigned int *)& xfloat;
	
	//copy sign
	unsigned short retval = (x & 0x80000000) >> 16; 
	
	//reshape exponent
	const int e = ((x & 0x7fffffff) >> 23) - 127;
	retval |= max(0, min(31, e + 15)) << 10;
	
	//reshape mantissa
	const unsigned int m = x & 0x007fffff;
	retval |= m >> 13;
	
	return retval;
}

float _cvtfromf16(const unsigned short f16)
{
	//copy sign
	int retval = (f16 & 0x8000) << 16;
	
	//reshape exponent
	const int e = (((0x1f << 10) & f16) >> 10) - 15;
	retval |= (e + 127) << 23;
	
	//reshape mantissa
	const int m = f16 & 0x3ff;
	retval |= m << 13;
	
	return *(float *)& retval;
}

size_t WaveletCompressor::compress(const Real threshold, const bool float16, const Real data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_])
{	
	WaveletsOnInterval::FullTransform<_BLOCKSIZE_, lifting_scheme> full;
		
	const Real * const src = &data[0][0][0];
	WaveletsOnInterval::FwtAp * const dst = &full.data[0][0][0];
	for(int i = 0; i < BS3; ++i) dst[i] = src[i];
	
	full.fwt();
	pair<vector<Real> , bitset<BS3> > survivors = full.threshold(threshold);

	size_t bytes_copied = 0;
	
	serialize_bitset<BS3>(survivors.second, bufcompression, BITSETSIZE);
	bytes_copied += BITSETSIZE;
	
	const int nelements = survivors.first.size();
	if (!float16)
	{
		memcpy(bufcompression + bytes_copied, &survivors.first.front(), sizeof(Real) * nelements);
		bytes_copied += sizeof(Real) * nelements;
	}
	else
		for(int i = 0; i < nelements; ++i, bytes_copied += 2)
		{
			unsigned short x = _cvt2f16(survivors.first[i]);
			
			bufcompression[bytes_copied + 0] = ((unsigned char *)&x)[0];
			bufcompression[bytes_copied + 1] = ((unsigned char *)&x)[1];			
		}
	
	return bytes_copied;
}

void WaveletCompressor::decompress(const bool float16, size_t bytes, Real data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_])
{
	assert((bytes - sizeof(bitset<BS3>)) % sizeof(Real) == 0);
	
	bitset<BS3> mask;
	const int expected = deserialize_bitset<BS3>(mask, bufcompression, BITSETSIZE);

	size_t bytes_read = BITSETSIZE;
	
	const int nelements = (bytes - bytes_read) / (float16 ? sizeof(unsigned short) : sizeof(Real));
	assert(expected == nelements);

	vector<Real> datastream(nelements);
	
	if (!float16)
		memcpy((void *)&datastream.front(), bufcompression + bytes_read, sizeof(Real) * nelements);	
	else
	{
		unsigned short elements[nelements];
		memcpy((void *)&elements, bufcompression + bytes_read, sizeof(unsigned short) * nelements);	
		
		for(int i = 0; i < nelements; ++i)
			datastream[i] =  _cvtfromf16(elements[i]);
	}
	
	WaveletsOnInterval::FullTransform<_BLOCKSIZE_, lifting_scheme> full;
	full.load(datastream, mask);
	full.iwt();

	Real * const dst = &data[0][0][0];
	const WaveletsOnInterval::FwtAp * const src = &full.data[0][0][0];
	for(int i = 0; i < BS3; ++i)
	{ 
		dst[i] = src[i];
		assert(!std::isnan(dst[i]));
	}
}
