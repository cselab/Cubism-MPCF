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

template<int DATASIZE1D, typename DataType>
size_t WaveletCompressorGeneric<DATASIZE1D, DataType>::compress(const float threshold, const bool float16, const DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
{			
	const DataType * const src = &data[0][0][0];
	WaveletsOnInterval::FwtAp * const dst = &full.data[0][0][0];
	for(int i = 0; i < BS3; ++i) dst[i] = src[i];
	
	full.fwt();
	pair<vector<DataType> , bitset<BS3> > survivors = full.template threshold<DataType>(threshold);

	size_t bytes_copied = 0;
	
	serialize_bitset<BS3>(survivors.second, bufcompression, BITSETSIZE);
	bytes_copied += BITSETSIZE;
	
	const int nelements = survivors.first.size();
	if (!float16)
	{
		memcpy(bufcompression + bytes_copied, &survivors.first.front(), sizeof(DataType) * nelements);
		bytes_copied += sizeof(DataType) * nelements;
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

template<int DATASIZE1D, typename DataType>
void WaveletCompressorGeneric<DATASIZE1D, DataType>::decompress(const bool float16, size_t bytes, DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
{
	assert((bytes - sizeof(bitset<BS3>)) % sizeof(DataType) == 0);
	
	bitset<BS3> mask;
	const int expected = deserialize_bitset<BS3>(mask, bufcompression, BITSETSIZE);

	size_t bytes_read = BITSETSIZE;
	
	const int nelements = (bytes - bytes_read) / (float16 ? sizeof(unsigned short) : sizeof(DataType));
	assert(expected == nelements);

	vector<DataType> datastream(nelements);
	
	if (!float16)
		memcpy((void *)&datastream.front(), bufcompression + bytes_read, sizeof(DataType) * nelements);	
	else
	{
		unsigned short elements[nelements];
		memcpy((void *)&elements, bufcompression + bytes_read, sizeof(unsigned short) * nelements);	
		
		for(int i = 0; i < nelements; ++i)
			datastream[i] =  _cvtfromf16(elements[i]);
	}
	
	full.load(datastream, mask);
	full.iwt();

	DataType * const dst = &data[0][0][0];
	const WaveletsOnInterval::FwtAp * const src = &full.data[0][0][0];
	for(int i = 0; i < BS3; ++i)
	{ 
		dst[i] = src[i];
		assert(!std::isnan(dst[i]));
	}
}

#ifdef _BLOCKSIZE_
template class WaveletCompressorGeneric<_BLOCKSIZE_, Real>;
template class WaveletCompressorGeneric_zlib<_BLOCKSIZE_, Real>;
#endif

#ifdef _VOXELS_ //mammamia whattahack
template class WaveletCompressorGeneric<_VOXELS_, float>;
template class WaveletCompressorGeneric_zlib<_VOXELS_, float>;
#endif
