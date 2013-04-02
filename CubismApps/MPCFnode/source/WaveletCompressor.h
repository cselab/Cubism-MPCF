/*
 *  WaveletCompressor.h
 *  
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <cassert>
#include <cstdlib>
#include <cstring>

#include <zlib.h>

#include "FullWaveletTransform.h"
static const bool lifting_scheme = false;

template<int DATASIZE1D, typename DataType>
class WaveletCompressorGeneric
{
protected:
	
	enum 
	{ 
		BS3 = DATASIZE1D * DATASIZE1D * DATASIZE1D,
		BITSETSIZE = (BS3 + 7) / 8,
		BUFMAXSIZE = BITSETSIZE + sizeof(DataType) * BS3
	};
	
	WaveletsOnInterval::FullTransform<DATASIZE1D, lifting_scheme> full;

private:
	__attribute__((__aligned__(16))) unsigned char bufcompression[BUFMAXSIZE];
	
	size_t bufsize;
	
public: 
			
	virtual void * data()
	{
		return bufcompression;
	}
	
	virtual	size_t compress(const float threshold, const bool float16, const DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D]);
	
	virtual void decompress(const bool float16, size_t bytes, DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D]);
};

template<int DATASIZE1D, typename DataType>
class WaveletCompressorGeneric_zlib: protected WaveletCompressorGeneric<DATASIZE1D, DataType>
{
	__attribute__((__aligned__(16))) unsigned char bufzlib[WaveletCompressorGeneric<DATASIZE1D, DataType>::BUFMAXSIZE];

public:

	void * data()
	{
		return bufzlib;
	}

	size_t compress(const float threshold, const bool float16, const DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
	{
		int compressedbytes = 0;
		
		const size_t ninputbytes = WaveletCompressorGeneric<DATASIZE1D, DataType>::compress(threshold, float16, data);

		z_stream datastream = { 0 };
		datastream.total_in = datastream.avail_in = ninputbytes;
		datastream.total_out = datastream.avail_out = WaveletCompressorGeneric<DATASIZE1D, DataType>::BUFMAXSIZE;
		datastream.next_in = (unsigned char*) WaveletCompressorGeneric<DATASIZE1D, DataType>::data();
		datastream.next_out = bufzlib;

		if (Z_OK == deflateInit(&datastream, Z_DEFAULT_COMPRESSION) && Z_STREAM_END == deflate(&datastream, Z_FINISH))
			compressedbytes = datastream.total_out;
		else
		{
			printf("ZLIB COMPRESSION FAILURE!!\n");
			abort();
		}

		deflateEnd( &datastream );


		return compressedbytes;
	}

	void decompress(const bool float16, size_t ninputbytes, DataType data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
	{
		int decompressedbytes = 0;

		z_stream datastream = { 0 };
		datastream.total_in = datastream.avail_in = ninputbytes;
		datastream.total_out = datastream.avail_out = WaveletCompressorGeneric<DATASIZE1D, DataType>::BUFMAXSIZE;
		datastream.next_in = bufzlib;
		datastream.next_out = (unsigned char*) WaveletCompressorGeneric<DATASIZE1D, DataType>::data();

		const int retval = inflateInit(&datastream);

		if (retval == Z_OK && inflate(&datastream, Z_FINISH))
			decompressedbytes = datastream.total_out;
                else
                {
                        printf("ZLIB DECOMPRESSION FAILURE!!\n");
                        abort();
                }

		inflateEnd(&datastream);
		
		WaveletCompressorGeneric<DATASIZE1D, DataType>::decompress(float16, decompressedbytes, data);
	}
};

#ifdef _BLOCKSIZE_
typedef WaveletCompressorGeneric<_BLOCKSIZE_, Real> WaveletCompressor;
typedef WaveletCompressorGeneric_zlib<_BLOCKSIZE_, Real> WaveletCompressor_zlib;
#endif
