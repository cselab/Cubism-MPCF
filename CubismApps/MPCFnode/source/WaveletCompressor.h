/*
 *  WaveletCompressor.h
 *  
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include <cstdlib>
#include <cstring>

#include <zlib.h>

template<int DATASIZE1D>
class WaveletCompressorGeneric
{
protected:
	enum 
	{ 
		BS3 = DATASIZE1D * DATASIZE1D * DATASIZE1D,
		BITSETSIZE = (BS3 + 7) / 8,
		BUFMAXSIZE = BITSETSIZE + sizeof(Real) * BS3
	};

private:
	__attribute__((__aligned__(16))) unsigned char bufcompression[BUFMAXSIZE];
	
	size_t bufsize;
	
public: 
			
	virtual void * data()
	{
		return bufcompression;
	}
	
	virtual	size_t compress(const Real threshold, const bool float16, const Real data[DATASIZE1D][DATASIZE1D][DATASIZE1D]);
	
	virtual void decompress(const bool float16, size_t bytes, Real data[DATASIZE1D][DATASIZE1D][DATASIZE1D]);
};

template<int DATASIZE1D>
class WaveletCompressorGeneric_zlib: protected WaveletCompressorGeneric<DATASIZE1D>
{
	__attribute__((__aligned__(16))) unsigned char bufzlib[WaveletCompressorGeneric<DATASIZE1D>::BUFMAXSIZE];

public:

	void * data()
	{
		return bufzlib;
	}

	size_t compress(const Real threshold, const bool float16, const Real data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
	{
		int compressedbytes = 0;
		
		const size_t ninputbytes = WaveletCompressorGeneric<DATASIZE1D>::compress(threshold, float16, data);

		z_stream datastream = {0};
		datastream.total_in = datastream.avail_in = ninputbytes;
		datastream.total_out = datastream.avail_out = WaveletCompressorGeneric<DATASIZE1D>::BUFMAXSIZE;
		datastream.next_in = (unsigned char*) WaveletCompressorGeneric<DATASIZE1D>::data();
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

	void decompress(const bool float16, size_t ninputbytes, Real data[DATASIZE1D][DATASIZE1D][DATASIZE1D])
	{
		int decompressedbytes = 0;

		z_stream datastream = {0};
		datastream.total_in = datastream.avail_in = ninputbytes;
		datastream.total_out = datastream.avail_out = WaveletCompressorGeneric<DATASIZE1D>::BUFMAXSIZE;
		datastream.next_in = bufzlib;
		datastream.next_out = (unsigned char*) WaveletCompressorGeneric<DATASIZE1D>::data();

		const int retval = inflateInit(&datastream);

		if (retval == Z_OK && inflate(&datastream, Z_FINISH))
			decompressedbytes = datastream.total_out;
                else
                {
                        printf("ZLIB DECOMPRESSION FAILURE!!\n");
                        abort();
                }

		inflateEnd(&datastream);
		
		WaveletCompressorGeneric<DATASIZE1D>::decompress(float16, decompressedbytes, data);
	}
};

typedef WaveletCompressorGeneric<_BLOCKSIZE_> WaveletCompressor;
typedef WaveletCompressorGeneric_zlib<_BLOCKSIZE_> WaveletCompressor_zlib;
