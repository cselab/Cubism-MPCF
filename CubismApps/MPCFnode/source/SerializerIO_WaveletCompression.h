/*
 *  SerializerIO_WaveletCompression.h
 *  
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <limits>

#include "SerializerIO.h"
#include "WaveletCompressor.h"

template<typename GridType, typename Streamer>
class SerializerIO_WaveletCompression 
{
	enum { 
		NCHANNELS = Streamer::channels,
		NPTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_
	};
	
	Real threshold;
	bool normalized;
	
	bool halffloat;
	
	Real minval[NCHANNELS], maxval[NCHANNELS];
	
public:	
	
	SerializerIO_WaveletCompression():threshold(std::numeric_limits<Real>::epsilon()), normalized(false),	halffloat(false){}
	
	void set_threshold(const Real threshold) { this->threshold = threshold; }
	
	void normalize(const Real minval[NCHANNELS], const Real maxval[NCHANNELS]) 
	{
		for(int i = 0; i < NCHANNELS; ++i)
			this->minval[i] = minval[i];
		
		for(int i = 0; i < NCHANNELS; ++i)
			this->maxval[i] = maxval[i];
		
		normalized = true;
	} 
	
	void float16()
	{
		halffloat = true;
	}
	
	void Write(GridType & inputGrid, string fileName, 
			   Streamer streamer = Streamer())
	{
		size_t written_bytes = 0;
		
		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		const int NBLOCKS = vInfo.size();
		
		FILE * file = fopen( (fileName + ".iw4").c_str(), "wb");
		
		assert(file != NULL);
		
		{
			int one = 1;
			bool isone = *(char *)(&one);
			
			fprintf(file, "Endianess: %s\n", isone ? "little" : "big");
		}
		
		fprintf(file, "Size of real: %d\n", sizeof(Real));
		fprintf(file, "Blocks: %d\n", (int)NBLOCKS);
		fprintf(file, "Block size: %d\n", _BLOCKSIZE_);
		
		if (normalized)
		{
			fprintf(file, "data is normalized\n");
			
			fwrite(minval, sizeof(minval), 1, file);
			fwrite(maxval, sizeof(maxval), 1, file);
		}
		
		for(int i = 0; i < NBLOCKS; i++)
		{
			FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
			
			Real myaosbuffer[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][NCHANNELS];
			
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
						streamer.operate(b(ix, iy, iz), myaosbuffer[iz][iy][ix]);
			
			const int info[] = { i, vInfo[i].index[0], vInfo[i].index[1], vInfo[i].index[2] };
			
			fwrite(info, sizeof(int), 4, file);
			
			for(int channel = 0; channel < NCHANNELS; ++channel)
			{
				Real mysoabuf[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
				
				for(int iz=0; iz<FluidBlock::sizeZ; iz++)
					for(int iy=0; iy<FluidBlock::sizeY; iy++)
						for(int ix=0; ix<FluidBlock::sizeX; ix++)
							mysoabuf[iz][iy][ix] = myaosbuffer[iz][iy][ix][channel];
				
				if (normalized)
				{		
					const Real a = 1./max((maxval[channel] - minval[channel]), numeric_limits<Real>::epsilon());
					const Real b = - minval[channel] * a;
					
					Real * const x = &mysoabuf[0][0][0];
					
					for(int i = 0; i < NPTS; ++i)
						x[i] = a * x[i] + b;
				}
				
				WaveletCompressor compressor;
				
				const size_t nbytes = compressor.compress(threshold, halffloat, mysoabuf);
				fwrite(&nbytes, sizeof(nbytes), 1, file);
				
				fwrite(compressor.data(), sizeof(char), nbytes, file);
				
				written_bytes += nbytes;
			}	
		}
		
		fclose(file);
		
		const double naive = (double)(sizeof(Real) * _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_ * NCHANNELS * NBLOCKS);
		const double compr = written_bytes;
		
		printf("overall compression-rate: %f\n", naive / compr);
	}
};