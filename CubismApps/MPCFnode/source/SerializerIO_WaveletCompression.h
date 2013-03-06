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
#include <fstream>

#include "SerializerIO.h"
#include "WaveletCompressor.h"

template<typename GridType, typename Streamer>
class SerializerIO_WaveletCompression 
{
	enum { 
		NCHANNELS = Streamer::channels,
		NPTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_
	};
	
	const string filextension;
	
	Real threshold;
	bool normalized;
	
	bool halffloat;
	
	Real minval[NCHANNELS], maxval[NCHANNELS];
	
public:	
	
	SerializerIO_WaveletCompression():threshold(std::numeric_limits<Real>::epsilon()), normalized(false), halffloat(false), filextension(".zi4"){}
	
	void set_threshold(const Real threshold) { this->threshold = threshold; }
	
	void normalize(const Real minval[NCHANNELS], const Real maxval[NCHANNELS]) 
	{
		for(int i = 0; i < NCHANNELS; ++i)
		{
			this->minval[i] = minval[i];
			this->maxval[i] = maxval[i];

			assert(minval[i] <= maxval[i]);
		}
		
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
		
		ofstream myfile ( (fileName + filextension).c_str() );		
		assert(myfile.good());
		
		{
			int one = 1;
			bool isone = *(char *)(&one);
			
			myfile << "Endianess: " << (isone ? "little" : "big") << "\n";
		}
		
		myfile << "sizeofReal: " << sizeof(Real) << "\n";
		myfile << "Blocks: " << NBLOCKS << "\n";
		myfile << "Blocksize: " << _BLOCKSIZE_ << "\n";
		myfile << "Channels: " << NCHANNELS << "\n";
		myfile << "HalfFloat: " << (halffloat ? "yes" : "no") << "\n";
		myfile << "Normalized: " << (normalized ? "yes" : "no") << "\n";
		myfile << "==============START-BINARY-DATA===============\n";
		
		if (normalized)
		{
			myfile.write((const char *)minval, sizeof(minval));
			myfile.write((const char *)maxval, sizeof(maxval));
			
			written_bytes += sizeof(minval) + sizeof(maxval);
		}
		
		for(int i = 0; i < NBLOCKS; i++)
		{
			FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
			
			Real myaosbuffer[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][NCHANNELS];
			
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
						streamer.operate(b(ix, iy, iz), myaosbuffer[iz][iy][ix]);
			
			const int info[4] = { i, vInfo[i].index[0], vInfo[i].index[1], vInfo[i].index[2] };
			
			myfile.write((const char *)info, sizeof(info));
			
			for(int channel = 0; channel < NCHANNELS; ++channel)
			{
				Real mysoabuf[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
				
				for(int iz=0; iz<FluidBlock::sizeZ; iz++)
					for(int iy=0; iy<FluidBlock::sizeY; iy++)
						for(int ix=0; ix<FluidBlock::sizeX; ix++)
							mysoabuf[iz][iy][ix] = myaosbuffer[iz][iy][ix][channel];
				
				if (normalized && maxval[channel] - minval[channel] > numeric_limits<Real>::epsilon() * 100)
				{		
					const Real a = 1./(maxval[channel] - minval[channel]);
					const  Real b = - minval[channel]/(maxval[channel] - minval[channel]);
						
					Real * const x = &mysoabuf[0][0][0];
						
					for(int i = 0; i < NPTS; ++i)
							x[i] = a * x[i] + b;
				}
				
				WaveletCompressor compressor;
				
				const int nbytes = (int)compressor.compress(threshold, halffloat, mysoabuf);
				myfile.write((const char *)&nbytes, sizeof(nbytes));
				myfile.write((const char *)compressor.data(), sizeof(char) * nbytes);
				
				written_bytes += nbytes + sizeof(nbytes);
			}	
		}
		
		const double naive = (double)(sizeof(Real) * _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_ * NCHANNELS * NBLOCKS);
		const double compr = written_bytes;
		
		printf("overall compression-rate: %f\n", naive / compr);
	}
	
	template<typename T>
	double discrep(T ref[], T val[], const int N)
	{
		double gerr = 0;
		
		for(int i=0; i<N; ++i)
		{
			assert(!std::isnan(ref[i]));
			assert(!std::isnan(val[i]));
			
			const double err = ref[i] - val[i];
			const double relerr = err/std::max((Real)1e-6, std::max(fabs(val[i]), fabs(ref[i]))); 
			const double lerr = min(fabs(err), fabs(relerr));
			
			gerr = max(lerr, gerr);
		}
		
		return gerr;
	}
	
	void Read(GridType & inputGrid, string fileName, Streamer streamer = Streamer())
	{
		printf("hallo read!\n");
		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		const int NBLOCKS = vInfo.size();
		
		size_t read_bytes = 0;
		
		ifstream myfile( (fileName + filextension).c_str());
		
		{
			string word;
			string strendianess, strsizeofreal, strnblocks, strblocksize, strnchannels, strhalffloat, strnormalized;
			myfile >> word >> strendianess;
			myfile >> word >> strsizeofreal;
			myfile >> word >> strnblocks;
			myfile >> word >> strblocksize;
			myfile >> word >> strnchannels;
			myfile >> word >> strhalffloat;
			myfile >> word >> strnormalized;
			myfile >> word ;
			myfile.ignore(1, '\n');
			
			assert(strendianess == "little");
			assert(atoi(strsizeofreal.c_str()) == sizeof(Real));
			assert(atoi(strnblocks.c_str()) == NBLOCKS);
			assert(atoi(strblocksize.c_str()) == _BLOCKSIZE_);
			assert(atoi(strnchannels.c_str()) == NCHANNELS);
			assert(strhalffloat == (halffloat ? "yes" : "no"));
			assert(strnormalized == (normalized ? "yes" : "no"));
			
			cout  << "Endianess: <" << strendianess << ">\n"; 
			cout  << "SizeOfReal: <" << strsizeofreal << ">\n"; 
			cout  << "Blocks: <" << strnblocks<< ">\n"; 
			cout  << "Blocksize: <" << strblocksize<< ">\n";
			cout  << "Channels: <" << strnchannels<< ">\n"; 
			cout  << "HalfFloat: <" << strhalffloat<< ">\n"; 
			cout  << "Normalized: <" << strnormalized<< ">\n";
			
			cout.flush();
			
			if (normalized)
			{
				myfile.read((char *)minval, sizeof(minval));
				myfile.read((char *)maxval, sizeof(maxval));
				
				read_bytes += sizeof(minval) + sizeof(maxval);
			}
			
			for(int _i = 0; _i < NBLOCKS; ++_i)
			{
				int info[4]; //= { i, vInfo[i].index[0], vInfo[i].index[1], vInfo[i].index[2] };
				myfile.read((char *)info, sizeof(info));
				
				const int i = info[0];
				
				assert(vInfo[i].index[0] == info[1]);
				assert(vInfo[i].index[1] == info[2]);
				assert(vInfo[i].index[2] == info[3]);
				
				Real myaosbuffer[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][NCHANNELS];
				
				for(int channel = 0; channel < NCHANNELS; ++channel)
				{
					Real mysoabuf[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
					
					WaveletCompressor compressor;
					
					int nbytes;
					myfile.read((char *)&nbytes, sizeof(nbytes));
					myfile.read((char *)compressor.data(), sizeof(char) * nbytes);
					read_bytes += nbytes + sizeof(nbytes);
					
					compressor.decompress(halffloat, nbytes, mysoabuf);
					
					if (normalized && maxval[channel] - minval[channel] > numeric_limits<Real>::epsilon() * 100)
					{								
						const Real b = minval[channel];
						const Real a = maxval[channel] - b; 
						
						Real * const x = &mysoabuf[0][0][0];
						for(int i = 0; i < NPTS; ++i)
							x[i] = a * x[i] + b;
					}
					
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
						for(int iy=0; iy<FluidBlock::sizeY; iy++)
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
								myaosbuffer[iz][iy][ix][channel] = mysoabuf[iz][iy][ix];
				}
				
				FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
				
				/* we are not allowed to write into it, but at least we can check the accuracy */
				Real myref[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][NCHANNELS];
				
				for(int iz=0; iz<FluidBlock::sizeZ; iz++)
					for(int iy=0; iy<FluidBlock::sizeY; iy++)
						for(int ix=0; ix<FluidBlock::sizeX; ix++)
							streamer.operate(b(ix, iy, iz), myref[iz][iy][ix]);
				
				const double mydiscrep = discrep(&myref[0][0][0][0], &myaosbuffer[0][0][0][0], NCHANNELS * NPTS);
				printf("discrepancy: %e\n", mydiscrep);
				//assert(mydiscrep < 10 * threshold);
			}	
		}
		
		printf("Alright. SerializerIO_WaveletCompression::Read triggers a smooth exit here.\n");
		exit(0);
	}
};