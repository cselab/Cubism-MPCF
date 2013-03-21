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
#include <omp.h>

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
	
	bool halffloat, isthreaded;
	
	Real minval[NCHANNELS], maxval[NCHANNELS];
	
public:	
	
	SerializerIO_WaveletCompression():threshold(std::numeric_limits<Real>::epsilon()), halffloat(false), isthreaded(true), filextension(".zi4"){}
	
	void set_threshold(const Real threshold) { this->threshold = threshold; }
	
	void singlethreaded() { isthreaded = false; } 
		
	void float16()
	{
		halffloat = true;
	}
	
	size_t _write_naive(ofstream& myfile, vector<BlockInfo>& vInfo, const int NBLOCKS, Streamer streamer)
	{
		size_t written_bytes = 0;

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
				
				WaveletCompressor compressor;
				
				const int nbytes = (int)compressor.compress(threshold, halffloat, mysoabuf);
				myfile.write((const char *)&nbytes, sizeof(nbytes));
				myfile.write((const char *)compressor.data(), sizeof(char) * nbytes);
				written_bytes += nbytes + sizeof(nbytes);
			}	
		}
	}
	
	size_t _write_threaded(ofstream& myfile, vector<BlockInfo>& vInfo, const int NBLOCKS, Streamer streamer)
	{
		size_t written_bytes = 0;
		
#pragma omp parallel 
		{
			enum 
			{ 
				DESIREDMEM = 2 << 20, //3 megs is also cool?
				ENTRYSIZE = NCHANNELS * sizeof(WaveletCompressor) + sizeof(int) * 5,
				ENTRIES_CANDIDATE = DESIREDMEM / ENTRYSIZE,
				ENTRIES = ENTRIES_CANDIDATE ? ENTRIES_CANDIDATE : 1,
				BUFFERSIZE = ENTRIES * ENTRYSIZE,
				ALERT = (ENTRIES - 1) * ENTRYSIZE
			};
										
			size_t mybytes = 0;
			char buffer[BUFFERSIZE];

#pragma omp for
			for(int i = 0; i < NBLOCKS; i++)
			{
				FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
				
				Real myaosbuffer[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][NCHANNELS];
				
				for(int iz=0; iz<FluidBlock::sizeZ; iz++)
					for(int iy=0; iy<FluidBlock::sizeY; iy++)
						for(int ix=0; ix<FluidBlock::sizeX; ix++)
							streamer.operate(b(ix, iy, iz), myaosbuffer[iz][iy][ix]);
				
				const int info[4] = { i, vInfo[i].index[0], vInfo[i].index[1], vInfo[i].index[2] };
				
				memcpy(buffer + mybytes, info, sizeof(info));
				mybytes += sizeof(info);
				
				for(int channel = 0; channel < NCHANNELS; ++channel)
				{
					Real mysoabuf[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
					
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
						for(int iy=0; iy<FluidBlock::sizeY; iy++)
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
								mysoabuf[iz][iy][ix] = myaosbuffer[iz][iy][ix][channel];
					
					WaveletCompressor compressor;

					const int nbytes = (int)compressor.compress(threshold, halffloat, mysoabuf);
					memcpy(buffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);
					
					memcpy(buffer + mybytes, compressor.data(), sizeof(char) * nbytes);
					mybytes += nbytes;
				}

				if (mybytes >= ALERT) //flush: aquire mutex and write to file
				#pragma omp critical
				{
					myfile.write((const char *)&buffer, mybytes);
					written_bytes += mybytes;
					mybytes = 0;
				}
			}

			if (mybytes > 0) //remaining data in the buffer must be flushed: aquire mutex and write to file
			#pragma omp critical
			{      
				myfile.write((const char *)&buffer, mybytes);
				written_bytes += mybytes;
				mybytes = 0;
			} 
		}

		return written_bytes;
	}
	
	void Write(GridType & inputGrid, string fileName, 
			   Streamer streamer = Streamer())
	{		
		Timer timer;
		timer.start();
		
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
		myfile << "==============START-BINARY-DATA===============\n";
		
		size_t written_bytes = 0;

		if (!isthreaded)
			written_bytes = _write_naive(myfile, vInfo, NBLOCKS, streamer);	
		else
			written_bytes = _write_threaded(myfile, vInfo, NBLOCKS, streamer);	
		
		const double naive = (double)(sizeof(Real) * NPTS * NCHANNELS * NBLOCKS);
		const double compr = written_bytes;
		
		printf("overall compression-rate: %f, time spent: %.2f s\n", naive / compr, timer.stop());
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
			const double maxval = std::max(fabs(val[i]), fabs(ref[i]));
			const double relerr = err/std::max((Real)std::numeric_limits<Real>::epsilon()*10, std::max(fabs(val[i]), fabs(ref[i]))); 
			
			const double lerr = min( fabs(relerr), fabs(err));
			
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
			myfile >> word ;
			myfile.ignore(1, '\n');
			
			assert(strendianess == "little");
			assert(atoi(strsizeofreal.c_str()) == sizeof(Real));
			assert(atoi(strnblocks.c_str()) == NBLOCKS);
			assert(atoi(strblocksize.c_str()) == _BLOCKSIZE_);
			assert(atoi(strnchannels.c_str()) == NCHANNELS);
			assert(strhalffloat == (halffloat ? "yes" : "no"));
			
			cout  << "Endianess: <" << strendianess << ">\n"; 
			cout  << "SizeOfReal: <" << strsizeofreal << ">\n"; 
			cout  << "Blocks: <" << strnblocks<< ">\n"; 
			cout  << "Blocksize: <" << strblocksize<< ">\n";
			cout  << "Channels: <" << strnchannels<< ">\n"; 
			cout  << "HalfFloat: <" << strhalffloat<< ">\n"; 
			
			cout.flush();
			
			for(int _i = 0; _i < NBLOCKS; ++_i)
			{
				int info[4];
				myfile.read((char *)info, sizeof(info));
				
				const int i = info[0];
				//printf("%d %d %d %d\n", info[0], info[1], info[2], info[3]);		
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
					
					/* i used this code for debugging
						if (i == 27 && channel == 6)
						{
							ofstream fileref("reference.txt");
							ofstream fileres("results.txt");
							
							const int mychannel = 6;
							
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
							{	
								bool myinit = false;
								const int iz = 14;
								for(int iy=0; iy<FluidBlock::sizeY; iy+=2)
									//for(int iz=0; iz<FluidBlock::sizeZ; iz++)
								{
									Real myref[NCHANNELS];
									
									streamer.operate(b(ix, iy, iz), myref);
									
									if (!myinit)
									{
										myinit = true;
										fileres << ix << " ";
										fileref << ix << " ";
									}
									
									fileres <<  setprecision(10) << mysoabuf[iz][iy][ix] << " ";
									fileref <<  setprecision(10) << myref[mychannel] << " ";
									
								}
								fileres << "\n";
								fileref << "\n";
							}
						}*/
					
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
						for(int iy=0; iy<FluidBlock::sizeY; iy++)
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
								myaosbuffer[iz][iy][ix][channel] = mysoabuf[iz][iy][ix];
				}
				
				{
					/* we are not allowed to write into it, but at least we can check the accuracy */
					Real myref[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][NCHANNELS];
					
					FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;

					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
						for(int iy=0; iy<FluidBlock::sizeY; iy++)
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
								streamer.operate(b(ix, iy, iz), myref[iz][iy][ix]);

					const double mytol = max (threshold, numeric_limits<Real>::epsilon() * 5000);
					const double mydiscrep = discrep(&myref[0][0][0][0], &myaosbuffer[0][0][0][0], NCHANNELS * NPTS);

					if (mydiscrep > mytol) printf("block %d, discrepancy: %e, tol %e\n", i, mydiscrep, mytol);
					//assert(mydiscrep <= mytol);
				}
			}	
		}
		
		printf("Alright. SerializerIO_WaveletCompression::Read triggers a smooth exit here.\n");
		exit(0);
	}
};
