/*
 *  SerializerIO_WaveletCompression_MPI.h
 *  
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#error This class (SerializerIO_WaveletCompression_MPI) is not ready. Please don't include it.

#pragma once


#include "SerializerIO_WaveletCompression.h"

/* this class will be provided by panos, into another header file */
#include "MPI_Streamer.h"

template<int DESIREDMEM, int NSLOTS>
class ChainedWriteBuffer_MPI: public ChainedWriteBuffer<DESIREDMEM, NSLOTS> 
{
protected:
	
	MPI_Streamer& mympistreamer;
	
	void _nonblocking_write(const int targetslot)
	{
		//panos, dont worry this method is not called concurrently
		mympistreamer.Write(this->buffers[targetslot], this->nbytes[targetslot] * sizeof(char));
	}
	
	void _wait_resolve(const int targetslot)
	{
		//panos, dont worry this method is not called concurrently
		mympistreamer.WaitResolve(this->buffers[targetslot]);
	}
	
public:
	
	ChainedWriteBuffer_MPI(MPI_Streamer& mympistreamer):
	ChainedWriteBuffer<DESIREDMEM, NSLOTS>(), mympistreamer(mympistreamer) { }
	
	void Reset() { this->_reset(); }
};

template<typename GridType, typename Streamer>
class SerializerIO_WaveletCompression_MPI: public SerializerIO_WaveletCompression<GridType, Streamer> 
{
	enum {
		NCHANNELS = Streamer::channels,
		NPTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_
	};

	int callscounter; //number of times the "Write" method has been called

	MPI_Streamer mympistreamer;

	//say we want a heap buffer of 50 MB, two slots
	typedef ChainedWriteBuffer_MPI<50, 2> MyChainedBuffer;
	MyChainedBuffer * chainedbuffer;

	template<typename CompressorType> size_t _to_mpi_file(vector<BlockInfo>& vInfo, const int NBLOCKS, Streamer streamer)
	{
		const int MYNTHREADS = omp_get_max_threads( );
		float t[MYNTHREADS];

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

			Timer timer;
			timer.start();

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

					CompressorType compressor;

					const int nbytes = (int)compressor.compress(this->threshold, this->halffloat, mysoabuf);
					memcpy(buffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);

					memcpy(buffer + mybytes, compressor.data(), sizeof(char) * nbytes);
					mybytes += nbytes;
				}

				if (mybytes >= ALERT)
				{
					void * dest = chainedbuffer->initiate_write(mybytes);
					memcpy(dest, buffer, mybytes);
					chainedbuffer->finalize_write(dest);
					mybytes = 0;
				}
			}

			if (mybytes > 0)
			{
				void * dest = chainedbuffer->initiate_write(mybytes);
				memcpy(dest, buffer, mybytes);
				chainedbuffer->finalize_write(dest);
				mybytes = 0;
			}

			t[omp_get_thread_num()] = timer.stop();
		}

		//print out the thread timings to check the imbalance
		if(false)
			for(int i=0; i<MYNTHREADS; ++i)
				printf("t[%d] = %.2e s\n", i, t[i]);

		return chainedbuffer->get_written_bytes();
	}

	public:

	SerializerIO_WaveletCompression_MPI(): callscounter(0), mympistreamer(), chainedbuffer(NULL)
	{
		//printf("This class is not ready");
		chainedbuffer = new MyChainedBuffer(mympistreamer);
	}

	~SerializerIO_WaveletCompression_MPI()
	{	
		if (callscounter > 0)
			mympistreamer.Flush();

		delete chainedbuffer;
	}

	void Write(GridType & inputGrid, string fileName, 
			Streamer streamer = Streamer())
	{		
		Timer timer;
		timer.start();

		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		const int NBLOCKS = vInfo.size();

		if (callscounter > 0)
		{
			mympistreamer.Flush();
			chainedbuffer->Reset();
		}

#if 1
		{
		ostringstream ss;
		
	        int one = 1;
	        bool isone = *(char *)(&one);

		ss << "Endianess: "  << (isone ? "little" : "big") << endl;
		ss << "sizeofReal: " << sizeof(Real) << endl;
		ss << "Blocks: " << NBLOCKS << endl;
		ss << "Blocksize: " << _BLOCKSIZE_ << endl;
		ss << "Channels: " <<  NCHANNELS << endl;
//		ss << "HalfFloat: " << (halffloat ? "yes" : "no")  << endl;
//		ss << "Normalized: " << (normalized ? "yes" : "no") << endl;
		ss << "Groupsize: " << mympistreamer.groupsize << endl;
		ss << "First_Rank: " << mympistreamer.first << endl;
		ss << "Last_Rank: " << mympistreamer.last << endl;
		ss << "==============START-BINARY-DATA==============" << endl;

//		printf("header\n"); 
//		cout << ss.str();
//		cout << ss.str().size();

		mympistreamer.SetHeader((char *)ss.str().c_str());
		}
#endif
		mympistreamer.Init(fileName);

		//send to panos gigantic compressed streams
		//what do we do with the written bytes, panos? we sum them all?
		const size_t written_bytes = _to_mpi_file<WaveletCompressor_zlib>(vInfo, NBLOCKS, streamer);	

		chainedbuffer->flush();

		const double naive = (double)(sizeof(Real) * NPTS * NCHANNELS * NBLOCKS);
		const double compr = written_bytes;

		printf("overall compression-rate: %f, time spent: %.2f s\n", naive / compr, timer.stop());

		callscounter++;
	}
};

