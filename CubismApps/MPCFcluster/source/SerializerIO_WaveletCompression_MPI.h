/*
 *  SerializerIO_WaveletCompression_MPI.h
 *  
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once


#include "SerializerIO_WaveletCompression.h"

/* this class will be provided by panos, into another header file */
class MPI_Streamer //DUMMY IMPL - TRUE IMPL PROVIDED BY PANOS
{
	public:
		MPI_Streamer(int _groupsize=-1) { }

		~MPI_Streamer() { }

		void Init(string _filename) { }

		void Write(char *buf, int nbytes) { }
		void WaitResolve(char * buf) { } //panos, this is the new method we agreed upon this morning
		void Flush() { }
};

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
			ChainedWriteBuffer<DESIREDMEM, NSLOTS>(), mympistreamer(mympistreamer)
	{
	}
};


template<typename GridType, typename Streamer>
class SerializerIO_WaveletCompression_MPI: public SerializerIO_WaveletCompression<GridType, Streamer> 
{
	enum {
		NCHANNELS = Streamer::channels,
		NPTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_
	};

	MPI_Streamer mympistreamer;

	template<typename CompressorType> size_t _to_mpi_file(vector<BlockInfo>& vInfo, const int NBLOCKS, Streamer streamer)
	{
		//say we want a heap buffer of 50 MB, two slots
		ChainedWriteBuffer_MPI<50, 2> * chainedbuffer = new ChainedWriteBuffer_MPI<50, 2>(mympistreamer);

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
		//if(false)
		for(int i=0; i<MYNTHREADS; ++i)
			printf("t[%d] = %.2e s\n", i, t[i]);

		size_t written_bytes = chainedbuffer->get_written_bytes();

		delete chainedbuffer;
		
		return written_bytes;
	}

	public:
	void Write(GridType & inputGrid, string fileName, 
			Streamer streamer = Streamer())
	{		
		Timer timer;
		timer.start();

		vector<BlockInfo> vInfo = inputGrid.getBlocksInfo();
		const int NBLOCKS = vInfo.size();

		mympistreamer.Init(fileName);
		
		//send to panos gigantic compressed streams
		//what do we do with the written bytes, panos? we sum them all?
		const size_t written_bytes = _to_mpi_file<WaveletCompressor_zlib>(vInfo, NBLOCKS, streamer);	
		
		mympistreamer.Flush();

		const double naive = (double)(sizeof(Real) * NPTS * NCHANNELS * NBLOCKS);
		const double compr = written_bytes;

		printf("overall compression-rate: %f, time spent: %.2f s\n", naive / compr, timer.stop());
	}
};

