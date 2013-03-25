/*
 *  SerializerIO_WaveletCompression_MPI_Simple.h
 *  ... te lo do io il simple, porco porco
 *
 *  Created by Diego Rossinelli on 3/4/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <typeinfo>
#include <sstream>
#include <omp.h>

using namespace std;

#include "WaveletCompressor.h"

inline int deflate_inplace(z_stream *strm, unsigned char *buf, unsigned len, unsigned *max);

template<typename GridType, typename IterativeStreamer>
class SerializerIO_WaveletCompression_MPI_Simple
{
	enum 
	{
		NCHANNELS = IterativeStreamer::channels,
		NPTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_
	};
	
	struct BlockMetadata { int idcompression, subid, ix, iy, iz; };
		
	vector< BlockMetadata > myblockindices; //tells in which compressed chunk is any block, nblocks
	vector< size_t > lut_compression; //tells the number of compressed chunk, and where do they start, nchunks + 2
	vector< unsigned char > allmydata; //buffer with the compressed data
	size_t written_bytes;
	
	Real threshold;
	bool halffloat, verbosity;
	int callscounter; //number of times the "Write" method has been called

	vector< float > workload; //per-thread cpu time for imbalance insight

	void _zcompress_and_flush(unsigned char inputbuffer[], int& bufsize, const int maxsize, BlockMetadata metablocks[], int& nblocks)
	{
		//0. setup
		//1. compress the data with zlib, obtain zptr, zbytes
		//2. obtain an offset from allmydata -> dstoffset
		//3. obtain a new entry in lut_compression -> idcompression
		//4. copy the [zptr,zptr+zbytes] into in allmydata, starting from dstoffset
		//5. for all blocks, set the myblockindices[blockids[i]] to a valid state
		//6. set nblocks to zero	
		
		const unsigned char * const zptr = inputbuffer;
		size_t zbytes = bufsize;
		int dstoffset = -1;
		int idcompression = -1;
		
		//1.
		{
			z_stream myzstream = {0};
			deflateInit(&myzstream, Z_DEFAULT_COMPRESSION);
			
			unsigned mah = maxsize;	
			int err = deflate_inplace(&myzstream, inputbuffer, bufsize, &mah);

			assert(err == Z_OK);
			
			zbytes = myzstream.total_out;
		}
				
		//2-3.
		#pragma critical
		{
			dstoffset = written_bytes;
			written_bytes += zbytes;

			idcompression = lut_compression.size();
			lut_compression.push_back(zbytes);
		}
		
		//4.
		assert(allmydata.size() >= written_bytes);
		memcpy(&allmydata.front(), zptr, zbytes);

		//5.
		for(int i = 0; i < nblocks; ++i)
		{
			const int entry = metablocks[i].idcompression;
			assert(entry >= 0 && entry < myblockindices.size());
			
			myblockindices[entry] = metablocks[i];
			myblockindices[entry].idcompression = idcompression;
		}
		
		//6.
		bufsize = 0;
		nblocks = 0;
	}
	
	void _compress(const vector<BlockInfo>& vInfo, const int NBLOCKS, IterativeStreamer streamer, const int channelid)
	{
#pragma omp parallel  
		{
			enum //some considerations about the per-thread working set 
			{ 
				DESIREDMEM = (4 * 1024) * 1024,
				ENTRYSIZE = sizeof(WaveletCompressor) + sizeof(int),
				ENTRIES_CANDIDATE = DESIREDMEM / ENTRYSIZE,
				ENTRIES = ENTRIES_CANDIDATE ? ENTRIES_CANDIDATE : 1,
				BUFFERSIZE = ENTRIES * ENTRYSIZE, //sonnentanz
				ALERT = (ENTRIES - 1) * ENTRYSIZE //got rain?
			};
						
			BlockMetadata hotblocks[ENTRIES];
			unsigned char compressedbuffer[BUFFERSIZE];

			int mybytes = 0, myhotblocks = 0;
			
			Timer timer;
			timer.start();
			
#pragma omp for
			for(int i = 0; i < NBLOCKS; ++i)
			{
				//wavelet compression
				{
					FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
					
					Real mysoabuffer[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
				
					for(int iz=0; iz<FluidBlock::sizeZ; iz++)
						for(int iy=0; iy<FluidBlock::sizeY; iy++)
							for(int ix=0; ix<FluidBlock::sizeX; ix++)
								mysoabuffer[iz][iy][ix] = streamer.operate(channelid, b(ix, iy, iz));
							
					WaveletCompressor compressor;

					//wavelet digestion
					const int nbytes = (int)compressor.compress(this->threshold, this->halffloat, mysoabuffer);
					memcpy(compressedbuffer + mybytes, &nbytes, sizeof(nbytes));
					mybytes += sizeof(nbytes);
					
					memcpy(compressedbuffer + mybytes, compressor.data(), sizeof(unsigned char) * nbytes);
					mybytes += nbytes;
				}
				
				//building the meta data
				{
					BlockMetadata curr = { i, myhotblocks, vInfo[i].index[0], vInfo[i].index[1], vInfo[i].index[2]};
					hotblocks[myhotblocks] = curr;
					myhotblocks++;
				}
				
				if (mybytes >= ALERT || myhotblocks >= ENTRIES)
					_zcompress_and_flush(compressedbuffer, mybytes, BUFFERSIZE, hotblocks, myhotblocks);
			}
			
			if (mybytes > 0)
				_zcompress_and_flush(compressedbuffer, mybytes, BUFFERSIZE, hotblocks, myhotblocks);
			
			workload[omp_get_thread_num()] += timer.stop();
		}
	}
	
	string _header(GridType & inputGrid)
	{
		const int xtotalbpd = inputGrid.getBlocksPerDimension(0);
		const int ytotalbpd = inputGrid.getBlocksPerDimension(1);
		const int ztotalbpd = inputGrid.getBlocksPerDimension(2);
		
		std::stringstream ss;
			
		ss << "\n==============START-ASCI-HEADER==============\n";

		{
			int one = 1;
			bool isone = *(char *)(&one);
			
			ss << "Endianess: " << (isone ? "little" : "big") << "\n";
		}
		
		ss << "sizeofReal: " << sizeof(Real) << "\n";
		ss << "Blocks: " << xtotalbpd << " x "  << ytotalbpd << " x " << ztotalbpd  << "\n";
		ss << "HalfFloat: " << (this->halffloat ? "yes" : "no") << "\n";
		ss << "==============START-BINARY-LUT==============\n";
		
		return ss.str();
	}
	
	void _write(GridType & inputGrid, string fileName, IterativeStreamer streamer, const int channel)
	{				
		const vector<BlockInfo> infos = inputGrid.getBlocksInfo();
		const int NBLOCKS = infos.size();
			
		//compress my data, prepare for serialization
		{			
			written_bytes = 0;
			
			if (allmydata.size() == 0)
				allmydata.resize(NBLOCKS * sizeof(Real) * NPTS);
			
			myblockindices.clear();
			myblockindices.resize(NBLOCKS);
			
			lut_compression.clear();

			_compress(infos, infos.size(), streamer, channel);

			lut_compression.insert(lut_compression.begin(), lut_compression.size()); //insert nchunks at the beginning
			lut_compression.push_back(allmydata.size()); //add the "#nchunk+1 start"
		}
		
		const MPI::Comm& mycomm = inputGrid.getCartComm();
		
		const size_t mygid = mycomm.Get_rank();
		
		string binaryocean_title = "\n==============START-BINARY-OCEAN==============\n";
		
		MPI::File myfile = MPI::File::Open(mycomm, fileName.c_str(),  MPI::MODE_WRONLY|MPI::MODE_CREATE, MPI::INFO_NULL);
		
		size_t current_displacement = 0;

		//write the mini-header
		{			
			size_t blank_address = -1;

			const int miniheader_bytes = sizeof(blank_address) + binaryocean_title.size();			
			
			current_displacement += miniheader_bytes;
		}
				
		//write the buffer - alias the binary ocean 
		{
			myfile.Seek_shared(current_displacement, MPI_SEEK_CUR);

			myfile.Write_ordered(&allmydata.front(), written_bytes, MPI_CHAR);
			
			current_displacement += myfile.Get_position_shared();			
		}
		
		//go back at the blank address and fill it with the displacement
		if (mygid == 0)
		{
			myfile.Write(&current_displacement, sizeof(current_displacement), MPI_CHAR);
			myfile.Write(binaryocean_title.c_str(), binaryocean_title.size(), MPI_CHAR);				
		}
		
		//write the header
		{
			const string header = _header(inputGrid);
			const size_t header_bytes = header.length();
			
			if (mygid == 0)
				myfile.Write_at(current_displacement, header.c_str(), header_bytes, MPI_CHAR);
			
			current_displacement += header_bytes;			
		}
			
		//write block metadata
		{
			const int metadata_bytes = myblockindices.size() * sizeof(BlockMetadata);
			
			myfile.Write_at_all(current_displacement + mygid * metadata_bytes, &myblockindices.front(), metadata_bytes, MPI_CHAR);
			
			current_displacement += metadata_bytes * mycomm.Get_size();			
		}
		
		//write the lut header
		{			
			string title = "\n==============START-BINARY-LUT==============\n";
			
			const int title_bytes = title.size();
			
			if (mygid == 0)
				myfile.Write_at(current_displacement, title.c_str(), title_bytes, MPI_CHAR);				
			
			current_displacement += title_bytes;
		}
		
		//write the local buffer entries 
		{
			myfile.Seek_shared(current_displacement, MPI_SEEK_CUR);
			
			const int lut_bytes = lut_compression.size() * sizeof(size_t);
		
			myfile.Write_ordered( &lut_compression.front(), lut_bytes, MPI_CHAR);
		}
		
		myfile.Close();
		
		//just a report now
		if (verbosity)
		{
			vector<float> ts = workload;
			
			sort(ts.begin(), ts.end());
			
			float tmin = ts.front();
			float tmed = ts[ts.size() / 2];
			float tmax = ts.back();

			size_t aggregate_written_bytes = -1;
			
			mycomm.Reduce(&written_bytes, &aggregate_written_bytes, 1, MPI_UINT64_T, MPI::SUM, 0);
			mycomm.Reduce(mygid ? &tmin : MPI::IN_PLACE, &tmin, 1, MPI_FLOAT, MPI::MIN, 0);
			mycomm.Reduce(mygid ? &tmed : MPI::IN_PLACE, &tmed, 1, MPI_FLOAT, MPI::SUM, 0);
			mycomm.Reduce(mygid ? &tmax : MPI::IN_PLACE, &tmax, 1, MPI_FLOAT, MPI::MAX, 0);
			
			tmed /= mycomm.Get_size();
			
			if (mygid == 0)
			{
				printf("Channel %d: %.2f kB, wavelet-threshold: %.1e, compr. rate: %.2f, imbalance (min, avg-median, max): %.2f %.2f %.2f\n", 
					channel, aggregate_written_bytes/1024., 
					threshold, NPTS * sizeof(Real) * NBLOCKS * mycomm.Get_size() / (float) aggregate_written_bytes, 
					tmin, tmed, tmax);
			}
		}
	}
	
public:
	
	void set_threshold(const Real threshold) 
	{ 
		//printf("setting the threshold to %e\n", threshold);
		this->threshold = threshold; 
	}

	void float16() { halffloat = true; }
	
	void verbose() { verbosity = true; }
	
	SerializerIO_WaveletCompression_MPI_Simple(): 
	threshold(0), halffloat(false), verbosity(false), workload(omp_get_max_threads()), written_bytes(0){ }
		
	void Write(GridType & inputGrid, string fileName, IterativeStreamer streamer = IterativeStreamer())
	{		
		for(int channel = 0; channel < NCHANNELS; ++channel)
		{
			std::stringstream ss;
			ss << "." << streamer.name() << ".channel"  << channel;
			
			_write(inputGrid, fileName + ss.str(), streamer, channel);
		}
	}
};

/* THIS CODE SERVES US TO COMPRESS IN-PLACE. TAKEN FROM THE WEB
 * http://stackoverflow.com/questions/12398377/is-it-possible-to-have-zlib-read-from-and-write-to-the-same-memory-buffer
 * 
 * Compress buf[0..len-1] in place into buf[0..*max-1].  *max must be greater
   than or equal to len.  Return Z_OK on success, Z_BUF_ERROR if *max is not
   enough output space, Z_MEM_ERROR if there is not enough memory, or
   Z_STREAM_ERROR if *strm is corrupted (e.g. if it wasn't initialized or if it
   was inadvertently written over).  If Z_OK is returned, *max is set to the
   actual size of the output.  If Z_BUF_ERROR is returned, then *max is
   unchanged and buf[] is filled with *max bytes of uncompressed data (which is
   not all of it, but as much as would fit).

   Incompressible data will require more output space than len, so max should
   be sufficiently greater than len to handle that case in order to avoid a
   Z_BUF_ERROR. To assure that there is enough output space, max should be
   greater than or equal to the result of deflateBound(strm, len).

   strm is a deflate stream structure that has already been successfully
   initialized by deflateInit() or deflateInit2().  That structure can be
   reused across multiple calls to deflate_inplace().  This avoids unnecessary
   memory allocations and deallocations from the repeated use of deflateInit()
   and deflateEnd(). */
inline int deflate_inplace(z_stream *strm, unsigned char *buf, unsigned len,
                    unsigned *max)
{
    int ret;                    /* return code from deflate functions */
    unsigned have;              /* number of bytes in temp[] */
    unsigned char *hold;        /* allocated buffer to hold input data */
    unsigned char temp[11];     /* must be large enough to hold zlib or gzip
                                   header (if any) and one more byte -- 11
                                   works for the worst case here, but if gzip
                                   encoding is used and a deflateSetHeader()
                                   call is inserted in this code after the
                                   deflateReset(), then the 11 needs to be
                                   increased to accomodate the resulting gzip
                                   header size plus one */
                                   

    /* initialize deflate stream and point to the input data */
    ret = deflateReset(strm);
    
    if (ret != Z_OK)
        return ret;
    strm->next_in = buf;
    strm->avail_in = len;

    /* kick start the process with a temporary output buffer -- this allows
       deflate to consume a large chunk of input data in order to make room for
       output data there */
    if (*max < len)
        *max = len;
    strm->next_out = temp;
    strm->avail_out = sizeof(temp) > *max ? *max : sizeof(temp);
    ret = deflate(strm, Z_FINISH);
          
    if (ret == Z_STREAM_ERROR)
        return ret;

    /* if we can, copy the temporary output data to the consumed portion of the
       input buffer, and then continue to write up to the start of the consumed
       input for as long as possible */
    have = strm->next_out - temp;
    if (have <= (strm->avail_in ? len - strm->avail_in : *max)) {
        memcpy(buf, temp, have);
        strm->next_out = buf + have;
        have = 0;
        while (ret == Z_OK) {
            strm->avail_out = strm->avail_in ? strm->next_in - strm->next_out :
                                               (buf + *max) - strm->next_out;
            ret = deflate(strm, Z_FINISH);
        }
        if (ret != Z_BUF_ERROR || strm->avail_in == 0) {
            *max = strm->next_out - buf;
            return ret == Z_STREAM_END ? Z_OK : ret;
        }
    }

    /* the output caught up with the input due to insufficiently compressible
       data -- copy the remaining input data into an allocated buffer and
       complete the compression from there to the now empty input buffer (this
       will only occur for long incompressible streams, more than ~20 MB for
       the default deflate memLevel of 8, or when *max is too small and less
       than the length of the header plus one byte) */

    hold = (unsigned char*)strm->zalloc(strm->opaque, strm->avail_in, 1);
    if (hold == Z_NULL)
        return Z_MEM_ERROR;
    memcpy(hold, strm->next_in, strm->avail_in);
    strm->next_in = hold;
    if (have) {
        memcpy(buf, temp, have);
        strm->next_out = buf + have;
    }
    strm->avail_out = (buf + *max) - strm->next_out;
    ret = deflate(strm, Z_FINISH);
    strm->zfree(strm->opaque, hold);
    *max = strm->next_out - buf;
    return ret == Z_OK ? Z_BUF_ERROR : (ret == Z_STREAM_END ? Z_OK : ret);
}

