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
inline size_t zdecompress(unsigned char * inputbuf, size_t ninputbytes, unsigned char * outputbuf, const size_t maxsize);

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
		#pragma omp critical
		{
			dstoffset = written_bytes;
			written_bytes += zbytes;

			idcompression = lut_compression.size();
			lut_compression.push_back(dstoffset);
		}
		
		//4.
		assert(allmydata.size() >= written_bytes);
		memcpy(&allmydata.front() + dstoffset, zptr, zbytes);
		
		//5.
		for(int i = 0; i < nblocks; ++i)
		{
			const int entry = metablocks[i].idcompression;
			assert(entry >= 0 && entry < myblockindices.size());
			
			myblockindices[entry] = metablocks[i];
			myblockindices[entry].idcompression = idcompression;
			myblockindices[entry].subid = i;
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
				DESIREDMEM = (1 * 1024) * 1024,
				ENTRYSIZE = sizeof(WaveletCompressor) + sizeof(int),
				ENTRIES_CANDIDATE = DESIREDMEM / ENTRYSIZE,
				ENTRIES = ENTRIES_CANDIDATE ? ENTRIES_CANDIDATE : 1,
				BUFFERSIZE = ENTRIES * ENTRYSIZE, //sonnentanz
				ALERT = (ENTRIES - 1) * ENTRYSIZE //got rain?
			};
						
			BlockMetadata hotblocks[ENTRIES];
			unsigned char compressedbuffer[BUFFERSIZE];
			//int * header = (int *)compressedbuffer;

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
				
				//*header = -1;//lets fake it
				
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
		
		const int xbpd = inputGrid.getResidentBlocksPerDimension(0);
		const int ybpd = inputGrid.getResidentBlocksPerDimension(1);
		const int zbpd = inputGrid.getResidentBlocksPerDimension(2);
		
		std::stringstream ss;
			
		ss << "\n==============START-ASCI-HEADER==============\n";

		{
			int one = 1;
			bool isone = *(char *)(&one);
			
			ss << "Endianess: " << (isone ? "little" : "big") << "\n";
		}
		
		ss << "sizeofReal: " << sizeof(Real) << "\n";
		ss << "Blocks: " << xtotalbpd << " x "  << ytotalbpd << " x " << ztotalbpd  << "\n";
		ss << "SubdomainBlocks: " << xbpd << " x "  << ybpd << " x " << zbpd  << "\n";
		ss << "HalfFloat: " << (this->halffloat ? "yes" : "no") << "\n";
		ss << "==============START-BINARY-METABLOCKS==============\n";
		
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

			const int nchunks = lut_compression.size();
			lut_compression.insert(lut_compression.begin(), nchunks); //insert nchunks at the beginning
			lut_compression.push_back(written_bytes); //add the "#nchunk+1 start"
			
		}
		
		const MPI::Intracomm& mycomm = inputGrid.getCartComm();
		
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
			myfile.Seek_shared(current_displacement, MPI_SEEK_SET);

			myfile.Write_ordered(&allmydata.front(), written_bytes, MPI_CHAR);
			
			current_displacement = myfile.Get_position_shared();			
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
			myfile.Seek_shared(current_displacement, MPI_SEEK_SET);
			
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
	
	struct CompressedBlock{ size_t start, extent; int subid; } ;
	
	void _read(string path)
	{
		//THE FIRST PART IS SEQUENTIAL
		//THE SECOND ONE IS RANDOM ACCESS
		
		size_t global_header_displacement = -1;
		int NBLOCKS = -1;
		int totalbpd[3] = {-1, -1, -1};
		int bpd[3] = { -1, -1, -1};
		string binaryocean_title = "\n==============START-BINARY-OCEAN==============\n";	
		const int miniheader_bytes = sizeof(size_t) + binaryocean_title.size();		
		
		vector<BlockMetadata> metablocks;
		
		//random access data structures: meta2subchunk, lutchunks;
		map<int, map<int, map<int, CompressedBlock > > > meta2subchunk;
		vector<size_t> lutchunks;
				
		//printf("reading %s\n", path.c_str());
		
		{
			FILE * file = fopen(path.c_str(), "rb");
		
			assert(file);
			
			//reading the header and mini header
			{
				size_t header_displacement = -1;
				fread(&header_displacement, sizeof(size_t), 1, file);
				
				fseek(file, header_displacement, SEEK_SET);
				global_header_displacement = header_displacement;
			
				char buf[1024];
				fgets(buf, sizeof(buf), file);
				fgets(buf, sizeof(buf), file);
				assert(string("==============START-ASCI-HEADER==============\n") == string(buf));
			
				fscanf(file, "Endianess:  %s\n", buf);
				assert(string(buf) == "little");
			
				int sizeofreal = -1;
				fscanf(file, "sizeofReal:  %d\n", &sizeofreal);
				assert(sizeof(Real) == sizeofreal);		
			
				fscanf(file, "Blocks: %d x %d x %d\n", totalbpd, totalbpd + 1, totalbpd + 2);
				fscanf(file, "SubdomainBlocks: %d x %d x %d\n", bpd, bpd + 1, bpd + 2);

				fscanf(file, "HalfFloat: %s\n", buf);
				this->halffloat = (string(buf) == "yes");
					
				fgets(buf, sizeof(buf), file);
				assert(string("==============START-BINARY-METABLOCKS==============\n") == string(buf));
				
				NBLOCKS = totalbpd[0] * totalbpd[1] * totalbpd[2];
				//printf("Blocks: %d -> %dx%dx%d -> subdomains of %dx%dx%d\n", 
				//	NBLOCKS, totalbpd[0], totalbpd[1], totalbpd[2], bpd[0], bpd[1], bpd[2]);
				
			}
			
			//reading the binary lut
			{
				metablocks.resize(NBLOCKS);
				
				for(int i = 0; i < NBLOCKS; ++i)
				{
					BlockMetadata entry;
					fread(&entry, sizeof(entry), 1, file);
					//printf("reading metablock %d -> %d %d %d  cid %d\n", i, entry.ix, entry.iy, entry.iz, entry.idcompression); 
					assert(entry.idcompression >= 0 && entry.idcompression < bpd[0] * bpd[1] * bpd[2]);
					metablocks[i] = entry;
				}
			}
			
			//reading the compression lut
			{				
				char buf[1024];
				fgetc(file);//buf, sizeof(buf), file);
				fgets(buf, sizeof(buf), file);
				assert(string("==============START-BINARY-LUT==============\n") == string(buf));

				//printf("reading compression lut..at location %d\n", ftell(file));
				
				bool done = false;
				
				size_t base = miniheader_bytes;
				
				const int BPS = bpd[0] * bpd[1] * bpd[2];
				assert(NBLOCKS % BPS == 0);
				const int SUBDOMAINS = NBLOCKS / BPS;
				
				for(int s = 0, currblock = 0; s < SUBDOMAINS; ++s)
				{
					const int nglobalchunks = lutchunks.size();
					
					assert(!feof(file));
					
					size_t nchunks = -1;
					fread(&nchunks, sizeof(nchunks), 1, file);
					//printf("this buffer has %d entries\n", nchunks);

					vector<size_t> start_lut(nchunks + 1);
					fread(&start_lut.front(), sizeof(size_t), start_lut.size(), file);
					const size_t myamount = start_lut.back();
					printf("my amount is %d\n", myamount);
					start_lut.pop_back();
					
					//compute the global positioning of the compressed chunks within the file				
					for(int i = 0; i < start_lut.size(); ++i)
					{
						assert(start_lut[i] < myamount);
						start_lut[i] += base;
					}
						
					assert(myamount > 0);
					base += myamount;
					assert(base <= global_header_displacement);
					
					//compute the base for this blocks
					for(int i = 0; i < BPS; ++i, ++currblock)
						metablocks[currblock].idcompression += nglobalchunks;
					
					lutchunks.insert(lutchunks.end(), start_lut.begin(), start_lut.end());					
				} 
				
				//printf("minheader takes %d bytes\n", miniheader_bytes);
				//printf("my base is now 0x%x whereas my header is at 0x%x -> discrepancy is %d bytes\n", base, global_header_displacement, global_header_displacement - base);
				assert(base == global_header_displacement);
				lutchunks.push_back(base);
			
				{
					int c = fgetc(file);
					
					do 
					{ 
						//printf("shouldnt be here! 0x%x\n", c); 
						c = fgetc(file);
					}
					while (! feof(file) );
				}			
			}
			
			fclose(file);
		}
		
		for(int i = 0; i < NBLOCKS ; ++i)
		{
			BlockMetadata entry = metablocks[i];
			
			assert(entry.idcompression >= 0);
			assert(entry.idcompression < lutchunks.size()-1);
			
			size_t start_address = lutchunks[entry.idcompression];
			size_t end_address = lutchunks[entry.idcompression + 1];
			
			assert(start_address < end_address);
			assert(end_address <= global_header_displacement);
			assert( start_address < global_header_displacement );
			
			CompressedBlock compressedblock = { start_address, end_address - start_address, entry.subid };
			
			meta2subchunk[entry.iz][entry.iy][entry.ix] = compressedblock;
		}
		
		//at this point we can access every block
		{
			//lets try with block 7, 3, 5
			FILE * f = fopen(path.c_str(), "rb");
			assert(f);
			
			const int ix = 0;
			const int iy = 2;
			const int iz = 0;
			
			assert(ix >= 0 && ix < totalbpd[0]);
			assert(iy >= 0 && iy < totalbpd[1]);
			assert(iz >= 0 && iz < totalbpd[2]);
			
			assert(meta2subchunk.find(iz) != meta2subchunk.end());
			assert(meta2subchunk[iz].find(iy) != meta2subchunk[iz].end());
			assert(meta2subchunk[iz][iy].find(ix) != meta2subchunk[iz][iy].end());
			
			CompressedBlock compressedchunk = meta2subchunk[2][0][3];
			
			size_t start = compressedchunk.start;
			assert(start >= miniheader_bytes);
			assert(start < global_header_displacement);
			assert(start + compressedchunk.extent < global_header_displacement);
			
			vector<unsigned char> compressedbuf(compressedchunk.extent);
			fseek(f, compressedchunk.start, SEEK_SET);
			fread(&compressedbuf.front(), compressedchunk.extent, 1, f);
			assert(!feof(f));
			//printf("my steart: 0x%x extent is %d my subid %d\n" , compressedchunk.start, compressedchunk.extent, compressedchunk.subid);

			
			vector<unsigned char> waveletbuf(4 << 20);
			const size_t decompressedbytes = zdecompress(&compressedbuf.front(), compressedbuf.size(), &waveletbuf.front(), waveletbuf.size());
			//printf("decompressed bytes is %d\n", decompressedbytes);
			int readbytes = 0;
			for(int i = 0; i<compressedchunk.subid; ++i)
			{
				int nbytes = *(int *)&waveletbuf[readbytes];
				readbytes += sizeof(int);
				assert(readbytes <= decompressedbytes);
				//printf("scanning nbytes...%d\n", nbytes);
				readbytes += nbytes;
								assert(readbytes <= decompressedbytes);

			}
			
			Real MYBLOCK[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
			
			{
				int nbytes = *(int *)&waveletbuf[readbytes];
				readbytes += sizeof(int);
				assert(readbytes <= decompressedbytes);

				//printf("OK MY BYTES ARE: %d\n", nbytes);
				
				WaveletCompressor compressor;
				memcpy(compressor.data(), &waveletbuf[readbytes], nbytes);
				readbytes += nbytes;

				compressor.decompress(halffloat, nbytes, MYBLOCK);
			}
			
			printf("OK FINAL TEST: THE DATA\n");
			for(int iz = 0; iz< _BLOCKSIZE_; ++iz)
				for(int iy = 0; iy< _BLOCKSIZE_; ++iy)
					for(int ix = 0; ix< _BLOCKSIZE_; ++ix)
						printf("%d %d %d: %e\n", ix, iy, iz, MYBLOCK[iz][iy][ix]);
				
			fclose(f);
		}
	}
	
public:
	
	void set_threshold(const Real threshold) { this->threshold = threshold; }

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
	
	void Read(string fileName, IterativeStreamer streamer = IterativeStreamer())
	{
		for(int channel = 0; channel < NCHANNELS; ++channel)
		{
			std::stringstream ss;
			ss << "." << streamer.name() << ".channel"  << channel;
			
			_read(fileName + ss.str());
		}
	}
};

inline size_t zdecompress(unsigned char * inputbuf, size_t ninputbytes, unsigned char * outputbuf, const size_t maxsize)
{
	int decompressedbytes = 0;

	z_stream datastream = {0};
	datastream.total_in = datastream.avail_in = ninputbytes;
	datastream.total_out = datastream.avail_out = maxsize;
	datastream.next_in = inputbuf;
	datastream.next_out = outputbuf;

	const int retval = inflateInit(&datastream);

	if (retval == Z_OK && inflate(&datastream, Z_FINISH))
	{
		printf("goooooood\n");
		decompressedbytes = datastream.total_out;
	}
	else
	{
			printf("ZLIB DECOMPRESSION FAILURE!!\n");
			abort();
	}

	inflateEnd(&datastream);
	printf("ADSDSO %d\n", decompressedbytes);
	return decompressedbytes;
}

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

