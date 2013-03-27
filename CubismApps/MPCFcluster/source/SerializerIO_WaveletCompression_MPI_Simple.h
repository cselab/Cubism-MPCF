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
#include <numeric>
#include <omp.h>

using namespace std;

#include "WaveletSerializationTypes.h"
#include "FullWaveletTransform.h"
#include "WaveletCompressor.h"

#include "CompressionEncoders.h"

template<typename GridType, typename IterativeStreamer>
class SerializerIO_WaveletCompression_MPI_SimpleBlocking
{
protected:
	
	enum 
	{
		NCHANNELS = IterativeStreamer::channels,
		NPTS = _BLOCKSIZE_ * _BLOCKSIZE_ * _BLOCKSIZE_
	};
	
	string binaryocean_title, binarylut_title, header;
	
	vector< BlockMetadata > myblockindices; //tells in which compressed chunk is any block, nblocks
	HeaderLUT lutheader;
	vector< size_t > lut_compression; //tells the number of compressed chunk, and where do they start, nchunks + 2
	vector< unsigned char > allmydata; //buffer with the compressed data
	size_t written_bytes, pending_writes, completed_writes;
	
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
			
			//exception: we have to resize allmydata
			if (written_bytes > allmydata.size())
			{
				//spin-wait until writes complete
				while (pending_writes != completed_writes);
				
				//safely resize
				allmydata.resize(written_bytes);
			}
			
			idcompression = lut_compression.size();
			lut_compression.push_back(dstoffset);
			
			++pending_writes;
		}
		
		//4.
		assert(allmydata.size() >= written_bytes);
		memcpy(&allmydata.front() + dstoffset, zptr, zbytes);
		
#pragma omp atomic
		++completed_writes;
		
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
				DESIREDMEM = (2 * 1024) * 1024,
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
	
	virtual void _to_file(const MPI::Intracomm& mycomm, const string fileName)
	{
		const int mygid = mycomm.Get_rank();
		const int nranks = mycomm.Get_size();
		
		MPI::Info myfileinfo = MPI::Info::Create();
		
		myfileinfo.Set("access_syle", "write_once");
		
		MPI::File myfile = MPI::File::Open(mycomm, fileName.c_str(),  MPI::MODE_WRONLY | MPI::MODE_CREATE, myfileinfo);
		
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
			const size_t header_bytes = header.size();
			
			if (mygid == 0)
				myfile.Write_at(current_displacement, header.c_str(), header_bytes, MPI_CHAR);
			
			current_displacement += header_bytes;			
		}
		
		//write block metadata
		{
			const int metadata_bytes = myblockindices.size() * sizeof(BlockMetadata);
			
			myfile.Write_at_all(current_displacement + mygid * metadata_bytes, &myblockindices.front(), metadata_bytes, MPI_CHAR);
			
			current_displacement += metadata_bytes * nranks;			
		}
		
		//write the lut title
		{			
			const int title_bytes = binarylut_title.size();
			
			if (mygid == 0)
				myfile.Write_at(current_displacement, binarylut_title.c_str(), title_bytes, MPI_CHAR);				
			
			current_displacement += title_bytes;
		}
		
		//write the local buffer entries 
		{			
			assert(lut_compression.size() == 0);
			
			const int lutheader_bytes = sizeof(lutheader);
			
			myfile.Write_at_all(current_displacement + mygid * lutheader_bytes, &lutheader, lutheader_bytes, MPI_CHAR);
		}
		
		myfile.Close(); //bon voila tu vois ou quoi
	}
	
	void _write(GridType & inputGrid, string fileName, IterativeStreamer streamer, const int channel)
	{				
		const vector<BlockInfo> infos = inputGrid.getBlocksInfo();
		const int NBLOCKS = infos.size();
		
		//prepare the headers
		{
			this->binaryocean_title = "\n==============START-BINARY-OCEAN==============\n";
			this->binarylut_title = "\n==============START-BINARY-LUT==============\n";
			
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
				ss << "Blocksize: " << _BLOCKSIZE_ << "\n";
				ss << "Blocks: " << xtotalbpd << " x "  << ytotalbpd << " x " << ztotalbpd  << "\n";
				ss << "SubdomainBlocks: " << xbpd << " x "  << ybpd << " x " << zbpd  << "\n";
				ss << "HalfFloat: " << (this->halffloat ? "yes" : "no") << "\n";
				ss << "Wavelets: " << WaveletsOnInterval::ChosenWavelets_GetName() << "\n";
				ss << "Encoder: " << "zlib" << "\n";
				ss << "==============START-BINARY-METABLOCKS==============\n";
				
				this->header = ss.str();
			}
		}
		
		//compress my data, prepare for serialization
		{			
			written_bytes = pending_writes = completed_writes = 0;
			
			if (allmydata.size() == 0)
			{
				const size_t speculated_compression_rate = 10;
				allmydata.resize(NBLOCKS * sizeof(Real) * NPTS);
			}
			
			myblockindices.clear();
			myblockindices.resize(NBLOCKS);
			
			lut_compression.clear();
			
			_compress(infos, infos.size(), streamer, channel);
			
			//manipulate the file data (allmydata, lut_compression, myblockindices)
			//so that they are file-friendly
			{
				const int nchunks = lut_compression.size();
				const size_t extrabytes = lut_compression.size() * sizeof(size_t);
				const char * const lut_ptr = (char *)&lut_compression.front();
				
				allmydata.insert(allmydata.begin() + written_bytes, lut_ptr, lut_ptr + extrabytes);
				lut_compression.clear();
				
				HeaderLUT newvalue = { written_bytes + extrabytes, nchunks };
				lutheader = newvalue;
				
				written_bytes += extrabytes;
			}
		}
		
		const MPI::Comm& mycomm = inputGrid.getCartComm();
		const size_t mygid = mycomm.Get_rank();
		
		//write into the file
		_to_file(mycomm, fileName);
		
		//just a report now
		if (verbosity)
		{			
			float tmin = *std::min_element(workload.begin(), workload.end());
			float tmax = *std::max_element(workload.begin(), workload.end());
			float tmed = std::accumulate(workload.begin(), workload.end(), 0) / workload.size();
			
			size_t aggregate_written_bytes = -1;
			
			mycomm.Reduce(&written_bytes, &aggregate_written_bytes, 1, MPI_UINT64_T, MPI::SUM, 0);
			mycomm.Reduce(mygid ? &tmin : MPI::IN_PLACE, &tmin, 1, MPI_FLOAT, MPI::MIN, 0);
			mycomm.Reduce(mygid ? &tmed : MPI::IN_PLACE, &tmed, 1, MPI_FLOAT, MPI::SUM, 0);
			mycomm.Reduce(mygid ? &tmax : MPI::IN_PLACE, &tmax, 1, MPI_FLOAT, MPI::MAX, 0);
			
			tmed /= mycomm.Get_size();
			
			if (mygid == 0)
			{
				printf("Channel %d: %.2f kB, wavelet-threshold: %.1e, compr. rate: %.2f, thread-workload (min, avg, max): %.2f %.2f %.2f\n", 
					   channel, aggregate_written_bytes/1024., 
					   threshold, NPTS * sizeof(Real) * NBLOCKS * mycomm.Get_size() / (float) aggregate_written_bytes, 
					   tmin, tmed, tmax);
			}
		}
	}
		
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
				
				int bsize = -1;
				fscanf(file, "Blocksize: %d\n", &bsize);
				assert(bsize == _BLOCKSIZE_);
				fscanf(file, "Blocks: %d x %d x %d\n", totalbpd, totalbpd + 1, totalbpd + 2);
				fscanf(file, "SubdomainBlocks: %d x %d x %d\n", bpd, bpd + 1, bpd + 2);
				
				fscanf(file, "HalfFloat: %s\n", buf);
				this->halffloat = (string(buf) == "yes");
				
				fscanf(file, "Wavelets: %s\n", buf);
				assert(buf == string(WaveletsOnInterval::ChosenWavelets_GetName()));
				
				fscanf(file, "Encoder: %s\n", buf);
				assert(buf == string("zlib"));
				
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
				
				vector<HeaderLUT> headerluts(SUBDOMAINS); //oh mamma mia
				fread(&headerluts.front(), sizeof(HeaderLUT), SUBDOMAINS, file);
				
				for(int s = 0, currblock = 0; s < SUBDOMAINS; ++s)
				{
					const int nglobalchunks = lutchunks.size();
					
					assert(!feof(file));
					
					/* this is history
					 size_t nchunks = -1;
					 fread(&nchunks, sizeof(nchunks), 1, file);
					 //printf("this buffer has %d entries\n", nchunks);
					 
					 vector<size_t> mylut(nchunks + 1);
					 fread(&mylut.front(), sizeof(size_t), mylut.size(), file);
					 const size_t myamount = mylut.back();
					 printf("my amount is %d\n", myamount);
					 mylut.pop_back();
					 */
					
					const int nchunks = headerluts[s].nchunks;
					const size_t myamount = headerluts[s].aggregate_bytes;
					const size_t lutstart = base + myamount - sizeof(size_t) * nchunks;
					//printf("my header lut is %d %d (0x%x 0x%x)\n", headerluts[s].aggregate_bytes, headerluts[s].nchunks,
					//	   headerluts[s].aggregate_bytes, headerluts[s].nchunks);
					//read the lut
					fseek(file, lutstart, SEEK_SET);
					vector<size_t> mylut(nchunks);
					fread(&mylut.front(), sizeof(size_t), nchunks, file);
					
					//for(int i=0; i< nchunks; ++i)
					//	printf("reading %d) 0x%x\n", i, mylut[i]);
					
					for(int i=0; i< nchunks; ++i)
						assert(mylut[i] < myamount);
					
					for(int i=1; i< nchunks; ++i)
						assert(mylut[i-1] < mylut[i]);
					
					//printf("goodu p to here\n");
					
					//bump the lut with the base
					//for(int i=0; i< nchunks; ++i)
					//	mylut[i] += base;
					
					//exit(0);
					//compute the global positioning of the compressed chunks within the file				
					for(int i = 0; i < mylut.size(); ++i)
					{
						assert(mylut[i] < myamount);
						mylut[i] += base;
					}
					
					assert(myamount > 0);
					base += myamount;
					//printf("new base is 0x%x\n", base);
					//printf("global_header_displacement is 0x%x\n", global_header_displacement);
					assert(base <= global_header_displacement);
					
					//compute the base for this blocks
					for(int i = 0; i < BPS; ++i, ++currblock)
						metablocks[currblock].idcompression += nglobalchunks;
					
					lutchunks.insert(lutchunks.end(), mylut.begin(), mylut.end());					
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
			
			CompressedBlock compressedchunk = meta2subchunk[iz][iy][ix];
			
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
	
	SerializerIO_WaveletCompression_MPI_SimpleBlocking(): 
	threshold(0), halffloat(false), verbosity(false), workload(omp_get_max_threads()), 
	written_bytes(0), pending_writes(0)
	{
	}
	
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

template<typename GridType, typename IterativeStreamer>
class SerializerIO_WaveletCompression_MPI_Simple : public SerializerIO_WaveletCompression_MPI_SimpleBlocking<GridType, IterativeStreamer>
{	
	size_t nofcalls;
	vector<MPI::Request> pending_requests;
	MPI::File myopenfile;
	
	void _wait_all_quiet()
	{
		//wait for pending requests
		assert(pending_requests.size() > 0); //bah ouaie ou quois, tu vois bon voila
		MPI::Request::Waitall(pending_requests.size(), &pending_requests.front());
		pending_requests.clear();
		
		//close the split collective io
		myopenfile.Write_ordered_end(&this->allmydata.front());
		
		//e buonanotte
		myopenfile.Close();
	}
	
	//et bon, la on ce lance le fleurs
	void _to_file(const MPI::Intracomm& mycomm, const string fileName)	
	{
		if (nofcalls)
			_wait_all_quiet(); 
		
		const int mygid = mycomm.Get_rank();
		const int nranks = mycomm.Get_size();
		
		MPI::Info myfileinfo = MPI::Info::Create();
		
		myfileinfo.Set("access_syle", "write_once");
		
		myopenfile = MPI::File::Open(mycomm, fileName.c_str(),  MPI::MODE_WRONLY | MPI::MODE_CREATE, myfileinfo);
		
		size_t current_displacement = 0;
		
		
		//write the mini-header
		{			
			size_t blank_address = -1;
			
			const int miniheader_bytes = sizeof(blank_address) + this->binaryocean_title.size();			
			
			current_displacement += miniheader_bytes;
		}
		
		//write the buffer - alias the binary ocean 
		{
			myopenfile.Seek_shared(current_displacement, MPI_SEEK_SET);
			
			myopenfile.Write_ordered_begin(&this->allmydata.front(), this->written_bytes, MPI_CHAR);
			
			current_displacement = myopenfile.Get_position_shared();			
		}
		
		//go back at the blank address and fill it with the displacement
		if (mygid == 0)
		{
			pending_requests.push_back( myopenfile.Iwrite(&current_displacement, sizeof(current_displacement), MPI_CHAR) );
			pending_requests.push_back( myopenfile.Iwrite(this->binaryocean_title.c_str(), this->binaryocean_title.size(), MPI_CHAR) );				
		}
		
		//write the header
		{
			const size_t header_bytes = this->header.size();
			
			if (mygid == 0)
				pending_requests.push_back( myopenfile.Iwrite_at(current_displacement, this->header.c_str(), header_bytes, MPI_CHAR) );
			
			current_displacement += header_bytes;			
		}
		
		//write block metadata
		{
			const int metadata_bytes = this->myblockindices.size() * sizeof(BlockMetadata);
			
			pending_requests.push_back( myopenfile.Iwrite_at(current_displacement + mygid * metadata_bytes, &this->myblockindices.front(), metadata_bytes, MPI_CHAR) );
			
			current_displacement += metadata_bytes * nranks;			
		}
		
		//write the lut title
		{			
			const int title_bytes = this->binarylut_title.size();
			
			if (mygid == 0)
				pending_requests.push_back( myopenfile.Iwrite_at(current_displacement, this->binarylut_title.c_str(), title_bytes, MPI_CHAR) );				
			
			current_displacement += title_bytes;
		}
		
		//write the local buffer entries 
		{			
			assert(this->lut_compression.size() == 0);
			
			const int lutheader_bytes = sizeof(this->lutheader);
			
			pending_requests.push_back( myopenfile.Iwrite_at(current_displacement + mygid * lutheader_bytes, &this->lutheader, lutheader_bytes, MPI_CHAR) );
		}
		
		++nofcalls;
	}
	
public:
	SerializerIO_WaveletCompression_MPI_Simple(): 
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<GridType, IterativeStreamer>(),
	nofcalls(0)
	{
	}
	
	void force_close() { _wait_all_quiet(); }
	
	~SerializerIO_WaveletCompression_MPI_Simple()
	{
		if (nofcalls)
			_wait_all_quiet();
	}
};
