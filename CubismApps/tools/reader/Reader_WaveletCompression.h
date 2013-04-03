/*
 *  Reader_WaveletCompression.h
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include "../../MPCFnode/source/WaveletCompressor.h"

#include "../../MPCFcluster/source/WaveletSerializationTypes.h"
#include "../../MPCFcluster/source/CompressionEncoders.h"
#include "../../MPCFnode/source/FullWaveletTransform.h"

//MACRO TAKEN FROM http://stackoverflow.com/questions/3767869/adding-message-to-assert
#   define MYASSERT(condition, message) \
do { \
if (! (condition)) { \
std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
<< " line " << __LINE__ << ": " << message << std::endl; \
std::exit(EXIT_FAILURE); \
} \
} while (false)

class Reader_WaveletCompression
{
protected:
	string path;
	
	size_t global_header_displacement;
	int miniheader_bytes;	
	int NBLOCKS;
	int totalbpd[3], bpd[3];
	bool halffloat;
	
	vector<CompressedBlock> idx2chunk;
	
	int _id(int ix, int iy, int iz) const
	{
		assert(ix >= 0 && ix < totalbpd[0]);
		assert(iy >= 0 && iy < totalbpd[1]);
		assert(iz >= 0 && iz < totalbpd[2]);
		
		return ix + totalbpd[0] * ( iy + totalbpd[1] * iz );
	}
	
public:
	
	Reader_WaveletCompression(const string path): NBLOCKS(-1), global_header_displacement(-1), path(path) {	}
	
	virtual void load_file()
	{		
		for(int i = 0; i < 3; ++i)
			totalbpd[i] = -1;
		
		for(int i = 0; i < 3; ++i)
			bpd[i] = -1;
		
		string binaryocean_title = "\n==============START-BINARY-OCEAN==============\n";	
		
		this->miniheader_bytes = sizeof(size_t) + binaryocean_title.size();		
		
		vector<BlockMetadata> metablocks;
		vector<size_t> lutchunks;
		
		{
			FILE * file = fopen(path.c_str(), "rb");
			
			MYASSERT(file, "\nAAATTENZIONE:\nOooops could not open the file. Path: " << path);
			
			//reading the header and mini header
			{
				size_t header_displacement = -1;
				fread(&header_displacement, sizeof(size_t), 1, file);
				
				fseek(file, header_displacement, SEEK_SET);
				global_header_displacement = header_displacement;
				
				char buf[1024];
				fgets(buf, sizeof(buf), file);
				fgets(buf, sizeof(buf), file);
				
				printf("\n%s", buf);
				assert(string("==============START-ASCI-HEADER==============\n") == string(buf));
				
				fscanf(file, "Endianess:  %s\n", buf);
				printf("Endianess: <%s>\n", buf);
				assert(string(buf) == "little");
				
				int sizeofreal = -1;
				fscanf(file, "sizeofReal:  %d\n", &sizeofreal);
				printf("sizeofReal: <%d>\n", sizeofreal);
				
				MYASSERT(sizeof(Real) == sizeofreal, 
						 "\nATTENZIONE:\nSizeof(Real) in the file is " << sizeofreal << " which is wrong\n");		
				
				int bsize = -1;
				fscanf(file, "Blocksize: %d\n", &bsize);
				printf("Blocksize: <%d>\n", bsize);
				
				MYASSERT(bsize == _BLOCKSIZE_,
						 "\nATTENZIONE:\nBlocksize in the file is " << bsize << 
						 " and i have " << _BLOCKSIZE_ << "\n");
				
				fscanf(file, "Blocks: %d x %d x %d\n", totalbpd, totalbpd + 1, totalbpd + 2);
				printf("Blocks: %d x %d x %d\n", totalbpd[0], totalbpd[1], totalbpd[2]);
				
				fscanf(file, "SubdomainBlocks: %d x %d x %d\n", bpd, bpd + 1, bpd + 2);
				printf("SubdomainBlocks: <%d x %d x %d>\n", bpd[0], bpd[1], bpd[2]);
				
				fscanf(file, "HalfFloat: %s\n", buf);
				printf("HalfFloat: <%s>\n", buf);
				this->halffloat = (string(buf) == "yes");
				
				fscanf(file, "Wavelets: %s\n", buf);
				printf("Wavelets: <%s>\n", buf);
				MYASSERT(buf == string(WaveletsOnInterval::ChosenWavelets_GetName()),
						 "\nATTENZIONE:\nWavelets in the file is " << buf << 
						 " and i have " << WaveletsOnInterval::ChosenWavelets_GetName() << "\n");
				
				fscanf(file, "Encoder: %s\n", buf);
				printf("Encoder: <%s>\n", buf);
				MYASSERT(buf == string("zlib"),
						 "\nATTENZIONE:\nWavelets in the file is " << buf << 
						 " and i have zlib.\n");
				
				fgets(buf, sizeof(buf), file);
				
				assert(string("==============START-BINARY-METABLOCKS==============\n") == string(buf));
				printf("==============END ASCI-HEADER==============\n\n");
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
			
			//reading the lut header
			{
				char buf[1024];
				
				fgetc(file);
				fgets(buf, sizeof(buf), file);
				
				assert(string("==============START-BINARY-LUT==============\n") == string(buf));
				
				bool done = false;
				
				size_t base = miniheader_bytes;
				
				const int BPS = bpd[0] * bpd[1] * bpd[2];
				assert(NBLOCKS % BPS == 0);
				const int SUBDOMAINS = NBLOCKS / BPS;
				
				vector<HeaderLUT> headerluts(SUBDOMAINS); //oh mamma mia
				fread(&headerluts.front(), sizeof(HeaderLUT), SUBDOMAINS, file);
				
				{
					int c = fgetc(file);
					
					do 
					{ 
						printf("shouldnt be here! 0x%x\n", c); 
						//abort();
						c = fgetc(file);
					}
					while (! feof(file) );
				}	
				//assert(feof(file));
				
				for(int s = 0, currblock = 0; s < SUBDOMAINS; ++s)
				{
					const int nglobalchunks = lutchunks.size();
					
					
					const int nchunks = headerluts[s].nchunks;
					const size_t myamount = headerluts[s].aggregate_bytes;
					const size_t lutstart = base + myamount - sizeof(size_t) * nchunks;
					
					fseek(file, lutstart, SEEK_SET);
					vector<size_t> mylut(nchunks);
					fread(&mylut.front(), sizeof(size_t), nchunks, file);
					
					for(int i=0; i< nchunks; ++i)
						assert(mylut[i] < myamount);
					
					for(int i=1; i< nchunks; ++i)
						assert(mylut[i-1] < mylut[i]);
					
					for(int i = 0; i < mylut.size(); ++i)
					{
						assert(mylut[i] < myamount);
						mylut[i] += base;
					}
					
					assert(myamount > 0);
					base += myamount;
					assert(base <= global_header_displacement);
					
					//compute the base for this blocks
					for(int i = 0; i < BPS; ++i, ++currblock)
						metablocks[currblock].idcompression += nglobalchunks;
					
					lutchunks.insert(lutchunks.end(), mylut.begin(), mylut.end());					
				} 
				
				assert(base == global_header_displacement);
				
				lutchunks.push_back(base);
			}
			
			fclose(file);
		}
		
		idx2chunk.resize(NBLOCKS);
		
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
			
			idx2chunk[_id(entry.ix, entry.iy, entry.iz)] = compressedblock;
		}
		
		const bool verbose = true;
		
		if (verbose)
		{
			const size_t size_idx2chunk = idx2chunk.size() * sizeof(CompressedBlock);
			const double footprint_mb =  size_idx2chunk / 1024. / 1024.;
			printf("the header data is taking %.2f MB\n", footprint_mb);
		}
	}
	
	int xblocks() { return totalbpd[0]; } 
	int yblocks() { return totalbpd[1]; } 
	int zblocks() { return totalbpd[2]; } 
	
	void load_block(int ix, int iy, int iz, Real MYBLOCK[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_])
	{
		//printf("trying to load <%s>\n", path.c_str());
		FILE * f = fopen(path.c_str(), "rb");
		
		assert(f);
		
		CompressedBlock compressedchunk = idx2chunk[_id(ix, iy, iz)];
		
		size_t start = compressedchunk.start;
		
		assert(start >= miniheader_bytes);
		assert(start < global_header_displacement);
		assert(start + compressedchunk.extent <= global_header_displacement);
		
		vector<unsigned char> compressedbuf(compressedchunk.extent);
		fseek(f, compressedchunk.start, SEEK_SET);
		fread(&compressedbuf.front(), compressedchunk.extent, 1, f);
		
		assert(!feof(f));
		
		vector<unsigned char> waveletbuf(4 << 20);
		const size_t decompressedbytes = zdecompress(&compressedbuf.front(), compressedbuf.size(), &waveletbuf.front(), waveletbuf.size());
		
		int readbytes = 0;
		for(int i = 0; i<compressedchunk.subid; ++i)
		{
			int nbytes = * (int *) & waveletbuf[readbytes];
			readbytes += sizeof(int);
			readbytes += nbytes;
			
			assert(readbytes <= decompressedbytes);			
		}
		
		{
			int nbytes = *(int *)&waveletbuf[readbytes];
			readbytes += sizeof(int);
			assert(readbytes <= decompressedbytes);
			//printf("decompressing %d bytes...\n", nbytes);
			WaveletCompressor compressor;
			
			memcpy(compressor.data(), &waveletbuf[readbytes], nbytes);
			readbytes += nbytes;
			
			compressor.decompress(halffloat, nbytes, MYBLOCK);
		}
		
		fclose(f);
	}
};

class Reader_WaveletCompressionMPI: public Reader_WaveletCompression
{
	const MPI::Comm& comm;
	
public: 
	
	Reader_WaveletCompressionMPI(const MPI::Comm& comm, const string path): 
	Reader_WaveletCompression(path), comm(comm)
	{
		
	}
	
	virtual void load_file()
	{
		const int myrank = comm.Get_rank();
		
		if (myrank == 0)
			Reader_WaveletCompression::load_file();
		
		//propagate primitive type data members
		{			
			comm.Bcast(&global_header_displacement, sizeof(global_header_displacement), MPI_CHAR, 0);
			comm.Bcast(&miniheader_bytes, sizeof(miniheader_bytes), MPI_CHAR, 0);
			comm.Bcast(&NBLOCKS, sizeof(NBLOCKS), MPI_CHAR, 0);
			comm.Bcast(totalbpd, sizeof(totalbpd), MPI_CHAR, 0);
			comm.Bcast(bpd, sizeof(bpd), MPI_CHAR, 0);
			comm.Bcast(&halffloat, sizeof(halffloat), MPI_CHAR, 0);
		}
		
		size_t nentries = idx2chunk.size();
		
		comm.Bcast(&nentries, sizeof(nentries), MPI_CHAR, 0);
		
		if (myrank)		
			idx2chunk.resize(nentries);
		
		const size_t nbytes = nentries * sizeof(CompressedBlock);
		char * const entries = (char *)&idx2chunk.front();
		
		comm.Bcast(entries, nbytes, MPI_CHAR, 0);
	}
};
