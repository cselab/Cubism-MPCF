/*
 *  Reader_WaveletCompression.h
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <map>
#include <string>
#include <vector>

#include "../MPCFnode/source/WaveletCompressor.h"


#include "../MPCFcluster/source/WaveletSerializationTypes.h"
#include "../MPCFcluster/source/CompressionEncoders.h"
#include "../MPCFnode/source/FullWaveletTransform.h"

using namespace std;

class Reader_WaveletCompression
{
	string path;
	int miniheader_bytes;
	size_t global_header_displacement;
	
	int NBLOCKS;
	int totalbpd[3], bpd[3];
	
	vector<BlockMetadata> metablocks;
	
	map<int, map<int, map<int, CompressedBlock > > > meta2subchunk;
	vector<size_t> lutchunks;
	
	bool halffloat;
	
public:
	
	Reader_WaveletCompression(string path): NBLOCKS(-1), global_header_displacement(-1), path(path)
	{
		for(int i = 0; i < 3; ++i)
			totalbpd[i] = -1;
		
		for(int i = 0; i < 3; ++i)
			bpd[i] = -1;
		
		//THE FIRST PART IS SEQUENTIAL
		//THE SECOND ONE IS RANDOM ACCESS
		
		string binaryocean_title = "\n==============START-BINARY-OCEAN==============\n";	
		miniheader_bytes = sizeof(size_t) + binaryocean_title.size();		
		
		//random access data structures: meta2subchunk, lutchunks;
		
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
	}
	
	int xblocks() { return totalbpd[0]; } 
	int yblocks() { return totalbpd[1]; } 
	int zblocks() { return totalbpd[2]; } 
	
	void load_block(int ix, int iy, int iz, Real MYBLOCK[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_])
	{
		FILE * f = fopen(path.c_str(), "rb");
		assert(f);
		
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
	
};
