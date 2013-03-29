/*
 *  WaveletTexture3D.h
 *  
 *
 *  Created by Diego Rossinelli on 3/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <string>
#include <vector>
#include <sstream>

using namespace std;

#include <mpi.h>

#include "../../MPCFnode/source/WaveletCompressor.h"

typedef WaveletCompressorGeneric_zlib<_VOXELS_, float> TextureCompressor;

struct WaveletTexture3D
{
	struct Geometry
	{
		float pos[3], size[3]; 
		float texcoordstart[3],texcoordend[3];
		
		template<int dim>
		void setup(const int gstart, const int gend, const int ghost1side, const double gridspacing)
		{			
			pos[dim] = gstart * gridspacing;
			size[dim] = (gend - gstart) * gridspacing;
			
			const double voxelsize = 1. / _VOXELS_;

			texcoordstart[dim] = ghost1side * voxelsize;
			texcoordend[dim] = (gend - gstart - ghost1side) * voxelsize; 
		}
		
	} geometry;
		
	float data[_VOXELS_][_VOXELS_][_VOXELS_];
	
	TextureCompressor wavcomp;

	//not anymore: void clear() { memset(data, 0, sizeof(data)); }
	//not yet: void decompress(vector<char>& data);
	
	void compress(const float threshold, const bool halffloat, const unsigned char *& compresseddata, size_t& nbytes)
	{
		nbytes = wavcomp.compress(threshold, halffloat, data);
		compresseddata = (unsigned char *) wavcomp.data();
	}	
};

class WaveletTexture3D_MPICollection
{
	struct CompressedTexData { size_t start; int nbytes; WaveletTexture3D::Geometry geometry; }; 
	
	string path, header;
	size_t lutfile_start, bufferfile_start, * file_offset;

	const int xtextures, ytextures, ztextures, ntextures;
	const bool halffloat;
	const float wavelet_threshold;
	
	MPI::Win rmawindow;
	MPI::File myfile;
	const MPI::Intracomm& mycomm;

	string _header()
	{
		std::stringstream ss;
		
		ss << "\n==============START-ASCI-HEADER==============\n";
		
		{
			int one = 1;
			bool isone = *(char *)(&one);
			
			ss << "Endianess: " << (isone ? "little" : "big") << "\n";
		}
		
		ss << "Voxelsperdimension: " << _VOXELS_ << "\n";
		ss << "Textures: " << xtextures << " x "  << ytextures << " x " << ztextures  << "\n";
		ss << "HalfFloat: " << (this->halffloat ? "yes" : "no") << "\n";
		ss << "Wavelets: " << WaveletsOnInterval::ChosenWavelets_GetName() << "\n";
		ss << "Threshold: " << wavelet_threshold << "\n";
		ss << "Encoder: " << "zlib" << "\n";
		
		ss << "==============START-BINARY-LUT==============\n";
	
		//printf("my header is: <%s>\n", ss.str().c_str());
		return ss.str();
	}
	
public:
	
	WaveletTexture3D_MPICollection(const MPI::Intracomm& comm, 
								   const string path, const int xtextures, const int ytextures, const int ztextures,
								   const float wavelet_threshold, const bool halffloat): 
	mycomm(comm), path(path), xtextures(xtextures), ytextures(ytextures), ztextures(ztextures), ntextures(xtextures * ytextures * ztextures),
	wavelet_threshold(wavelet_threshold), halffloat(halffloat), file_offset(NULL)
	{		
		const int mygid = comm.Get_rank();

		//file setup
		{
			header = _header();
			
			string data_title = "==============START-BINARY-DATA==============\n";
			
			lutfile_start = header.size();
			bufferfile_start = lutfile_start + ntextures * sizeof(CompressedTexData);
			
			myfile = MPI::File::Open(mycomm, path.c_str(), MPI::MODE_CREATE | MPI::MODE_WRONLY, MPI::INFO_NULL);
						
			if (mygid == 0)
			{
				myfile.Write_at(0, header.c_str(), header.size(), MPI::CHAR);
				myfile.Write_at(bufferfile_start, data_title.c_str(), data_title.size(), MPI::CHAR);				
			}
			
			bufferfile_start += data_title.size();
		}
		
		//one sided communication setup
		{
			file_offset = (size_t *) MPI::Alloc_mem(sizeof(size_t), MPI::INFO_NULL);
			*file_offset = bufferfile_start; 
			assert(*file_offset != 0);
			
			if (mygid == 0)
				printf("at the beginning my offset was %d\n", *file_offset);
			
			rmawindow = MPI::Win::Create(file_offset, sizeof(size_t), sizeof(size_t), MPI::INFO_NULL, mycomm);
		}
	}
		
	void write(const int ix, const int iy, const int iz, WaveletTexture3D& texture)
	{
		//check that we are not totally nuts
		assert(ix >= 0 && ix < xtextures);
		assert(iy >= 0 && iy < ytextures);
		assert(iz >= 0 && iz < ztextures);

		//compress the texture
		const unsigned char * ptr = NULL;
		size_t nbytes = 0;
		texture.compress(wavelet_threshold, halffloat, ptr, nbytes);

		//spit some output for now
		{
			const size_t uncompressedbytes = sizeof(Real) * _VOXELS_ * _VOXELS_ * _VOXELS_;
			
			printf("Texture %d %d %d. Compression ratio: %.2f (threshold:%.2e)\n", 
				   ix, iy, iz, uncompressedbytes * 1. / nbytes, wavelet_threshold);	
		}

		//obtain file offset
		size_t myoffset = 0;
		{
			rmawindow.Lock(MPI::LOCK_EXCLUSIVE, 0, 0);
			rmawindow.Get(&myoffset, 1, MPI_UINT64_T, 0, 0, 1, MPI_UINT64_T); 
			rmawindow.Accumulate(&nbytes, 1, MPI_UINT64_T, 0, 0, 1, MPI_UINT64_T, MPI::SUM);
			rmawindow.Unlock(0);
		}
		assert(myoffset != 0);
		
		//write lut
		CompressedTexData entry = { myoffset, nbytes, texture.geometry };
		const size_t mylutoffset = sizeof(entry) * (ix + xtextures * (iy + ytextures * iz));
		myfile.Write_at(lutfile_start + mylutoffset, &entry, sizeof(entry), MPI::CHAR);

		//write data
		assert(ptr != NULL);
		assert(nbytes != 0);
		myfile.Write_at(myoffset, ptr, nbytes, MPI::CHAR);
	}
	
	~WaveletTexture3D_MPICollection() 
	{		
		rmawindow.Free();

		if (mycomm.Get_rank() == 0)
			printf("Terminating...closing file (%.2f kB) now.\n", *file_offset/1024.);
		
		MPI::Free_mem(file_offset);

		myfile.Close(); 		
	}
};
