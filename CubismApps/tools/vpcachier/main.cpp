/*
 *  main.cpp
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>
#include <string>
#include <mpi.h>

#include <ArgumentParser.h>

#include "Reader_WaveletCompression.h"
#include "WaveletTexture3D.h"

int main(int argc, const char **  argv)
{
	MPI::Init_thread(MPI_THREAD_SERIALIZED);
	
	//create my cartesian communicator
	int pesize[3] = {0, 0, 0};
	
	MPI::Compute_dims(MPI::COMM_WORLD.Get_size(), 3, pesize);
	
	const bool isperiodic[3] = { false, false, false };	
	
	MPI::Cartcomm mycartcomm = MPI::COMM_WORLD.Create_cart(3, pesize, isperiodic, true);
	
	MPI::Intracomm& mycomm = mycartcomm;
	
	const int xpe = pesize[0];
	const int ype = pesize[1];
	const int zpe = pesize[2];
	
	const int myrank = mycomm.Get_rank();
	const bool isroot = !myrank;
	
	if (isroot)
		printf("PE are organized as:  %d %d %d\n", xpe, ype, zpe);
	
	//open the file
	//std::string path = "../../MPCFcluster/makefiles/datawavelet00000.StreamerGridPointIterative.channel0";
	std::string path = "datawavelet00000.StreamerGridPointIterative.channel0";
	
	Reader_WaveletCompressionMPI myreader(mycomm, path);
	
	myreader.load_file();
	
	//figure out the amount of work per rank
	const int xblocks = myreader.xblocks();
	const int yblocks = myreader.yblocks();
	const int zblocks = myreader.zblocks();
	
	const double gridspacing = 1. / max(max(xblocks, yblocks), zblocks);
	
	MYASSERT(_BLOCKSIZE_ <= _VOXELS_, "_BLOCKSIZE_ = " << _BLOCKSIZE_ << " is bigger than _VOXELS_ =" << _VOXELS_ << "\n");
	
	const int ghosts1side = 2;
	const int puredata1d = _VOXELS_ - 2 * ghosts1side;
	
	const int xtextures = (xblocks * _BLOCKSIZE_ - 2 * ghosts1side) / puredata1d;
	const int ytextures = (yblocks * _BLOCKSIZE_ - 2 * ghosts1side) / puredata1d;
	const int ztextures = (zblocks * _BLOCKSIZE_ - 2 * ghosts1side) / puredata1d;
	
	MYASSERT(xtextures > 0 && ytextures > 0 && ztextures > 0, "NUMBER OF VP BLOCKS IS ZERO!");
	
	const int xrankwork = std::max(1, xtextures / xpe);
	const int yrankwork = std::max(1, ytextures / ype);
	const int zrankwork = std::max(1, ztextures / zpe);
	
	int peindex[3] = {-1, -1, -1};
	mycartcomm.Get_coords(myrank, 3, peindex);
	
	assert(peindex[0] >= 0 && peindex[1] >= 0 && peindex[2] >= 0);
	
	const int myxstart = xrankwork * peindex[0];
	const int myystart = yrankwork * peindex[1];
	const int myzstart = zrankwork * peindex[2];
	
	const int myxend = std::min(xtextures, xrankwork * (peindex[0] + 1));
	const int myyend = std::min(ytextures, yrankwork * (peindex[1] + 1));
	const int myzend = std::min(ztextures, zrankwork * (peindex[2] + 1));
	const int mytotalwork = (myxend - myxstart) * (myyend - myystart) * (myzend - myzstart);
	
	//spit a warning in case we starve
	{
		int leastwork = -1;
		mycomm.Reduce(&mytotalwork, &leastwork, 1, MPI_INT, MPI_SUM, 0);
		
		assert(leastwork >= 0 || !isroot);
		
		if (isroot && leastwork == 0)
			printf("WARNING: watchout because some of the ranks have zero work.\n");
	}
	
	for(int gz = myzstart ;  gz < myzend; ++gz)
		for(int gy = myystart ;  gy < myyend; ++gy)
			for(int gx = myxstart ;  gx < myxend; ++gx)
			{
				WaveletTexture3D * texture = new WaveletTexture3D;
								
				const int xdatastart = gx * puredata1d;
				const int ydatastart = gy * puredata1d;
				const int zdatastart = gz * puredata1d;
				
				const int xdataend = (gx + 1) * puredata1d + 2 * ghosts1side;
				const int ydataend = (gy + 1) * puredata1d + 2 * ghosts1side;
				const int zdataend = (gz + 1) * puredata1d + 2 * ghosts1side;
				
				texture->geometry.setup<0>(xdatastart, xdataend, ghosts1side, gridspacing);
				texture->geometry.setup<1>(ydatastart, ydataend, ghosts1side, gridspacing);
				texture->geometry.setup<2>(zdatastart, zdataend, ghosts1side, gridspacing);
				
				assert(xdataend - xdatastart == _VOXELS_);
				assert(ydataend - ydatastart == _VOXELS_);
				assert(zdataend - zdatastart == _VOXELS_);
				
				const int xblockstart = xdatastart / _BLOCKSIZE_;
				const int yblockstart = ydatastart / _BLOCKSIZE_;
				const int zblockstart = zdatastart / _BLOCKSIZE_;
				
				const int xblockend = (xdataend - 1) / _BLOCKSIZE_ + 1;
				const int yblockend = (ydataend - 1) / _BLOCKSIZE_ + 1;
				const int zblockend = (zdataend - 1) / _BLOCKSIZE_ + 1;
				
				for(int bz = zblockstart; bz < zblockend; ++bz)
					for(int by = yblockstart; by < yblockend; ++by)
						for(int bx = xblockstart; bx < xblockend; ++bx)
						{
							Real data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
							
							assert(bx >= 0 && bx < xblocks);
							assert(by >= 0 && by < yblocks);
							assert(bz >= 0 && bz < zblocks);
							
							myreader.load_block(bx, by, bz, data);
							
							const int xdststart = std::max(xdatastart, bx * _BLOCKSIZE_) - xdatastart;
							const int ydststart = std::max(ydatastart, by * _BLOCKSIZE_) - ydatastart;
							const int zdststart = std::max(zdatastart, bz * _BLOCKSIZE_) - zdatastart;

							const int xdstend = std::min(xdataend, (bx + 1) * _BLOCKSIZE_) - xdatastart;
							const int ydstend = std::min(ydataend, (by + 1) * _BLOCKSIZE_) - ydatastart;
							const int zdstend = std::min(zdataend, (bz + 1) * _BLOCKSIZE_) - zdatastart;
							
							const int xsrcoffset = xdatastart - bx * _BLOCKSIZE_;
							const int ysrcoffset = ydatastart - by * _BLOCKSIZE_;
							const int zsrcoffset = zdatastart - bz * _BLOCKSIZE_;
							
							for(int dz = zdststart; dz < zdstend; ++dz)
								for(int dy = ydststart; dy < ydstend; ++dy)
									for(int dx = xdststart; dx < xdstend; ++dx)
									{
										assert(dx >= 0 && dx < _VOXELS_);
										assert(dy >= 0 && dy < _VOXELS_);
										assert(dz >= 0 && dz < _VOXELS_);
										
										assert(dx + xsrcoffset >= 0 && dx + xsrcoffset < _BLOCKSIZE_);
										assert(dy + ysrcoffset >= 0 && dy + ysrcoffset < _BLOCKSIZE_);
										assert(dz + zsrcoffset >= 0 && dz + zsrcoffset < _BLOCKSIZE_);
										
										texture->data[dz][dy][dx] = data[dz + zsrcoffset][dy + ysrcoffset][dx + xsrcoffset];
										//texture->data[dz][dy][dx] = 1;
									}
						}
				
				/*for(int dz = 0; dz < _VOXELS_; ++dz)
					for(int dy = 0; dy < _VOXELS_; ++dy)
						for(int dx = 0; dx < _VOXELS_; ++dx)
							assert(texture->data[dz][dy][dx] == 1);*/
				
				assert(xblockstart >= 0 && xblockstart < xblocks);
				assert(yblockstart >= 0 && yblockstart < yblocks);
				assert(zblockstart >= 0 && zblockstart < zblocks);
				
				assert(xblockend > 0 && xblockend <= xblocks);
				assert(yblockend > 0 && yblockend <= yblocks);
				assert(zblockend > 0 && zblockend <= zblocks);
				
				vector<unsigned char> datacompressed = texture->compress();
				
				delete texture;
			}
	
	printf("good up to here\n");
	
	MPI::Finalize();
	
	if (myrank == 0)
		std::cout << "Also tchuess zame.\n";
	
	return 0;
}
