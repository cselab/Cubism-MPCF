/*
 *  main.cpp
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#include "Reader_WaveletCompression.h"

int main()
{
	Real targetdata[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
	
	Reader_WaveletCompression asd("../../MPCFcluster/makefiles/datawavelet00000.StreamerGridPointIterative.channel0");

	printf("I found in total %dx%dx%d blocks.\n", asd.xblocks(), asd.yblocks(), asd.zblocks());
	
	asd.load_block(3, 0, 1, targetdata);
	
	if (false)
	{
		printf("OK FINAL TEST: THE DATA\n");
		for(int iz = 0; iz< _BLOCKSIZE_; ++iz)
			for(int iy = 0; iy< _BLOCKSIZE_; ++iy)
				for(int ix = 0; ix< _BLOCKSIZE_; ++ix)
					printf("%d %d %d: %e\n", ix, iy, iz, targetdata[iz][iy][ix]);
	}
	
	return 0;
}