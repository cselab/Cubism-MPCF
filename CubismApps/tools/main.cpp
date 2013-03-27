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
	
	Reader_WaveletCompression asd("../MPCFcluster/makefiles/datawavelet00000.StreamerGridPointIterative.channel0");

	printf("I found in total %dx%dx%d blocks.\n", asd.xblocks(), asd.yblocks(), asd.zblocks());
	
	asd.load_block(3, 0, 1, targetdata);
	
	return 0;
}