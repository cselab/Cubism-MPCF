/*
 *  WaveletSerializationTypes.h
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

struct BlockMetadata { int idcompression, subid, ix, iy, iz; };
struct HeaderLUT { size_t aggregate_bytes; int nchunks; };
struct CompressedBlock{ size_t start, extent; int subid; } ;
