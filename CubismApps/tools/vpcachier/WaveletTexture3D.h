/*
 *  WaveletTexture3D.h
 *  
 *
 *  Created by Diego Rossinelli on 3/28/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <vector>

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

	void clear() { memset(data, 0, sizeof(data)); }
		
	vector<unsigned char> compress();
	
	void decompress(vector<char>& data);
};


