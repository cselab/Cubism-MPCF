/*
 *  Types.h
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/19/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <cassert>

#include "common.h"

struct StateVector { Real r, u, v, w, s, levelset;};

//union GP { StateVector s; StateVector dsdt; };
struct GP {StateVector s; StateVector dsdt;};
template < int _LX, int _LY, int _LZ> 
struct MatrixGP
{
	static const int LX = _LX;
	static const int LY = _LY;
	static const int LZ = _LZ;
	
	//char screwup_alignment;
	GP __attribute__((aligned(16))) data[_LZ][_LY][_LX];
	
	
	inline GP& operator()(const int ix, const int iy, const int iz)
	{
		assert(ix >= 0); assert(ix < _LX);
		assert(iy >= 0); assert(iy < _LY);
		assert(iy >= 0); assert(iy < _LY);
		
		return data[iz][iy][ix];
	}
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(MatrixGP)/1024.;
	}
};

typedef MatrixGP<_BLOCKSIZE_, _BLOCKSIZE_, _BLOCKSIZE_> Block;

struct TestLab : MatrixGP<_BLOCKSIZE_+6, _BLOCKSIZE_+6, _BLOCKSIZE_+6>
{
	inline GP& operator()(const int ix, const int iy, const int iz)
	{
		assert(ix >= -3); assert(ix < _BLOCKSIZE_+3);
		assert(iy >= -3); assert(iy < _BLOCKSIZE_+3);
		assert(iz >= -3); assert(iz < _BLOCKSIZE_+3);
		
		return data[iz+3][iy+3][ix+3];
	}
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(TestLab)/1024.;
	}
};

struct TestLab_S2 : MatrixGP<_BLOCKSIZE_+2, _BLOCKSIZE_+2, _BLOCKSIZE_+2>
{
	inline GP& operator()(const int ix, const int iy, const int iz)
	{
		assert(ix >= -1); assert(ix < _BLOCKSIZE_+1);
		assert(iy >= -1); assert(iy < _BLOCKSIZE_+1);
		assert(iz >= -1); assert(iz < _BLOCKSIZE_+1);
		
		return data[iz+1][iy+1][ix+1];
	}
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(TestLab_S2)/1024.;
	}
};
