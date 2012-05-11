/*
 *  PUPKernelsMPI.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 10/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once
#include <vector>
#include <cassert>

void pack(const Real * const srcbase, Real * const dst, 
			   const unsigned int gptfloats,
			   int * selected_components, const int ncomponents,
			   const int xstart, const int ystart, const int zstart,
			   const int xend, const int yend, const int zend)
{	
	for(int idst=0, iz=zstart; iz<zend; ++iz)
		for(int iy=ystart; iy<yend; ++iy)
			for(int ix=xstart; ix<xend; ++ix)
			{
				const Real * src = srcbase + gptfloats*(ix + _BLOCKSIZEX_*(iy + _BLOCKSIZEY_*iz));
				
				for(int ic=0; ic<ncomponents; ic++, idst++)
					dst[idst] = src[selected_components[ic]];
			}
}

void pack_stripes(const Real * const srcbase, Real * const dst, 
					   const unsigned int gptfloats, 
					   const int selstart, const int selend, 
					   const int xstart, const int ystart, const int zstart,
					   const int xend, const int yend, const int zend)
{	
	for(int idst=0, iz=zstart; iz<zend; ++iz)
		for(int iy=ystart; iy<yend; ++iy)
			for(int ix=xstart; ix<xend; ++ix)
			{
				const Real * src = srcbase + gptfloats*(ix + _BLOCKSIZEX_*(iy + _BLOCKSIZEY_*iz));
				
				for(int ic=selstart; ic<selend; ic++, idst++)
					dst[idst] = src[ic];
			}
}

void unpack(const Real * const pack, Real * const dstbase, 
		  const unsigned int gptfloats,
		  const int * const selected_components, const int ncomponents,
		  const int nsrc,
		  const int dstxstart, const int dstystart, const int dstzstart,
		  const int dstxend, const int dstyend, const int dstzend,
		  const int xsize, const int ysize, const int zsize)
{
	for(int s=0, zd=dstzstart; zd<dstzend; ++zd)
		for(int yd=dstystart; yd<dstyend; ++yd)
			for(int xd=dstxstart; xd<dstxend; ++xd)
			{
				Real * const dst = dstbase + gptfloats * (xd + xsize * (yd + ysize * zd));
				for(int c=0; c<ncomponents; ++c, ++s)
					dst[selected_components[c]] = pack[s];
			}
}

void unpack_subregion(const Real * const pack, Real * const dstbase, 
			const unsigned int gptfloats,
			const int * const selected_components, const int ncomponents,
			const int srcxstart, const int srcystart, const int srczstart,
			const int LX, const int LY,
			const int dstxstart, const int dstystart, const int dstzstart,
			const int dstxend, const int dstyend, const int dstzend,
			const int xsize, const int ysize, const int zsize)
{
	for(int zd=dstzstart; zd<dstzend; ++zd)
		for(int yd=dstystart; yd<dstyend; ++yd)
			for(int xd=dstxstart; xd<dstxend; ++xd)
			{
				Real * const dst = dstbase + gptfloats * (xd + xsize * (yd + ysize * zd));
				const Real * src = pack + ncomponents*(xd - dstxstart + srcxstart + LX * (yd - dstystart + srcystart + LY * (zd - dstzstart + srczstart)));
				
				for(int c=0; c<ncomponents; ++c)
					dst[selected_components[c]] = src[c];
			}
}
