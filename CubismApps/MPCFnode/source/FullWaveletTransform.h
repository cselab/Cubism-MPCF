/*
 *  FullWaveletTransform.h
 *  
 *
 *  Created by Diego Rossinelli on 3/27/13.
 *  Copyright 2013 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <vector>
#include <algorithm>
#include <bitset>

#include "WaveletsOnInterval.h"

using namespace std;

namespace WaveletsOnInterval 
{
	typedef WI4<false> ChosenWavelets;
	
	template<typename W>  
	inline const char * _name() { abort();  return "None"; }
	
	template<>  
	inline const char * _name<WI4<false> >() { return "InterpWavelet4thOrder"; }
	
	template<>  
	inline const char * _name<WI4<true> >() { return "LiftedInterpWavelet4thOrder"; }

	inline const char * ChosenWavelets_GetName() { return _name<ChosenWavelets>(); }
	
	template<int BS, bool lifting>
	struct FullTransform : WaveletSweep< ChosenWavelets >
	{
		FwtAp data[BS][BS][BS];
		
		FullTransform<BS/2, lifting> child;
		
		inline void fwt()
		{
			this->template sweep3D<BS, true>(data);
			
			static const int BSH = BS/2;
			for(int iz = 0; iz < BSH; ++iz)
				for(int iy = 0; iy < BSH; ++iy)
					for(int ix = 0; ix < BSH; ++ix)
						child.data[iz][iy][ix] = data[iz][iy][ix];
			
			child.fwt();
		}
		
		inline void iwt()
		{
			child.iwt();
			
			static const int BSH = BS/2;
			for(int iz = 0; iz < BSH; ++iz)
				for(int iy = 0; iy < BSH; ++iy)
					for(int ix = 0; ix < BSH; ++ix)
						data[iz][iy][ix] = child.data[iz][iy][ix];
			
			this->template sweep3D<BS, false>(data);
		}
		
		template<typename DataType>
		pair<vector<DataType>, bitset<BS * BS * BS> > threshold(const FwtAp eps)
		{
			static const int BSH = BS / 2;
			
			pair<vector<DataType>, bitset<BS * BS * BS> > retval;
			
			//code 0
			{
				pair<vector<DataType>, bitset<BSH * BSH * BSH> > childretval = child.template threshold<DataType>(eps);
				
				retval.first = childretval.first;
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
						{
							const int dst = ix + BS * (iy + BS * iz);
							const int src = ix + BSH * (iy + BSH * iz);
							
							retval.second[dst] = childretval.second[src];
						}
			}
			
			for(int code = 1; code < 8; ++code)
			{
				const int xstart = BSH * (code & 1);
				const int ystart = BSH * (code / 2 & 1);
				const int zstart = BSH * (code / 4 & 1);
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
						{
							const int xsrc = xstart + ix;
							const int ysrc = ystart + iy;
							const int zsrc = zstart + iz;
							
							const int dst = xsrc + BS * (ysrc + BS * zsrc);
							
							const DataType mydata = (DataType)data[zsrc][ysrc][xsrc];
							
							const bool accepted = fabs(mydata) > eps; 
							
							retval.second[dst] = accepted;
							
							if (accepted)
								retval.first.push_back(mydata);
							else 
								data[zsrc][ysrc][xsrc] = 0;
						}
			}
			
			return retval;
		}
		
		template<typename DataType>
		void load(vector<DataType>& datastream, bitset<BS * BS * BS> mask)
		{			
			static const int BSH = BS / 2;
			
			for(int code = 7; code >= 1; --code)
			{
				const int xstart = BSH * (code & 1);
				const int ystart = BSH * (code / 2 & 1);
				const int zstart = BSH * (code / 4 & 1);
				
				for(int iz = BSH - 1; iz >= 0; --iz)
					for(int iy = BSH - 1; iy >= 0; --iy)
						for(int ix = BSH - 1; ix >= 0; --ix)
						{
							const int myx = xstart + ix;
							const int myy = ystart + iy;
							const int myz = zstart + iz;
							
							const int srcidx = myx + BS * (myy + BS * myz);
							
							assert(srcidx >= 0);
							assert(srcidx < BS * BS * BS);
							
							const bool eat = mask[srcidx];
							
							const DataType mydata = eat ? datastream.back() : 0;
							
							data[myz][myy][myx] = mydata;
							
							if (eat) 
								datastream.pop_back();						
						}
			}
			
			//code 0
			{
				bitset<BSH * BSH * BSH> childmask;
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
						{
							const int src = ix + BS * (iy + BS * iz);
							const int dst = ix + BSH * (iy + BSH * iz);
							
							childmask[dst] = mask[src];
						}
				
				child.load(datastream, childmask);
			}
		}
	};
	
	template<bool lifting>
	struct FullTransform<4, lifting>
	{
		enum { BS = 4 } ;
		
		FwtAp data[BS][BS][BS];
		
		void fwt() { }
		
		void iwt() { }
		
		template<typename DataType>
		pair<vector<DataType>, bitset<BS * BS * BS> > threshold(const FwtAp eps) 
		{
			enum { N = BS * BS * BS };
			
			const FwtAp * const e = &data[0][0][0];
			
			vector<DataType> v(N);
			
			for(int i = 0; i < N; ++i)
				v[i] = e[i];
			
			pair<vector<DataType>, bitset<BS * BS * BS> > retval;
			retval.first = v;
			retval.second.set();
			
			return retval; 
		}
		
		template<typename DataType>
		void load(vector<DataType>& datastream, bitset<BS * BS * BS> mask)
		{
			assert(datastream.size() == BS * BS * BS);
			
			for(int iz = 0, c = 0; iz < BS; ++iz)
				for(int iy = 0; iy < BS; ++iy)
					for(int ix = 0; ix < BS; ++ix, ++c)
						data[iz][iy][ix] = datastream[c];
			
			datastream.clear();
		}
	};
}