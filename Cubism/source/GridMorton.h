/*
 *  GridMorton.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/24/11.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#include <iostream>

using namespace std;

#include "Grid.h"
#include "Indexers.h"

template <typename TGrid>
class GridMorton: public TGrid
{
	IndexerMorton indexer;
	vector<BlockInfo> cached_infos;

	vector<BlockInfo> _getBlocksInfo() const
	{
		std::vector<BlockInfo> r;
		r.reserve(this->N);

		const unsigned int nX = this->NX;
		const unsigned int nY = this->NY;
		const unsigned int nZ = this->NZ;

		const double h = (1./nX);

		unsigned int ix, iy, iz;
		for(int i=0;i<this->N; ++i)
		{
			indexer.decode(i, ix, iy, iz);

			if (ix>=nX || iy>=nY || iz>=nZ) continue;

			const int idx[3] = {ix, iy, iz};
			const double origin[3] = {ix*h, iy*h, iz*h};

			r.push_back(BlockInfo(i, idx, origin, h, h/TBlock::sizeX, this->_linaccess(i)));
		}

		return r;
	}

	bool _is_p2(int x)
	{
		return ( (x > 0) && ((x & (x - 1)) == 0) );
	}

	public:

	typedef typename TGrid::BlockType TBlock;

	GridMorton(unsigned int nX, unsigned int nY=1, unsigned int nZ=1):
		TGrid(nX, nY, nZ), indexer(nX,nY,nZ)
	{

		if (! _is_p2(nX) || !_is_p2(nY) || !_is_p2(nZ) || nX!=nY || nX!=nZ)
		{
			std::cout << "GridMorton::GridMorton: Ooops: GridMorton works only when BPDX=DPBY=DPBZ AND BPDX is a power of 2. Aborting now." << std::endl;
			abort();
		}

		cached_infos = _getBlocksInfo();

		std::cout << "____________________________________ YOU ARE USING MORTON GRIDS _______________________________" << std::endl;
	}

	TBlock& operator()(unsigned int ix, unsigned int iy=0, unsigned int iz=0) const
	{
		const unsigned int nX = this->NX;
		const unsigned int nY = this->NY;
		const unsigned int nZ = this->NZ;

		assert(ix>= - nX && ix<2*nX);
		assert(iy>= - nY && iy<2*nY);
		assert(iz>= - nZ && iz<2*nZ);

		const unsigned int idx = indexer.encode((ix+nX) % nX, (iy+nY) % nY, (iz+nZ) % nZ);

		return *(this->_linaccess(idx));
	}

	vector<BlockInfo> getBlocksInfo() const
	{
		return cached_infos;
	}
};
