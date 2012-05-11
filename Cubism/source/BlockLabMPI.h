/*
 *  BlockLabMPI.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 10/21/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once
#include "GridMPI.h"
#include "SynchronizerMPI.h"

template<typename MyBlockLab>
class BlockLabMPI : public MyBlockLab
{
	const SynchronizerMPI * refSynchronizerMPI;
	typedef typename MyBlockLab::BlockType BlockType;
	
protected:
	int mypeindex[3], pesize[3], mybpd[3];
	int gLastX, gLastY, gLastZ;
	
public:
	template< typename TGrid >
	void prepare(GridMPI<TGrid>& grid, const SynchronizerMPI& SynchronizerMPI)
	{
		refSynchronizerMPI = &SynchronizerMPI;
		refSynchronizerMPI->getpedata(mypeindex, pesize, mybpd);
		StencilInfo stencil = refSynchronizerMPI->getstencil();
		assert(stencil.isvalid());
		MyBlockLab::prepare(grid, stencil.sx,  stencil.ex,  stencil.sy,  stencil.ey,  stencil.sz,  stencil.ez, stencil.tensorial);
		gLastX = grid.getBlocksPerDimension(0)-1;
		gLastY = grid.getBlocksPerDimension(1)-1;
		gLastZ = grid.getBlocksPerDimension(2)-1;
	}
	
	void load(const BlockInfo& info, const Real t=0, const bool applybc=true)
	{		
		MyBlockLab::load(info, t, false);
		
		assert(refSynchronizerMPI != NULL);
		
		const int xorigin = mypeindex[0]*mybpd[0];		
		const int yorigin = mypeindex[1]*mybpd[1];
		const int zorigin = mypeindex[2]*mybpd[2];
		
		const bool xskin = (info.index[0] == xorigin || info.index[0] == xorigin + mybpd[0]-1);
		const bool yskin = (info.index[1] == yorigin || info.index[1] == yorigin + mybpd[1]-1);
		const bool zskin = (info.index[2] == zorigin || info.index[2] == zorigin + mybpd[2]-1);
		
		const bool xboundary = info.index[0]==0 || info.index[0]==gLastX;
		const bool yboundary = info.index[1]==0 || info.index[1]==gLastY;
		const bool zboundary = info.index[2]==0 || info.index[2]==gLastZ;		
		
		const bool any_periodic = this->is_xperiodic() || this->is_yperiodic() || this->is_zperiodic();
		const bool any_boundary = xboundary || yboundary || zboundary;
		
		if ((xskin || yskin || zskin))// && (!any_boundary || any_boundary && any_periodic))
		{
			const bool xperiodic = true;//this->is_xperiodic();
			const bool yperiodic = true;//this->is_yperiodic();
			const bool zperiodic = true;//this->is_zperiodic();
			
			const int rsx = !xperiodic && info.index[0]==0? 0 : this->m_stencilStart[0];
			const int rex = !xperiodic && info.index[0]==gLastX ? BlockType::sizeX : (BlockType::sizeX+this->m_stencilEnd[0]-1);
			const int rsy = !yperiodic && info.index[1]==0? 0 : this->m_stencilStart[1];
			const int rey = !yperiodic && info.index[1]==gLastY ? BlockType::sizeY : (BlockType::sizeY+this->m_stencilEnd[1]-1);
			const int rsz = !zperiodic && info.index[2]==0? 0 : this->m_stencilStart[2];
			const int rez = !zperiodic && info.index[2]==gLastZ ? BlockType::sizeZ : (BlockType::sizeZ+this->m_stencilEnd[2]-1);
			
			//printf("%d %d %d %d %d %d\n", rsx, rex, rsy, rey, rsz, rez);
			
			Real * const dst = (Real *)&this->m_cacheBlock->LinAccess(0);
			
			typedef typename MyBlockLab::ElementType ET;
			
			refSynchronizerMPI->fetch((const Real*)info.ptrBlock, dst, 
									  this->m_stencilStart[0], this->m_stencilStart[1], this->m_stencilStart[2],
									  this->m_cacheBlock->getSize()[0], this->m_cacheBlock->getSize()[1], this->m_cacheBlock->getSize()[2],
									  sizeof(ET)/sizeof(Real),
									  rsx, rex, rsy, rey, rsz, rez);
		}
		
		if (applybc) MyBlockLab::_apply_bc(info, t);
	}
};

//..and this ends the story.
/*
template<typename Lab, typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabMPI_Absorbing: public BlockLabMPI<Lab>
{
    typedef typename BlockType::ElementType ElementTypeBlock;
	
public:
	BlockLabMPI_Absorbing(): BlockLabMPI<Lab>(){}
	
	void load(const BlockInfo& info, const Real t=0)
	{
		BlockLabMPI<Lab>::load(info);
		
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
        const int xorigin = this->mypeindex[0]*this->mybpd[0];		
		const int yorigin = this->mypeindex[1]*this->mybpd[1];
		const int zorigin = this->mypeindex[2]*this->mybpd[2];
		
		const bool x0skin = info.index[0] == xorigin && this->mypeindex[0]==0;
        const bool x1skin = info.index[0] == xorigin + this->mybpd[0]-1 && this->mypeindex[0]==this->pesize[0]-1;
		const bool y0skin = info.index[1] == yorigin && this->mypeindex[1]==0;
        const bool y1skin = info.index[1] == yorigin + this->mybpd[1]-1 && this->mypeindex[1]==this->pesize[1]-1;
		const bool z0skin = info.index[2] == zorigin && this->mypeindex[2]==0;
        const bool z1skin = info.index[2] == zorigin + this->mybpd[2]-1 && this->mypeindex[2]==this->pesize[2]-1;
        
		if (x0skin)    bc.template applyBC_absorbing<0,0>();
		if (x1skin)    bc.template applyBC_absorbing<0,1>();
		if (y0skin)    bc.template applyBC_absorbing<1,0>();
		if (y1skin)    bc.template applyBC_absorbing<1,1>();
		if (z0skin)    bc.template applyBC_absorbing<2,0>();
		if (z1skin)    bc.template applyBC_absorbing<2,1>();
	}
};

template<typename Lab, typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabMPI_Reflecting: public BlockLabMPI<Lab>
{
    typedef typename BlockType::ElementType ElementTypeBlock;
    
public:
	BlockLabMPI_Reflecting(): BlockLabMPI<Lab>(){}
	
	void load(const BlockInfo& info, const Real t=0)
	{
		BlockLabMPI<Lab>::load(info);
		
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
        
        const int xorigin = this->mypeindex[0]*this->mybpd[0];		
		const int yorigin = this->mypeindex[1]*this->mybpd[1];
		const int zorigin = this->mypeindex[2]*this->mybpd[2];
		
		const bool x0skin = info.index[0] == xorigin && this->mypeindex[0]==0;
        const bool x1skin = info.index[0] == xorigin + this->mybpd[0]-1 && this->mypeindex[0]==this->pesize[0]-1;
		const bool y0skin = info.index[1] == yorigin && this->mypeindex[1]==0;
        const bool y1skin = info.index[1] == yorigin + this->mybpd[1]-1 && this->mypeindex[1]==this->pesize[1]-1;
		const bool z0skin = info.index[2] == zorigin && this->mypeindex[2]==0;
        const bool z1skin = info.index[2] == zorigin + this->mybpd[2]-1 && this->mypeindex[2]==this->pesize[2]-1;
        
        if (x0skin)    bc.template applyBC_reflecting<0,0>();
        if (x1skin)    bc.template applyBC_reflecting<0,1>();
        if (y0skin)    bc.template applyBC_reflecting<1,0>();
        if (y1skin)    bc.template applyBC_reflecting<1,1>();
        if (z0skin)    bc.template applyBC_reflecting<2,0>();
        if (z1skin)    bc.template applyBC_reflecting<2,1>();
	}
};

template<typename Lab, typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabMPI_CVT: public BlockLabMPI<Lab>
{
    typedef typename BlockType::ElementType ElementTypeBlock;
    
public:
	BlockLabMPI_CVT(): BlockLabMPI<Lab>(){}
	
	void load(const BlockInfo& info, const Real t=0)
	{
		BlockLabMPI<Lab>::load(info);
		
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
        
        const int xorigin = this->mypeindex[0]*this->mybpd[0];		
		const int yorigin = this->mypeindex[1]*this->mybpd[1];
		const int zorigin = this->mypeindex[2]*this->mybpd[2];
		
		const bool x0skin = info.index[0] == 0 && this->mypeindex[0]==0;
        const bool x1skin = info.index[0] == xorigin + this->mybpd[0]-1 && this->mypeindex[0]==this->pesize[0]-1;
		const bool z0skin = info.index[2] == 0 && this->mypeindex[2]==0;
        const bool z1skin = info.index[2] == zorigin + this->mybpd[2]-1 && this->mypeindex[2]==this->pesize[2]-1;
        
        if (x0skin)    bc.template applyBC_reflecting<0,0>();
        if (x1skin)    bc.template applyBC_reflecting<0,1>();
        if (z0skin)    bc.template applyBC_reflecting<2,0>();
        if (z1skin)    bc.template applyBC_reflecting<2,1>();
	}
};
*/
