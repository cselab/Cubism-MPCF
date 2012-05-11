/*
 *  BlockLab.h
 *  Cubism
 *
 *  Created by Diego Rossinelli on 5/24/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */
//
#pragma once

#include "Matrix3D.h"
#include "Grid.h"
//#include "Concepts.h"

/**
 * Working copy of Block + Ghosts.
 * Data of original block is copied (!) here. So when changing something in
 * the lab we are not changing the original data.
 * Requirements:
 * - Concepts::BlockLabElementType<ElementTypeT>
 * - Concepts::Castable<typename BlockType::ElementType, ElementTypeT>
 */
template<typename TBlock, template<typename X> class allocator = std::allocator, typename ElementTypeT = typename TBlock::ElementType>
class BlockLab
{
    // check concepts
	//   CONCEPT_CHECK(Concepts::BlockLabElementType<ElementTypeT>);
	//   CONCEPT_CHECK(Concepts::Castable<typename BlockType::ElementType, ElementTypeT>);
	
public:
	typedef ElementTypeT ElementType;
	
protected:
	typedef TBlock BlockType;
	typedef typename BlockType::ElementType ElementTypeBlock;
	enum eBlockLab_State {eMRAGBlockLab_Prepared, eMRAGBlockLab_Loaded, eMRAGBlockLab_Uninitialized};
	
	eBlockLab_State m_state;
	Matrix3D<ElementType, true, allocator> * m_cacheBlock;
	
	int m_stencilStart[3], m_stencilEnd[3];
	int NX, NY, NZ;
	
	bool istensorial;
	
	const Grid<BlockType, allocator>* m_refGrid;
	
	virtual void _apply_bc(const BlockInfo& info, const Real t=0) { }
	
	template<typename T>
	void _release(T *& t) 
	{ 
		if (t != NULL)
		{
			allocator<T>().destroy(t);
			allocator<T>().deallocate(t,1); 
		}
		t = NULL; 
	}
	
public:
	
	BlockLab():
	m_state(eMRAGBlockLab_Uninitialized),
	m_cacheBlock(NULL),
	m_refGrid(NULL)
	{
		m_stencilStart[0] = m_stencilStart[1] = m_stencilStart[2] = 0;
		m_stencilEnd[0] = m_stencilEnd[1] = m_stencilEnd[2] = 0;
	}
	
	virtual bool is_xperiodic() { return true; }
	virtual bool is_yperiodic() { return true; }
	virtual bool is_zperiodic() { return true; }
	
	~BlockLab()
	{
		_release(m_cacheBlock);
	}
	
	template <int dim>
	int getActualSize() const
	{
		assert(dim>=0 && dim<3);
		return m_cacheBlock->getSize()[dim];
	}
	
	inline ElementType * getBuffer() const { return &m_cacheBlock->LinAccess(0);}
	
	void prepare(Grid<BlockType,allocator>& grid, int startX, int endX, int startY, int endY, int startZ, int endZ, const bool istensorial)
	{
		const int ss[3] = {startX, startY, startZ};
		const int se[3] = {endX, endY, endZ};
		prepare(grid, ss, se, istensorial);
	}
	
	/**
	 * Prepare the extended block.
	 * @param collection    Collection of blocks in the grid (e.g. result of Grid::getBlockCollection()).
	 * @param boundaryInfo  Info on the boundaries of the grid (e.g. result of Grid::getBoundaryInfo()).
	 * @param stencil_start Maximal stencil used for computations at lower boundary.
	 *                      Defines how many ghosts we will get in extended block.
	 * @param stencil_end   Maximal stencil used for computations at lower boundary.
	 *                      Defines how many ghosts we will get in extended block.
	 */
	
	void prepare(Grid<BlockType,allocator>& grid, const int stencil_start[3], const int stencil_end[3], const bool istensorial)
	{
		NX = grid.getBlocksPerDimension(0);
		NY = grid.getBlocksPerDimension(1);
		NZ = grid.getBlocksPerDimension(2);
		
		this->istensorial = istensorial;
		
		m_refGrid = &grid;
		
		assert(stencil_start[0]>= -BlockType::sizeX);
		assert(stencil_start[1]>= -BlockType::sizeY);
		assert(stencil_start[2]>= -BlockType::sizeZ);
		assert(stencil_end[0] < BlockType::sizeX*2);
		assert(stencil_end[1] < BlockType::sizeY*2);
		assert(stencil_end[2] < BlockType::sizeZ*2);
		
		m_stencilStart[0] = stencil_start[0];
		m_stencilStart[1] = stencil_start[1];
		m_stencilStart[2] = stencil_start[2];
		
		m_stencilEnd[0] = stencil_end[0];
		m_stencilEnd[1] = stencil_end[1];
		m_stencilEnd[2] = stencil_end[2];
		
		assert(m_stencilStart[0]<=m_stencilEnd[0]);
		assert(m_stencilStart[1]<=m_stencilEnd[1]);
		assert(m_stencilStart[2]<=m_stencilEnd[2]);
		
		if (m_cacheBlock == NULL || 
			m_cacheBlock->getSize()[0]!= BlockType::sizeX + m_stencilEnd[0] - m_stencilStart[0] -1 ||
			m_cacheBlock->getSize()[1]!= BlockType::sizeY + m_stencilEnd[1] - m_stencilStart[1] -1 ||
			m_cacheBlock->getSize()[2]!= BlockType::sizeZ + m_stencilEnd[2] - m_stencilStart[2] -1 )
		{
			if (m_cacheBlock != NULL)
				_release(m_cacheBlock);
			
			m_cacheBlock = allocator< Matrix3D<ElementType,  true, allocator> >().allocate(1);
			
			allocator< Matrix3D<ElementType,  true, allocator> >().construct(m_cacheBlock, Matrix3D<ElementType,  true, allocator> ());
			
			m_cacheBlock->_Setup(BlockType::sizeX + m_stencilEnd[0] - m_stencilStart[0] -1,
								 BlockType::sizeY + m_stencilEnd[1] - m_stencilStart[1] -1,
								 BlockType::sizeZ + m_stencilEnd[2] - m_stencilStart[2] -1);
			
		}
		
		m_state = eMRAGBlockLab_Prepared;
		//	printf("ss: %d %d %d  se: %d %d %d\n", m_stencilStart[0], m_stencilStart[1], m_stencilStart[2], 
		//			m_stencilEnd[0], m_stencilEnd[1], m_stencilEnd[2]);
	}
	
	/**
	 * Load a block (incl. ghosts for it).
	 * This is not called internally but by the BlockProcessing-class. Hence a new version of BlockLab,
	 * can just overwrite it and through template-passing to BlockProcessing, the right version will be
	 * called.
	 * @param info  Reference to info of block to be loaded.
	 */
	void load(const BlockInfo& info, const Real t=0, const bool applybc=true)
	{
		const Grid<BlockType,allocator>& grid = *m_refGrid;
		
		//0. couple of checks
		//1. load the block into the cache
		//2. put the ghosts into the cache
		
		//0.
		assert(m_state == eMRAGBlockLab_Prepared || m_state==eMRAGBlockLab_Loaded);
		assert(m_cacheBlock != NULL);
		
		const int nX = BlockType::sizeX;
		const int nY = BlockType::sizeY;
		const int nZ = BlockType::sizeZ;
		
		//1.
		{
			assert(sizeof(ElementType) == sizeof(typename BlockType::ElementType));
			
			BlockType& block = *(BlockType *)info.ptrBlock;
			
			ElementTypeBlock * ptrSource = &block(0);
			
			for(int iz=0; iz<nZ; iz++)
				for(int iy=0; iy<nY; iy++)
				{					
					ElementType * ptrDestination = &m_cacheBlock->Access(0-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]);
					
					//for(int ix=0; ix<nX; ix++, ptrSource++, ptrDestination++)
					//	*ptrDestination = (ElementType)*ptrSource;
					memcpy(ptrDestination, ptrSource, sizeof(ElementType)*nX);
					
					ptrSource+= nX;
				}
		}
		
		//2.
		{
			const bool xperiodic = is_xperiodic();
			const bool yperiodic = is_yperiodic();
			const bool zperiodic = is_zperiodic();
			
			const bool xskin = info.index[0]==0 || info.index[0]==grid.getBlocksPerDimension(0)-1;
			const bool yskin = info.index[1]==0 || info.index[1]==grid.getBlocksPerDimension(1)-1;
			const bool zskin = info.index[2]==0 || info.index[2]==grid.getBlocksPerDimension(2)-1;
			
			const int xskip = info.index[0]==0 ? -1 : 1;
			const int yskip = info.index[1]==0 ? -1 : 1;
			const int zskip = info.index[2]==0 ? -1 : 1;
			
			for(int icode=0; icode<27; icode++)
			{
				if (icode == 1*1 + 3*1 + 9*1) continue;
				
				const int code[3] = { icode%3-1, (icode/3)%3-1, (icode/9)%3-1};
				
				if (!xperiodic && code[0] == xskip && xskin) continue;
				if (!yperiodic && code[1] == yskip && yskin) continue;
				if (!zperiodic && code[2] == zskip && zskin) continue;
				
				if (!istensorial && abs(code[0])+abs(code[1])+abs(code[2])>1) continue;
				
				const int s[3] = { 
					code[0]<1? (code[0]<0 ? m_stencilStart[0]:0 ) : nX, 
					code[1]<1? (code[1]<0 ? m_stencilStart[1]:0 ) : nY, 
					code[2]<1? (code[2]<0 ? m_stencilStart[2]:0 ) : nZ };
				
				const int e[3] = {
					code[0]<1? (code[0]<0 ? 0:nX ) : nX+m_stencilEnd[0]-1, 
					code[1]<1? (code[1]<0 ? 0:nY ) : nY+m_stencilEnd[1]-1, 
					code[2]<1? (code[2]<0 ? 0:nZ ) : nZ+m_stencilEnd[2]-1};
				
				if (!grid.avail(info.index[0] + code[0], info.index[1] + code[1], info.index[2] + code[2])) continue;
				
				BlockType& b = grid(info.index[0] + code[0], info.index[1] + code[1], info.index[2] + code[2]);
				
				for(int iz=s[2]; iz<e[2]; iz++)
					for(int iy=s[1]; iy<e[1]; iy++)
					{
						/*	char * ptrDest = (char*)&m_cacheBlock->Access(0-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]);
						 const char * ptrSrc = (const char*)&b(0 - code[0]*BlockType::sizeX, iy - code[1]*BlockType::sizeY, iz - code[2]*BlockType::sizeZ);
						 const int bytes = (e[0]-s[0])*sizeof(ElementType);
						 memcpy(ptrDest, ptrSrc, bytes);*/
						for(int ix=s[0]; ix<e[0]; ix++)
							m_cacheBlock->Access(ix-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]) = 
							(ElementType)b(ix - code[0]*BlockType::sizeX, iy - code[1]*BlockType::sizeY, iz - code[2]*BlockType::sizeZ);
					}
			}
			
			if (applybc) _apply_bc(info, t);			
			
			m_state = eMRAGBlockLab_Loaded;
		}
	}
	
	/**
	 * Get a single element from the block.
	 * stencil_start and stencil_end refer to the values passed in BlockLab::prepare().
	 *
	 * @param ix    Index in x-direction (stencil_start[0] <= ix < BlockType::sizeX + stencil_end[0] - 1).
	 * @param iy    Index in y-direction (stencil_start[1] <= iy < BlockType::sizeY + stencil_end[1] - 1).
	 * @param iz    Index in z-direction (stencil_start[2] <= iz < BlockType::sizeZ + stencil_end[2] - 1).
	 */
	ElementType& operator()(int ix, int iy=0, int iz=0)
	{
		assert(m_state == eMRAGBlockLab_Loaded);
		
		const int nX = m_cacheBlock->getSize()[0];
		const int nY = m_cacheBlock->getSize()[1];
		const int nZ = m_cacheBlock->getSize()[2];
		
		assert(ix-m_stencilStart[0]>=0 && ix-m_stencilStart[0]<nX);
		assert(iy-m_stencilStart[1]>=0 && iy-m_stencilStart[1]<nY);
		assert(iz-m_stencilStart[2]>=0 && iz-m_stencilStart[2]<nZ);
		
		return m_cacheBlock->Access(ix-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]);
	}
	
	/** Just as BlockLab::operator() but returning a const. */
	const ElementType& read(int ix, int iy=0, int iz=0) const
	{
		assert(m_state == eMRAGBlockLab_Loaded);
		
		const int nX = m_cacheBlock->getSize()[0];
		const int nY = m_cacheBlock->getSize()[1];
		const int nZ = m_cacheBlock->getSize()[2];
		
		assert(ix-m_stencilStart[0]>=0 && ix-m_stencilStart[0]<nX);
		assert(iy-m_stencilStart[1]>=0 && iy-m_stencilStart[1]<nY);
		assert(iz-m_stencilStart[2]>=0 && iz-m_stencilStart[2]<nZ);
		
		return m_cacheBlock->Access(ix-m_stencilStart[0], iy-m_stencilStart[1], iz-m_stencilStart[2]);
	}
	
private:
	
	//forbidden
	BlockLab(const BlockLab&):
	m_state(eMRAGBlockLab_Uninitialized),
	m_cacheBlock(NULL){abort();}
	
	BlockLab& operator=(const BlockLab&){abort(); return *this;}
};


template<typename BlockType, template<typename X> class allocator = std::allocator, typename ElementTypeT = typename BlockType::ElementType >
class EmptyBlockLab:  BlockLab<BlockType, allocator, ElementTypeT>
{
	int current_size[3];
	
public:	
	
	EmptyBlockLab():BlockLab<BlockType, allocator, ElementTypeT>()
	{
		current_size[0] = current_size[1] = current_size[2] = 0;
	}
	
	typedef ElementTypeT ElementType;
	inline void prepare(Grid<BlockType,allocator>& grid, const int stencil_start[3], const int stencil_end[3])
	{
		this->m_refGrid = &grid;
		
		assert(stencil_start[0]>= -BlockType::sizeX);
		assert(stencil_start[1]>= -BlockType::sizeY);
		assert(stencil_start[2]>= -BlockType::sizeZ);
		assert(stencil_end[0] < BlockType::sizeX*2);
		assert(stencil_end[1] < BlockType::sizeY*2);
		assert(stencil_end[2] < BlockType::sizeZ*2);
		
		this->m_stencilStart[0] = stencil_start[0];
		this->m_stencilStart[1] = stencil_start[1];
		this->m_stencilStart[2] = stencil_start[2];
		
		this->m_stencilEnd[0] = stencil_end[0];
		this->m_stencilEnd[1] = stencil_end[1];
		this->m_stencilEnd[2] = stencil_end[2];
		
		assert(this->m_stencilStart[0]<=this->m_stencilEnd[0]);
		assert(this->m_stencilStart[1]<=this->m_stencilEnd[1]);
		assert(this->m_stencilStart[2]<=this->m_stencilEnd[2]);
		
		if (current_size[0]!= BlockType::sizeX + this->m_stencilEnd[0] - this->m_stencilStart[0] -1 ||
			current_size[1]!= BlockType::sizeY + this->m_stencilEnd[1] - this->m_stencilStart[1] -1 ||
			current_size[2]!= BlockType::sizeZ + this->m_stencilEnd[2] - this->m_stencilStart[2] -1 )
		{
			current_size[0] = BlockType::sizeX + this->m_stencilEnd[0] - this->m_stencilStart[0] -1;
			current_size[1] = BlockType::sizeY + this->m_stencilEnd[1] - this->m_stencilStart[1] -1;
			current_size[2] = BlockType::sizeZ + this->m_stencilEnd[2] - this->m_stencilStart[2] -1;
		}
		
		this->m_state = BlockLab<BlockType, allocator, ElementTypeT>::eMRAGBlockLab_Prepared;		
	}
};

