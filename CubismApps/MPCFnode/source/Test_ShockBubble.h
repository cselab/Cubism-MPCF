/*
 *  Test_ShockBubble.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 6/16/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Test_SteadyState.h"

class Test_ShockBubble: public Test_SteadyState
{
    friend class Test_ShockBubbleMPI;
    
	Real radius, bubble_pos[3];
    
	void _ic(FluidGrid& grid);
	void _setup_constants();
    void _dumpStatistics(FluidGrid& grid, const int counter, const Real t, const Real dt);
    
public:	
	Test_ShockBubble(const int argc, const char ** argv): Test_SteadyState(argc, argv) { }
    
	void run();
	void setup();
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabBubble: public BlockLab<BlockType,allocator>
{		
	typedef typename BlockType::ElementType ElementTypeBlock;

public:
	BlockLabBubble(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{	
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
            
		if (info.index[0]==0)           bc.template applyBC_absorbing_better_faces<0,0>();		
		if (info.index[0]==this->NX-1)  bc.template applyBC_absorbing_better_faces<0,1>();
        if (info.index[1]==0)			bc.template applyBC_absorbing_better_faces<1,0>();
		if (info.index[1]==this->NY-1)	bc.template applyBC_absorbing_better_faces<1,1>();
		if (info.index[2]==0)			bc.template applyBC_absorbing_better_faces<2,0>();
		if (info.index[2]==this->NZ-1)	bc.template applyBC_absorbing_better_faces<2,1>();
        
        if (info.index[0]==0 || info.index[0]==this->NX-1 || info.index[1]==0 || info.index[1]==this->NY-1 || info.index[2]==0 || info.index[2]==this->NZ-1) 
            bc.applyBC_absorbing_better_tensorials();
	}
};
