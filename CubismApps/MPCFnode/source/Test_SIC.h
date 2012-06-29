/*
 *  Test_SIC.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 5/30/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Test_ShockBubble.h"

class Test_SIC: public Test_ShockBubble
{
    friend class Test_SICMPI;
    
	void _ic(FluidGrid& grid);
    
public:	
	Test_SIC(const int argc, const char ** argv): Test_ShockBubble(argc, argv) { }

	  void setup();
};


template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabCollapse: public BlockLab<BlockType,allocator>
{		
	typedef typename BlockType::ElementType ElementTypeBlock;
    
protected:
	bool is_xperiodic() {return false;}
	bool is_yperiodic() {return false;}
	bool is_zperiodic() {return false;}
    
public:
	BlockLabCollapse(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{	
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
        
        ElementTypeBlock b;
        b.clear();
        b.rho = 10;
        b.G = 1./(6.59-1);
#ifdef _LIQUID_
        b.P = 4.049e4*b.G*6.59;
        b.energy = 1e4*b.G + b.P;
#else
        b.energy = 1e4*b.G;
#endif
        
        //if (info.index[0]==0)           bc.template applyBC_spaceDirichlet<0,0>(b,t,info.h_gridpoint);//pressure pulse		
        if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(b);
        if (info.index[0]==this->NX-1)  bc.template applyBC_reflecting<0,1>(); //Wall
        if (info.index[1]==0)			bc.template applyBC_absorbing_better_faces<1,0>();
        if (info.index[1]==this->NY-1)	bc.template applyBC_absorbing_better_faces<1,1>();
        if (info.index[2]==0)			bc.template applyBC_absorbing_better_faces<2,0>();
        if (info.index[2]==this->NZ-1)	bc.template applyBC_absorbing_better_faces<2,1>();
        
        const bool bEdgeXY = (info.index[0]==0 || info.index[0]==this->NX-1) && (info.index[1]==0 || info.index[1]==this->NY-1);
        const bool bEdgeYZ = (info.index[1]==0 || info.index[1]==this->NY-1) && (info.index[2]==0 || info.index[2]==this->NZ-1);
        const bool bEdgeZX = (info.index[2]==0 || info.index[2]==this->NZ-1) && (info.index[0]==0 || info.index[0]==this->NX-1);
        
        const bool bCorner = (info.index[0]==0 || info.index[0]==this->NX-1) && (info.index[1]==0 || info.index[1]==this->NY-1) && (info.index[2]==0 || info.index[2]==this->NZ-1);
        
        if (this->istensorial)
        {
            if (bEdgeXY || bEdgeYZ || bEdgeZX && !bCorner) 
                bc.applyBC_absorbing_better_tensorials_edges();
            
            if (bCorner)
                bc.applyBC_absorbing_better_tensorials_corners();
        }
    }
};
