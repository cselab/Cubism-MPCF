/*
 *  Test_Cloud.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 2/24/13.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "omp.h"

#include "Test_ShockBubble.h"
#include "Types.h"
#include <fstream>

namespace CloudData
{
    extern double seed_s[3], seed_e[3];
    extern Real min_rad, max_rad;
    extern int n_shapes, n_small, small_count;
};

//base class is a sphere
class shape
{
    Real center[3], radius;
    Real bbox_s[3], bbox_e[3];
    
public:
    shape()
    {
        const Real x_c = (Real)drand48()*(CloudData::seed_e[0]-CloudData::seed_s[0])+CloudData::seed_s[0];
        const Real y_c = (Real)drand48()*(CloudData::seed_e[1]-CloudData::seed_s[1])+CloudData::seed_s[1];
        
        const Real thickness = CloudData::seed_e[2]-CloudData::seed_s[2];        
        const Real z_c = (Real)drand48() * thickness + CloudData::seed_s[2];
        const Real _center[3] = {x_c, y_c, z_c};
        
        const Real _radius = (Real)drand48()*(CloudData::max_rad-CloudData::min_rad)+CloudData::min_rad;
        
        radius = _radius;
        
        for(int i=0; i<3; ++i)
        {
            center[i] = _center[i];
            bbox_s[i] = _center[i]-_radius-1.5*Simulation_Environment::EPSILON;
            bbox_e[i] = _center[i]+_radius+1.5*Simulation_Environment::EPSILON;
        }
    }
    
    shape(const Real _center[3], const Real _radius): radius(_radius)
    {
        for(int i=0; i<3; ++i)
        {
            center[i] = _center[i];
            bbox_s[i] = _center[i]-_radius-1.5*Simulation_Environment::EPSILON;
            bbox_e[i] = _center[i]+_radius+1.5*Simulation_Environment::EPSILON;
        }
    }
    
    void get_bbox(double s[3], double e[3]) const
    {
        for(int i=0; i<3; ++i)
        {
            s[i] = bbox_s[i];
            e[i] = bbox_e[i];
        }
    }
    
    bool inside_my_box(const Real pos[3]) const
    {
        const bool bXin = pos[0]>bbox_s[0] && pos[0]<bbox_e[0];
        const bool bYin = pos[1]>bbox_s[1] && pos[1]<bbox_e[1];
        const bool bZin = pos[2]>bbox_s[2] && pos[2]<bbox_e[2];
        
        return bXin && bYin && bZin;
    }
    
    void get_center(Real c[3]) const
    {
        for(int i=0; i<3; ++i)
            c[i] = center[i];
    }
    
    Real get_rad() const
    {
        return radius;
    }
    
    Real eval(const Real pos[3]) const
    {
        return sqrt(pow(pos[0]-center[0],2)+pow(pos[1]-center[1],2)+pow(pos[2]-center[2],2))-radius;
    }
    
    //Every other derived shape should implement this method.
    bool rejection_check(shape * this_shape, const Real start[3], const Real end[3]) const
    {
        double s[3], e[3];
        this->get_bbox(s,e);
        
        //this rule checks that the buble is inside the bounding box
        const bool bOut = s[0]<start[0] || s[1]<start[1] || s[2]<start[2] ||
        e[0]>end[0] || e[1]>end[1] || e[2]>end[2];
                
        if (bOut)
            return true;
        
        if(this!=this_shape)
        {
            double this_s[3], this_e[3];
            
            this_shape->get_bbox(this_s,this_e);
            
            const Real overlap_start[3] =
            {
                max(this_s[0], s[0]),
                max(this_s[1], s[1]),
                max(this_s[2], s[2])
            };
            
            const Real overlap_end[3] =
            {
                min(this_e[0], e[0]),
                min(this_e[1], e[1]),
                min(this_e[2], e[2])
            };
            
            const bool bOverlap = overlap_end[0] > overlap_start[0] && overlap_end[1] > overlap_start[1] && overlap_end[2] > overlap_start[2];
            
            if (bOverlap)
                return true;
        }
        
        return false;
    }
    
    Real * get_pointer_to_beginninig() {return &center[0];}
};

class Seed
{
    double start[3], end[3];
    
    vector<shape> v_shapes;
    
public:

    Seed(const double _start[3], const double _end[3])
    {
        start[0] = _start[0];
        start[1] = _start[1];
        start[2] = _start[2];
        
        end[0] = _end[0];
        end[1] = _end[1];
        end[2] = _end[2];
    }
        
    void make_shapes(int ntocheck)
    {        			
		ifstream f_read_cloud("cloud.dat");
		
		if (!f_read_cloud.good())
		{
			cout << "Watchout! cant read the file <cloud.dat>. Aborting now...\n";
			abort();
		}
		
		while(true)
		{
			if (!f_read_cloud.good())
				abort();
			
			int idx;
			Real c[3], rad;
			f_read_cloud >> idx >> c[0] >> c[1] >> c[2] >> rad;			
			
			cout << "shape " << idx << " " <<  c[0] << " " << c[1] << " " << c[2] << " " << rad << endl;
			
			shape cur_shape(c,rad);
			v_shapes.push_back(cur_shape);
			
			if (f_read_cloud.eof()) break;
		}
		
		f_read_cloud.close();
	
		cout << "number of shapes are " << v_shapes.size() << endl;
		
		if (v_shapes.size() != ntocheck)
		{
			 cout << "PROBLEM! ntocheck is " << ntocheck << " which does not correspond to the number of shapes!!!\n";
		}       
    }
    
    Seed retain_shapes(const double mystart[3], const double extent[3]) const
    {
		assert(v_shapes.size() > 0);
		
		const double myend[3] = {
			mystart[0] + extent[0],
			mystart[1] + extent[1],
			mystart[2] + extent[2]
		};
		
		Seed retval(mystart, myend);
		
		 for(int i=0; i<v_shapes.size(); ++i)
         {
				
                   shape curr_shape = v_shapes[i];
        
                    double s[3],e[3];
                    curr_shape.get_bbox(s,e);
                    
                    const Real xrange =  min(myend[0], e[0]) - max(mystart[0], s[0]);
                    const Real yrange =  min(myend[1], e[1]) - max(mystart[1], s[1]);
                    const Real zrange =  min(myend[2], e[2]) - max(mystart[2], s[2]);
                    const bool bOverlap = (xrange > 0) && (yrange > 0) && (zrange > 0);

					if (bOverlap) retval.v_shapes.push_back(curr_shape);
        }
        
        return retval;
    }
    
    vector<shape> get_shapes() const { return v_shapes; }
    
    int get_shapes_size() const { return v_shapes.size(); }
    
    void set_shapes(vector<shape> v) { v_shapes = v; }
};

inline double eval(const vector<shape>& v_shapes, const Real pos[3])
{
    Real d = HUGE_VAL;
    
    for( int i=0; i<v_shapes.size(); ++i)
    {
		const Real newdistance = v_shapes[i].eval(pos);
		
        d = min(d, newdistance);
	}
	
    const double alpha = M_PI*min(1., max(0., (double)(d + 4 * Simulation_Environment::EPSILON)/(4*Simulation_Environment::EPSILON)));
    
    return 0.5 + 0.5 * cos(alpha);
}

class Test_Cloud: public Test_ShockBubble
{
    friend class Test_CloudMPI;
    
    bool bRestartedSeed;
    
    void _ic(FluidGrid & grid);
    
    void _my_ic(FluidGrid& grid, const Seed myseed);
    void _my_ic_quad(FluidGrid& grid, const Seed myseed);
    void _set_energy(FluidGrid& grid);
    void _initialize_cloud();
    
public:
    Test_Cloud(const int argc, const char ** argv);
    
    void setup();
};

//ALERT: Energy channel forces pressure for solving
//the Laplace p = 0 at initial condition
template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLabCloudLaplace: public BlockLab<BlockType,allocator>
{
	typedef typename BlockType::ElementType ElementTypeBlock;
    
protected:
	bool is_xperiodic() {return false;}
	bool is_yperiodic() {return false;}
	bool is_zperiodic() {return false;}
    
public:
	BlockLabCloudLaplace(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
        BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
        
        ElementTypeBlock b;
        b.clear();
        b.rho = 1000;
        b.u = b.v = b.w = 0;
        b.G = 1./(6.59-1);
        b.P = 4.049e3*b.G*6.59;
        b.energy = 100.0;
        
        if (info.index[0]==0)           bc.template applyBC_dirichlet<0,0>(b);
        if (info.index[0]==this->NX-1)  bc.template applyBC_dirichlet<0,1>(b);
        if (info.index[1]==0)			bc.template applyBC_dirichlet<1,0>(b);
        if (info.index[1]==this->NY-1)	bc.template applyBC_dirichlet<1,1>(b);
        if (info.index[2]==0)			bc.template applyBC_dirichlet<2,0>(b);
        if (info.index[2]==this->NZ-1)	bc.template applyBC_dirichlet<2,1>(b);
    }
};
