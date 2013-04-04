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
    extern Real seed_s[3], seed_e[3], min_rad, max_rad;
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
        const Real x_c = (Real)drand48();
        const Real y_c = (Real)drand48();
        const Real z_c = (Real)drand48();
        const Real _center[3] = {x_c, y_c, z_c};
        
        const Real _radius = (Real)drand48()*(CloudData::max_rad-CloudData::min_rad)+CloudData::min_rad;
        
        radius = _radius;
        
        for(int i=0; i<3; ++i)
        {
            center[i] = _center[i];
            bbox_s[i] = _center[i]-_radius-2.0*Simulation_Environment::EPSILON;
            bbox_e[i] = _center[i]+_radius+2.0*Simulation_Environment::EPSILON;
        }
    }
    
    shape(const Real _center[3], const Real _radius): radius(_radius)
    {
        for(int i=0; i<3; ++i)
        {
            center[i] = _center[i];
            bbox_s[i] = _center[i]-_radius-2.0*Simulation_Environment::EPSILON;
            bbox_e[i] = _center[i]+_radius+2.0*Simulation_Environment::EPSILON;
        }
    }
    
    void get_bbox(Real s[3], Real e[3]) const
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
        //this rule is to stop making smaller bubbles after n_small smalle bubbles.
        const Real rad_av_rough = 0.5*(CloudData::min_rad+CloudData::max_rad);
        const Real rad_cut = 1.2*rad_av_rough;
        //if (CloudData::small_count==CloudData::n_small && this_shape->get_rad()<rad_cut)
        //   return true;
        
        Real s[3], e[3];
        this->get_bbox(s,e);
        
        //this rule checks that the buble is inside the bounding box
        const bool bOut = s[0]<start[0] || s[1]<start[1] || s[2]<start[2] ||
        e[0]>end[0] || e[1]>end[1] || e[2]>end[2];
        
        //this rule is to keep larger bubbles outside
        //larger w.r.t. rad_avg_rough
        const Real c_box[3] = {end[0]+start[0], end[1]+start[1], end[2]+start[2]};
        Real this_c[3];
        this_shape->get_center(this_c);
        const Real d_c_box = pow(this_c[0]-0.5*c_box[0],2)+pow(this_c[1]-0.5*c_box[1],2)+pow(this_c[2]-0.5*c_box[2],2);
        const Real R_cutoff = pow( pow((end[0]-start[0])*(end[1]-start[1])*(end[2]-start[2]),(Real)3)*0.75/(8*M_PI),1./3);
        const bool bLargeOuter = (d_c_box>=R_cutoff && this_shape->get_rad()>=rad_av_rough);
        
        if (bOut || bLargeOuter)
            return true;
        
        if(this!=this_shape)
        {
            Real this_s[3], this_e[3];
            
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
        
        if (this_shape->get_rad()<rad_cut)
            CloudData::small_count++;
        
        return false;
    }
};

class Seed
{
    Real start[3], end[3];
    int n_shapes;
    vector<shape *> v_shapes;
    
public:
    Seed(const Real _start[3], const Real _end[3], const int _n_shapes): n_shapes(_n_shapes)
    {
        start[0] = _start[0];
        start[1] = _start[1];
        start[2] = _start[2];
        
        end[0] = _end[0];
        end[1] = _end[1];
        end[2] = _end[2];
    }
    
    void make_shapes(const bool bRestartedSeed=false)
    {
        bool bFull = false;
        
        unsigned int restarted_seed;
        
        //read seed if needed
        if (bRestartedSeed)
        {
            //For the moment we just read the cloud.dat
            //and ignore seed and cloud_config.dat
            {
                ifstream f_read_cloud("cloud.dat");
                for(int i=0; i<CloudData::n_shapes; ++i)
                {
                    int idx;
                    Real c[3], rad;
                    f_read_cloud >> idx >> c[0] >> c[1] >> c[2] >> rad;
                    cout << "shape " << idx << " " <<  c[0] << " " << c[1] << " " << c[2] << " " << rad << endl;
                    
                    shape * cur_shape = new shape(c,rad);
                    v_shapes.push_back(cur_shape);
                }
                
                f_read_cloud.close();
                
                cout << "number of shapes are " << v_shapes.size() << endl;
            }
            
            return;
            
            ifstream f_read("seed.dat");
            if(f_read)
            {
                cout << "seed file is there" << endl;
                f_read >> restarted_seed;
                f_read.close();
            }
            else
            {
                cout << "seed file not there...aborting" << endl;
                abort();
            }
        }
        
        const unsigned int seed = bRestartedSeed? restarted_seed : time(0);
        
        //save seed anyway
        {
            ofstream f_save("seed.dat");
            f_save<<seed;
            f_save.close();
        }
        
        srand48(seed);
        
        while(!bFull)
        {
            shape * cur_shape = new shape();
            
            if (reject_check(cur_shape))
            {
                delete cur_shape;
                continue;
            }
            
            v_shapes.push_back(cur_shape);
            
            printf("size is %ld out of %d\n", v_shapes.size(), n_shapes);
            
            bFull = v_shapes.size()==n_shapes;
        }
        
        //We are even more paranoic so we save the bubble centers and radii
        //If the first method does not work on different platforms,
        //a method must be implemented to read the following data.
        {
            ofstream f_save("cloud.dat");
            for(int i=0; i< v_shapes.size(); i++)
            {
                Real c[3], rad;
                rad = v_shapes[i]->get_rad();
                v_shapes[i]->get_center(c);
                f_save<<i<< " " << c[0] << " " << c[1] << " " << c[2] << " " << rad << endl;
            }
            
            f_save.close();
        }
    }
    
    void make_shapes(const Real mystart[3], const Real myend[3], const bool bRestartedSeed=false, const bool isroot=false)
    {
        bool bFull = false;
        
        unsigned int restarted_seed;
        
        //read seed if needed
        if (bRestartedSeed)
        {
            v_shapes.clear();
            //For the moment we just read the cloud.dat
            //and ignore seed and cloud_config.dat
            {
                ifstream f_read_cloud("cloud.dat");
                for(int i=0; i<CloudData::n_shapes; ++i)
                {
                    int idx;
                    Real c[3], rad;
                    f_read_cloud >> idx >> c[0] >> c[1] >> c[2] >> rad;
                    if(isroot)
                        cout << "shape " << idx << " " <<  c[0] << " " << c[1] << " " << c[2] << " " << rad << endl;
                    
                    shape * cur_shape = new shape(c,rad);
                    
                    Real s[3],e[3];
                    cur_shape->get_bbox(s,e);
                    
                    const Real overlap_start[3] =
                    {
                        max(mystart[0]-CloudData::max_rad, s[0]),
                        max(mystart[1]-CloudData::max_rad, s[1]),
                        max(mystart[2]-CloudData::max_rad, s[2])
                    };
                    
                    const Real overlap_end[3] =
                    {
                        min(myend[0]+CloudData::max_rad, e[0]),
                        min(myend[1]+CloudData::max_rad, e[1]),
                        min(myend[2]+CloudData::max_rad, e[2])
                    };
                    
                    const bool bOverlap = overlap_end[0] > overlap_start[0] && overlap_end[1] > overlap_start[1] && overlap_end[2] > overlap_start[2];
                    
                    if (bOverlap)
                        v_shapes.push_back(cur_shape);
                    else
                        delete cur_shape;
                }
                
                f_read_cloud.close();
                
                if(isroot) cout << "number of shapes are " << v_shapes.size() << endl;
            }
            
            return;
            
            ifstream f_read("seed.dat");
            if(f_read)
            {
                if(isroot) cout << "seed file is there" << endl;
                f_read >> restarted_seed;
                f_read.close();
            }
            else
            {
                if(isroot) cout << "seed file not there...aborting" << endl;
                abort();
            }
        }
        
        //save seed anyway
        if(isroot)
        {
            const unsigned int seed = bRestartedSeed? restarted_seed : time(0);
            
            ofstream f_save("seed.dat");
            f_save<<seed;
            f_save.close();
            
            srand48(seed);
            
            while(!bFull)
            {
                shape * cur_shape = new shape();
                
                if (reject_check(cur_shape))
                {
                    delete cur_shape;
                    continue;
                }
                
                v_shapes.push_back(cur_shape);
                
                printf("size is %ld out of %d\n", v_shapes.size(), n_shapes);
                
                bFull = v_shapes.size()==n_shapes;
            }
            
            //We are even more paranoic so we save the bubble centers and radii
            //If the first method does not work on different platforms,
            //a method must be implemented to read the following data.
            {
                ofstream f_save("cloud.dat");
                for(int i=0; i< v_shapes.size(); i++)
                {
                    Real c[3], rad;
                    rad = v_shapes[i]->get_rad();
                    v_shapes[i]->get_center(c);
                    f_save<<i<< " " << c[0] << " " << c[1] << " " << c[2] << " " << rad << endl;
                }
                
                f_save.close();
            }
        }

        //ok let them all find their own subset
        this->make_shapes(mystart, myend, true, isroot);
    }
    
    bool reject_check(shape * cur_shape) const
    {
        if (v_shapes.size()==0)
            return cur_shape->rejection_check(cur_shape, start, end);
        
        for(int i=0; i<v_shapes.size(); ++i)
        {
            shape * this_shape = v_shapes[i];
            if (cur_shape->rejection_check(this_shape, start, end))
                return true;
        }
        
        return false;
    }
    
    vector<shape*> get_vshapes() const
    {
        return v_shapes;
    }
};

inline double eval(const vector<shape*> v_shapes, const Real pos[3])
{
    Real d = HUGE_VAL;
    
    for( int i=0; i<v_shapes.size(); ++i)
        d = min(d, v_shapes[i]->eval(pos));
    
    const double alpha = M_PI*min(1., max(0., (d+0.5*Simulation_Environment::EPSILON)/Simulation_Environment::EPSILON));
    return 0.5+0.5*cos(alpha);
}

class Test_Cloud: public Test_ShockBubble
{
    friend class Test_CloudMPI;
    
    bool bRestartedSeed;
    
    void _ic(FluidGrid & grid);
    
    void _my_ic(FluidGrid& grid, const vector< shape * > v_shapes);
    void _my_ic_quad(FluidGrid& grid, const vector< shape * > v_shapes);
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
