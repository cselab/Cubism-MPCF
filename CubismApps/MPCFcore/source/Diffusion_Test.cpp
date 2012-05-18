//
//  Diffusion_Test.cpp
//  MPCFcore
//
//  Created by Babak Hejazialhosseini on 9/16/11.
//  Copyright 2011 ETH Zurich. All rights reserved.
//

#include <cstring>
#include <algorithm>

using namespace std;

#include "Diffusion_Test.h"

//This code is NOT optimized and is a second-order accurate
//shear-stress tensor naive implementation.

typedef double RealTemp;

class Diffusion_CPP_M1
{	
	RealTemp dtinvh, nu1, nu2, g1, g2;
    
    RealTemp _project(const RealTemp phi, const RealTemp rho);
    
public:
    struct VectorElement
    {
        RealTemp a[6];
        
        VectorElement(const RealTemp _a, const RealTemp _b, const RealTemp _c, const RealTemp _d, const RealTemp _e, const RealTemp _f)
        {
            a[0]=_a;
            a[1]=_b;
            a[2]=_c;
            a[3]=_d;
            a[4]=_e;
            a[5]=_f;
        }
        
        RealTemp operator[](int i) const
        {
            return a[i];
        }
    };
    
    VectorElement _toPrimitive(StateVector & p);
    
	Diffusion_CPP_M1(const RealTemp _dtinvh, const RealTemp _nu1=0, const RealTemp _nu2=0, const RealTemp _g1=0, const RealTemp _g2=0): dtinvh(_dtinvh), nu1(_nu1), nu2(_nu2), g1(_g1), g2(_g2) { }
    
    void compute_conservative(TestLab_S2& lab, const BlockInfo& info, Block& o, const RealTemp a);
};

void Diffusion_Test::_gold(TestLab_S2& lab, Block& block, const Real _nu1, const Real _nu2, const Real _g1, const Real _g2, const Real a, const Real _dtinvh, const Real _h)
{	
	Diffusion_CPP_M1 gold(_dtinvh, _nu1, _nu2, _g1, _g2);
	
	const int idx[3] = {0,0,0};
	const double pos[3] = {0, 0, 0};
	BlockInfo info(0, idx, pos, 1.,  _h, NULL);
	gold.compute_conservative(lab, info, block, a);
}

void Diffusion_Test::_compare(Block& _a, Block& _b, double accuracy, string kernelname)
{
	double maxe[6] = {0,0,0,0,0,0};
	double sume[6] = {0,0,0,0,0,0};
	
	double relmaxe[6] = {0,0,0,0,0,0};
	double relsume[6] = {0,0,0,0,0,0};
	
	for(int iz = 0; iz<_BLOCKSIZE_; iz++)
		for(int iy = 0; iy<_BLOCKSIZE_; iy++)
			for(int ix = 0; ix<_BLOCKSIZE_; ix++)
			{
				StateVector a = _a(ix, iy, iz).dsdt;
				StateVector b = _b(ix, iy, iz).dsdt;
				const double s[6]  = {
					b.r ,
					b.u ,
					b.v ,
					b.w ,
					b.s ,
					b.levelset
				};
				
				const double e[6]  = {
					b.r - a.r,
					b.u - a.u,
					b.v - a.v,
					b.w - a.w,
					b.s - a.s,
					b.levelset - a.levelset
				};
				
				for(int i=0; i<6; i++)
					if (fabs(e[i])/fabs(s[i])>accuracy && fabs(e[i])>accuracy) printf("significant error at %d %d %d %d -> e=%e (rel is %e, values are %e %e)\n", ix, iy, iz, i, e[i], e[i]/s[i],s[i],s[i]-e[i]);
				
				for(int i=0; i<6; i++)
					maxe[i] = max(fabs(e[i]), maxe[i]);
				
				for(int i=0; i<6; i++)
					relmaxe[i] = max(fabs(e[i])/max(accuracy, fabs(s[i])), relmaxe[i]);
				
				for(int i=0; i<6; i++)
					sume[i] += fabs(e[i]);
				
				for(int i=0; i<6; i++)
					relsume[i] += fabs(e[i])/max(accuracy, fabs(s[i]));
			}
	
	printf("\tLinf discrepancy:\t");
	for(int i=0; i<6; i++)
		printf("%.2e ", maxe[i]);
	
	cout << endl;
	printf("\tL1 (dh=1):       \t");
	for(int i=0; i<6;i++)
		printf("%.2e ", sume[i]);
	cout<<endl;
	
	printf("\tLinf RELATIVE discrepancy:\t");
	for(int i=0; i<6; i++)
		printf("%.2e ", relmaxe[i]);
	
	cout << endl;
	printf("\tL1 RELATIVE (dh=1):       \t");
	for(int i=0; i<6;i++)
		printf("%.2e ", relsume[i]);
	cout<<endl;
}

//Assumed Input GridPoint for Diffusion (AIGPD)
struct AIGP_D { RealTemp r, u, v, w, s, l; }; 

typedef Diffusion_CPP_M1::VectorElement VectorElement;

VectorElement Diffusion_CPP_M1::_toPrimitive(StateVector & p)
{
    const RealTemp pressure = (p.s-0.5/p.r*(p.u*p.u+p.v*p.v+p.w*p.w))/p.levelset;
    
    VectorElement v(p.r, p.u/p.r, p.v/p.r, p.w/p.r, pressure, p.levelset);
    return v;
}

RealTemp Diffusion_CPP_M1::_project(const RealTemp phi, const RealTemp rho)
{
    const RealTemp hs = 1-min(max((phi-g2)/(g1-g2),(RealTemp)0),(RealTemp)1);
    const RealTemp nu = (nu1==nu2)? nu1 : nu2*hs + nu1*(1-hs);
    return rho*nu;//
}

void Diffusion_CPP_M1::compute_conservative(TestLab_S2& lab, const BlockInfo& info, Block& o, const RealTemp a)
{
    const int r = 0;
    const int u = 1;
    const int v = 2;
    const int w = 3;
    
    const RealTemp h = info.h_gridpoint;
    
    for(int iz=0; iz<_BLOCKSIZE_; iz++)
        for(int iy=0; iy<_BLOCKSIZE_; iy++)
            for(int ix=0; ix<_BLOCKSIZE_; ix++)
            {                
                RealTemp c000 = _toPrimitive(lab(ix,iy,iz).s)[u];
                RealTemp c100 = _toPrimitive(lab(ix+1,iy,iz).s)[u];
                RealTemp c110 = _toPrimitive(lab(ix+1,iy+1,iz).s)[u];
                RealTemp c010 = _toPrimitive(lab(ix,iy+1,iz).s)[u];
                RealTemp c101 = _toPrimitive(lab(ix+1,iy,iz+1).s)[u];
                RealTemp c001 = _toPrimitive(lab(ix,iy,iz+1).s)[u];
                RealTemp c111 = _toPrimitive(lab(ix+1,iy+1,iz+1).s)[u];
                RealTemp c011 = _toPrimitive(lab(ix,iy+1,iz+1).s)[u];
                RealTemp cm100 = _toPrimitive(lab(ix-1,iy,iz).s)[u];
                RealTemp cm110 = _toPrimitive(lab(ix-1,iy+1,iz).s)[u];
                RealTemp cm101 = _toPrimitive(lab(ix-1,iy,iz+1).s)[u];
                RealTemp cm111 = _toPrimitive(lab(ix-1,iy+1,iz+1).s)[u];
                RealTemp c0m10 = _toPrimitive(lab(ix,iy-1,iz).s)[u];
                RealTemp cm1m10 = _toPrimitive(lab(ix-1,iy-1,iz).s)[u];
                RealTemp c0m11 = _toPrimitive(lab(ix,iy-1,iz+1).s)[u];
                RealTemp cm1m11 = _toPrimitive(lab(ix-1,iy-1,iz+1).s)[u];
                RealTemp c1m10 = _toPrimitive(lab(ix+1,iy-1,iz).s)[u];
                RealTemp c1m11 = _toPrimitive(lab(ix+1,iy-1,iz+1).s)[u];
                RealTemp c10m1 = _toPrimitive(lab(ix+1,iy,iz-1).s)[u];
                RealTemp c00m1 = _toPrimitive(lab(ix,iy,iz-1).s)[u];
                RealTemp c11m1 = _toPrimitive(lab(ix+1,iy+1,iz-1).s)[u];
                RealTemp c01m1 = _toPrimitive(lab(ix,iy+1,iz-1).s)[u];
                RealTemp cm10m1 = _toPrimitive(lab(ix-1,iy,iz-1).s)[u];
                RealTemp cm11m1 = _toPrimitive(lab(ix-1,iy+1,iz-1).s)[u];
                RealTemp c0m1m1 = _toPrimitive(lab(ix,iy-1,iz-1).s)[u];
                RealTemp cm1m1m1 = _toPrimitive(lab(ix-1,iy-1,iz-1).s)[u];
                RealTemp c1m1m1 = _toPrimitive(lab(ix+1,iy-1,iz-1).s)[u];
                
                const RealTemp deriv_u_OnBoundaryX[6] =
                {
                    (4*c100-4*c000+2*c1m10-2*c0m10+2*c101-2*c001+c1m11-c0m11) + (2*c110-2*c010+c111-c011) + (2*c10m1-2*c00m1+c1m1m1-c0m1m1) + (c11m1-c01m1),
                    (2*c100+2*c110+c101+c111) + (-2*cm100-2*cm110-cm101-cm111) + (c10m1+c11m1) + (-cm10m1-cm11m1),
                    (4*c000-4*cm100+2*c010-2*cm110+2*c001-2*cm101+c011-cm111) + (2*c0m10-2*cm1m10+c0m11-cm1m11) + (2*c00m1-2*cm10m1+c01m1-cm11m1) + (c0m1m1-cm1m1m1),
                    (-2*cm100-2*cm1m10-cm101-cm1m11) + (2*c100+2*c1m10+c101+c1m11) + (-cm10m1-cm1m1m1) + (c10m1+c1m1m1),
                    (2*c100+c110+2*c101+c111) + (-2*cm100-cm110-2*cm101-cm111) + (-cm1m10-cm1m11) + (c1m10+c1m11),                    
                    (2*c100+c110+2*c10m1+c11m1) + (-2*cm100-cm110-2*cm10m1-cm11m1) + (-cm1m10-cm1m1m1) + (c1m10+c1m1m1)					
                };
                
                const RealTemp deriv_u_OnBoundaryY[6] =
                {
                    (-2*c0m10-2*c1m10-c0m11-c1m11) + (2*c010+2*c110+c011+c111) + (-c0m1m1-c1m1m1) + (c01m1+c11m1),                    
                    (4*c010-4*c000+2*c110-2*c100+2*c011-2*c001+c111-c101) + (2*cm110-2*cm100+cm111-cm101) + (2*c01m1-2*c00m1+c11m1-c10m1) + (cm11m1-cm10m1),                    
                    (2*cm110+2*c010+cm111+c011) + (-2*cm1m10-2*c0m10-cm1m11-c0m11) + (cm11m1+c01m1) + (-cm1m1m1-c0m1m1),                    
                    (2*cm100-2*cm1m10+4*c000-4*c0m10+cm101-cm1m11+2*c001-2*c0m11) + (2*c100-2*c1m10+c101-c1m11) + (cm10m1-cm1m1m1+2*c00m1-2*c0m1m1) + (c10m1-c1m1m1),                                        
                    (2*c010+c110+2*c011+c111) + (cm110+cm111) + (-cm1m10-2*c0m10-cm1m11-2*c0m11) + (-c1m10-c1m11),                                        
                    (2*c010+c110+2*c01m1+c11m1) + (cm110+cm11m1) + (-cm1m10-2*c0m10-cm1m1m1-2*c0m1m1) + (-c1m10-c1m1m1)
                };
                
                const RealTemp NormalVectZ[8] =                  
                {
                    c001-c000      + c011-c010 +
                    c101-c100  + c111-c110 ,						
                    cm101-cm100      + cm111-cm110 +
                    c001-c000  + c011-c010 ,						
                    cm1m11-cm1m10      + cm101-cm100 +
                    c0m11-c0m10  + c001-c000 ,						
                    c0m11-c0m10	     + c001-c000 +
                    c1m11-c1m10	 + c101-c100,
                    
                    c000-c00m1      + c010-c01m1 +
                    c100-c10m1  + c110-c11m1 ,						
                    cm100-cm10m1      + cm110-cm11m1 +
                    c000-c00m1  + c010-c01m1 ,						
                    cm1m10-cm1m1m1      + cm100-cm10m1 +
                    c0m10-c0m1m1  + c000-c00m1 ,						
                    c0m10-c0m1m1	     + c000-c00m1 +
                    c1m10-c1m1m1	 + c100-c10m1
                };
                
                const RealTemp deriv_u_OnBoundaryZ[6] =
                {
                    NormalVectZ[3] + NormalVectZ[0] + NormalVectZ[7] + NormalVectZ[4],
                    NormalVectZ[0] + NormalVectZ[1] + NormalVectZ[4] + NormalVectZ[5],
                    NormalVectZ[1] + NormalVectZ[2] + NormalVectZ[5] + NormalVectZ[6],
                    NormalVectZ[2] + NormalVectZ[3] + NormalVectZ[6] + NormalVectZ[7],
                    NormalVectZ[0] + NormalVectZ[1] + NormalVectZ[2] + NormalVectZ[3],
                    NormalVectZ[4] + NormalVectZ[5] + NormalVectZ[6] + NormalVectZ[7]
                };
                
                c000 = _toPrimitive(lab(ix,iy,iz).s)[v];
                c100 = _toPrimitive(lab(ix+1,iy,iz).s)[v];
                c110 = _toPrimitive(lab(ix+1,iy+1,iz).s)[v];
                c010 = _toPrimitive(lab(ix,iy+1,iz).s)[v];
                c101 = _toPrimitive(lab(ix+1,iy,iz+1).s)[v];
                c001 = _toPrimitive(lab(ix,iy,iz+1).s)[v];
                c111 = _toPrimitive(lab(ix+1,iy+1,iz+1).s)[v];
                c011 = _toPrimitive(lab(ix,iy+1,iz+1).s)[v];
                cm100 = _toPrimitive(lab(ix-1,iy,iz).s)[v];
                cm110 = _toPrimitive(lab(ix-1,iy+1,iz).s)[v];
                cm101 = _toPrimitive(lab(ix-1,iy,iz+1).s)[v];
                cm111 = _toPrimitive(lab(ix-1,iy+1,iz+1).s)[v];
                c0m10 = _toPrimitive(lab(ix,iy-1,iz).s)[v];
                cm1m10 = _toPrimitive(lab(ix-1,iy-1,iz).s)[v];
                c0m11 = _toPrimitive(lab(ix,iy-1,iz+1).s)[v];
                cm1m11 = _toPrimitive(lab(ix-1,iy-1,iz+1).s)[v];
                c1m10 = _toPrimitive(lab(ix+1,iy-1,iz).s)[v];
                c1m11 = _toPrimitive(lab(ix+1,iy-1,iz+1).s)[v];
                c10m1 = _toPrimitive(lab(ix+1,iy,iz-1).s)[v];
                c00m1 = _toPrimitive(lab(ix,iy,iz-1).s)[v];
                c11m1 = _toPrimitive(lab(ix+1,iy+1,iz-1).s)[v];
                c01m1 = _toPrimitive(lab(ix,iy+1,iz-1).s)[v];
                cm10m1 = _toPrimitive(lab(ix-1,iy,iz-1).s)[v];
                cm11m1 = _toPrimitive(lab(ix-1,iy+1,iz-1).s)[v];
                c0m1m1 = _toPrimitive(lab(ix,iy-1,iz-1).s)[v];
                cm1m1m1 = _toPrimitive(lab(ix-1,iy-1,iz-1).s)[v];
                c1m1m1 = _toPrimitive(lab(ix+1,iy-1,iz-1).s)[v];
                
                
                const RealTemp deriv_v_OnBoundaryX[6] =
                {
                    (4*c100-4*c000+2*c1m10-2*c0m10+2*c101-2*c001+c1m11-c0m11) + (2*c110-2*c010+c111-c011) + (2*c10m1-2*c00m1+c1m1m1-c0m1m1) + (c11m1-c01m1),
                    (2*c100+2*c110+c101+c111) + (-2*cm100-2*cm110-cm101-cm111) + (c10m1+c11m1) + (-cm10m1-cm11m1),
                    (4*c000-4*cm100+2*c010-2*cm110+2*c001-2*cm101+c011-cm111) + (2*c0m10-2*cm1m10+c0m11-cm1m11) + (2*c00m1-2*cm10m1+c01m1-cm11m1) + (c0m1m1-cm1m1m1),
                    (-2*cm100-2*cm1m10-cm101-cm1m11) + (2*c100+2*c1m10+c101+c1m11) + (-cm10m1-cm1m1m1) + (c10m1+c1m1m1),
                    (2*c100+c110+2*c101+c111) + (-2*cm100-cm110-2*cm101-cm111) + (-cm1m10-cm1m11) + (c1m10+c1m11),                    
                    (2*c100+c110+2*c10m1+c11m1) + (-2*cm100-cm110-2*cm10m1-cm11m1) + (-cm1m10-cm1m1m1) + (c1m10+c1m1m1)					
                };
                
                const RealTemp deriv_v_OnBoundaryY[6] =
                {
                    (-2*c0m10-2*c1m10-c0m11-c1m11) + (2*c010+2*c110+c011+c111) + (-c0m1m1-c1m1m1) + (c01m1+c11m1),                    
                    (4*c010-4*c000+2*c110-2*c100+2*c011-2*c001+c111-c101) + (2*cm110-2*cm100+cm111-cm101) + (2*c01m1-2*c00m1+c11m1-c10m1) + (cm11m1-cm10m1),                    
                    (2*cm110+2*c010+cm111+c011) + (-2*cm1m10-2*c0m10-cm1m11-c0m11) + (cm11m1+c01m1) + (-cm1m1m1-c0m1m1),                    
                    (2*cm100-2*cm1m10+4*c000-4*c0m10+cm101-cm1m11+2*c001-2*c0m11) + (2*c100-2*c1m10+c101-c1m11) + (cm10m1-cm1m1m1+2*c00m1-2*c0m1m1) + (c10m1-c1m1m1),                                        
                    (2*c010+c110+2*c011+c111) + (cm110+cm111) + (-cm1m10-2*c0m10-cm1m11-2*c0m11) + (-c1m10-c1m11),                                        
                    (2*c010+c110+2*c01m1+c11m1) + (cm110+cm11m1) + (-cm1m10-2*c0m10-cm1m1m1-2*c0m1m1) + (-c1m10-c1m1m1)
                };
                
                const RealTemp vNormalVectZ[8] =                  
                {
                    c001-c000      + c011-c010 +
                    c101-c100  + c111-c110 ,						
                    cm101-cm100      + cm111-cm110 +
                    c001-c000  + c011-c010 ,						
                    cm1m11-cm1m10      + cm101-cm100 +
                    c0m11-c0m10  + c001-c000 ,						
                    c0m11-c0m10	     + c001-c000 +
                    c1m11-c1m10	 + c101-c100,
                    
                    c000-c00m1      + c010-c01m1 +
                    c100-c10m1  + c110-c11m1 ,						
                    cm100-cm10m1      + cm110-cm11m1 +
                    c000-c00m1  + c010-c01m1 ,						
                    cm1m10-cm1m1m1      + cm100-cm10m1 +
                    c0m10-c0m1m1  + c000-c00m1 ,						
                    c0m10-c0m1m1	     + c000-c00m1 +
                    c1m10-c1m1m1	 + c100-c10m1
                };
                
                const RealTemp deriv_v_OnBoundaryZ[6] =
                {
                    vNormalVectZ[3] + vNormalVectZ[0] + vNormalVectZ[7] + vNormalVectZ[4],
                    vNormalVectZ[0] + vNormalVectZ[1] + vNormalVectZ[4] + vNormalVectZ[5],
                    vNormalVectZ[1] + vNormalVectZ[2] + vNormalVectZ[5] + vNormalVectZ[6],
                    vNormalVectZ[2] + vNormalVectZ[3] + vNormalVectZ[6] + vNormalVectZ[7],
                    vNormalVectZ[0] + vNormalVectZ[1] + vNormalVectZ[2] + vNormalVectZ[3],
                    vNormalVectZ[4] + vNormalVectZ[5] + vNormalVectZ[6] + vNormalVectZ[7]
                };
                
                c000 = _toPrimitive(lab(ix,iy,iz).s)[w];
                c100 = _toPrimitive(lab(ix+1,iy,iz).s)[w];
                c110 = _toPrimitive(lab(ix+1,iy+1,iz).s)[w];
                c010 = _toPrimitive(lab(ix,iy+1,iz).s)[w];
                c101 = _toPrimitive(lab(ix+1,iy,iz+1).s)[w];
                c001 = _toPrimitive(lab(ix,iy,iz+1).s)[w];
                c111 = _toPrimitive(lab(ix+1,iy+1,iz+1).s)[w];
                c011 = _toPrimitive(lab(ix,iy+1,iz+1).s)[w];
                cm100 = _toPrimitive(lab(ix-1,iy,iz).s)[w];
                cm110 = _toPrimitive(lab(ix-1,iy+1,iz).s)[w];
                cm101 = _toPrimitive(lab(ix-1,iy,iz+1).s)[w];
                cm111 = _toPrimitive(lab(ix-1,iy+1,iz+1).s)[w];
                c0m10 = _toPrimitive(lab(ix,iy-1,iz).s)[w];
                cm1m10 = _toPrimitive(lab(ix-1,iy-1,iz).s)[w];
                c0m11 = _toPrimitive(lab(ix,iy-1,iz+1).s)[w];
                cm1m11 = _toPrimitive(lab(ix-1,iy-1,iz+1).s)[w];
                c1m10 = _toPrimitive(lab(ix+1,iy-1,iz).s)[w];
                c1m11 = _toPrimitive(lab(ix+1,iy-1,iz+1).s)[w];
                c10m1 = _toPrimitive(lab(ix+1,iy,iz-1).s)[w];
                c00m1 = _toPrimitive(lab(ix,iy,iz-1).s)[w];
                c11m1 = _toPrimitive(lab(ix+1,iy+1,iz-1).s)[w];
                c01m1 = _toPrimitive(lab(ix,iy+1,iz-1).s)[w];
                cm10m1 = _toPrimitive(lab(ix-1,iy,iz-1).s)[w];
                cm11m1 = _toPrimitive(lab(ix-1,iy+1,iz-1).s)[w];
                c0m1m1 = _toPrimitive(lab(ix,iy-1,iz-1).s)[w];
                cm1m1m1 = _toPrimitive(lab(ix-1,iy-1,iz-1).s)[w];
                c1m1m1 = _toPrimitive(lab(ix+1,iy-1,iz-1).s)[w];
                
                const RealTemp deriv_w_OnBoundaryX[6] =
                {
                    (4*c100-4*c000+2*c1m10-2*c0m10+2*c101-2*c001+c1m11-c0m11) + (2*c110-2*c010+c111-c011) + (2*c10m1-2*c00m1+c1m1m1-c0m1m1) + (c11m1-c01m1),
                    (2*c100+2*c110+c101+c111) + (-2*cm100-2*cm110-cm101-cm111) + (c10m1+c11m1) + (-cm10m1-cm11m1),
                    (4*c000-4*cm100+2*c010-2*cm110+2*c001-2*cm101+c011-cm111) + (2*c0m10-2*cm1m10+c0m11-cm1m11) + (2*c00m1-2*cm10m1+c01m1-cm11m1) + (c0m1m1-cm1m1m1),
                    (-2*cm100-2*cm1m10-cm101-cm1m11) + (2*c100+2*c1m10+c101+c1m11) + (-cm10m1-cm1m1m1) + (c10m1+c1m1m1),
                    (2*c100+c110+2*c101+c111) + (-2*cm100-cm110-2*cm101-cm111) + (-cm1m10-cm1m11) + (c1m10+c1m11),                    
                    (2*c100+c110+2*c10m1+c11m1) + (-2*cm100-cm110-2*cm10m1-cm11m1) + (-cm1m10-cm1m1m1) + (c1m10+c1m1m1)					
                };
                
                const RealTemp deriv_w_OnBoundaryY[6] =
                {
                    (-2*c0m10-2*c1m10-c0m11-c1m11) + (2*c010+2*c110+c011+c111) + (-c0m1m1-c1m1m1) + (c01m1+c11m1),                    
                    (4*c010-4*c000+2*c110-2*c100+2*c011-2*c001+c111-c101) + (2*cm110-2*cm100+cm111-cm101) + (2*c01m1-2*c00m1+c11m1-c10m1) + (cm11m1-cm10m1),                    
                    (2*cm110+2*c010+cm111+c011) + (-2*cm1m10-2*c0m10-cm1m11-c0m11) + (cm11m1+c01m1) + (-cm1m1m1-c0m1m1),                    
                    (2*cm100-2*cm1m10+4*c000-4*c0m10+cm101-cm1m11+2*c001-2*c0m11) + (2*c100-2*c1m10+c101-c1m11) + (cm10m1-cm1m1m1+2*c00m1-2*c0m1m1) + (c10m1-c1m1m1),                                        
                    (2*c010+c110+2*c011+c111) + (cm110+cm111) + (-cm1m10-2*c0m10-cm1m11-2*c0m11) + (-c1m10-c1m11),                                        
                    (2*c010+c110+2*c01m1+c11m1) + (cm110+cm11m1) + (-cm1m10-2*c0m10-cm1m1m1-2*c0m1m1) + (-c1m10-c1m1m1)
                };
                
                const RealTemp wNormalVectZ[8] =                  
                {
                    c001-c000      + c011-c010 +
                    c101-c100  + c111-c110 ,						
                    cm101-cm100      + cm111-cm110 +
                    c001-c000  + c011-c010 ,						
                    cm1m11-cm1m10      + cm101-cm100 +
                    c0m11-c0m10  + c001-c000 ,						
                    c0m11-c0m10	     + c001-c000 +
                    c1m11-c1m10	 + c101-c100,
                    
                    c000-c00m1      + c010-c01m1 +
                    c100-c10m1  + c110-c11m1 ,						
                    cm100-cm10m1      + cm110-cm11m1 +
                    c000-c00m1  + c010-c01m1 ,						
                    cm1m10-cm1m1m1      + cm100-cm10m1 +
                    c0m10-c0m1m1  + c000-c00m1 ,						
                    c0m10-c0m1m1	     + c000-c00m1 +
                    c1m10-c1m1m1	 + c100-c10m1
                };
                
                const RealTemp deriv_w_OnBoundaryZ[6] =
                {
                    wNormalVectZ[3] + wNormalVectZ[0] + wNormalVectZ[7] + wNormalVectZ[4],
                    wNormalVectZ[0] + wNormalVectZ[1] + wNormalVectZ[4] + wNormalVectZ[5],
                    wNormalVectZ[1] + wNormalVectZ[2] + wNormalVectZ[5] + wNormalVectZ[6],
                    wNormalVectZ[2] + wNormalVectZ[3] + wNormalVectZ[6] + wNormalVectZ[7],
                    wNormalVectZ[0] + wNormalVectZ[1] + wNormalVectZ[2] + wNormalVectZ[3],
                    wNormalVectZ[4] + wNormalVectZ[5] + wNormalVectZ[6] + wNormalVectZ[7]
                };                
                
                RealTemp tau[6][9];
                for(int edge=0;edge<6;edge++)
                {
                    tau[edge][0] = (4./3*deriv_u_OnBoundaryX[edge])-2./3*(deriv_v_OnBoundaryY[edge]+deriv_w_OnBoundaryZ[edge]);
                    tau[edge][1] = deriv_v_OnBoundaryX[edge]+deriv_u_OnBoundaryY[edge];
                    tau[edge][2] = deriv_w_OnBoundaryX[edge]+deriv_u_OnBoundaryZ[edge];
                    
                    tau[edge][3] =  tau[edge][1];
                    tau[edge][4] = (4./3*deriv_v_OnBoundaryY[edge])-2./3*(deriv_u_OnBoundaryX[edge]+deriv_w_OnBoundaryZ[edge]);
                    tau[edge][5] = deriv_w_OnBoundaryY[edge]+deriv_v_OnBoundaryZ[edge];
                    
                    tau[edge][6] =  tau[edge][2];
                    tau[edge][7] =  tau[edge][5];
                    tau[edge][8] = (4./3*deriv_w_OnBoundaryZ[edge])-2./3*(deriv_u_OnBoundaryX[edge]+deriv_v_OnBoundaryY[edge]);
                }
                
                RealTemp muOnBoundary[6];
                for(int i=1;i<4;i++)
                {
                    muOnBoundary[0] = 0.5*(_project(_toPrimitive(lab(ix,iy,iz).s)[5], _toPrimitive(lab(ix,iy,iz).s)[0])+_project(_toPrimitive(lab(ix+1,iy,iz).s)[5], _toPrimitive(lab(ix+1,iy,iz).s)[0]));
                    muOnBoundary[1] = 0.5*(_project(_toPrimitive(lab(ix,iy,iz).s)[5], _toPrimitive(lab(ix,iy,iz).s)[0])+_project(_toPrimitive(lab(ix,iy+1,iz).s)[5], _toPrimitive(lab(ix,iy+1,iz).s)[0]));
                    muOnBoundary[2] = 0.5*(_project(_toPrimitive(lab(ix,iy,iz).s)[5], _toPrimitive(lab(ix,iy,iz).s)[0])+_project(_toPrimitive(lab(ix-1,iy,iz).s)[5], _toPrimitive(lab(ix-1,iy,iz).s)[0]));
                    muOnBoundary[3] = 0.5*(_project(_toPrimitive(lab(ix,iy,iz).s)[5], _toPrimitive(lab(ix,iy,iz).s)[0])+_project(_toPrimitive(lab(ix,iy-1,iz).s)[5], _toPrimitive(lab(ix,iy-1,iz).s)[0]));
                    muOnBoundary[4] = 0.5*(_project(_toPrimitive(lab(ix,iy,iz).s)[5], _toPrimitive(lab(ix,iy,iz).s)[0])+_project(_toPrimitive(lab(ix,iy,iz+1).s)[5], _toPrimitive(lab(ix,iy,iz+1).s)[0]));
                    muOnBoundary[5] = 0.5*(_project(_toPrimitive(lab(ix,iy,iz).s)[5], _toPrimitive(lab(ix,iy,iz).s)[0])+_project(_toPrimitive(lab(ix,iy,iz-1).s)[5], _toPrimitive(lab(ix,iy,iz-1).s)[0]));
                }
                
                const RealTemp divTau[3] =                
                {
                    muOnBoundary[0]*tau[0][0]- muOnBoundary[2]*tau[2][0] + muOnBoundary[1]*tau[1][3]-muOnBoundary[3]*tau[3][3] + muOnBoundary[4]*tau[4][6]-muOnBoundary[5]*tau[5][6],
                    muOnBoundary[0]*tau[0][1]- muOnBoundary[2]*tau[2][1] + muOnBoundary[1]*tau[1][4]-muOnBoundary[3]*tau[3][4] + muOnBoundary[4]*tau[4][7]-muOnBoundary[5]*tau[5][7],
                    muOnBoundary[0]*tau[0][2]- muOnBoundary[2]*tau[2][2] + muOnBoundary[1]*tau[1][5]-muOnBoundary[3]*tau[3][5] + muOnBoundary[4]*tau[4][8]-muOnBoundary[5]*tau[5][8]
                };
                
                //add to rhs of momenta
                {
					o(ix, iy, iz).dsdt.u = a*o(ix, iy, iz).dsdt.u + divTau[0]*dtinvh/(16*h);//a*o(ix, iy, iz).dsdt.u - dtinvh/(16*h);//*divTau[0];
                    o(ix, iy, iz).dsdt.v = a*o(ix, iy, iz).dsdt.v + dtinvh/(16*h)*divTau[1];
                    o(ix, iy, iz).dsdt.w = a*o(ix, iy, iz).dsdt.w + dtinvh/(16*h)*divTau[2];
                }
                
                RealTemp uOnBoundary[6][3];
                for(int i=1;i<4;i++)
                {
                    uOnBoundary[0][i-1] = _toPrimitive(lab(ix,iy,iz).s)[i]+_toPrimitive(lab(ix+1,iy,iz).s)[i];
                    uOnBoundary[1][i-1] = _toPrimitive(lab(ix,iy,iz).s)[i]+_toPrimitive(lab(ix,iy+1,iz).s)[i];
                    uOnBoundary[2][i-1] = _toPrimitive(lab(ix,iy,iz).s)[i]+_toPrimitive(lab(ix-1,iy,iz).s)[i];
                    uOnBoundary[3][i-1] = _toPrimitive(lab(ix,iy,iz).s)[i]+_toPrimitive(lab(ix,iy-1,iz).s)[i];
                    uOnBoundary[4][i-1] = _toPrimitive(lab(ix,iy,iz).s)[i]+_toPrimitive(lab(ix,iy,iz+1).s)[i];
                    uOnBoundary[5][i-1] = _toPrimitive(lab(ix,iy,iz).s)[i]+_toPrimitive(lab(ix,iy,iz-1).s)[i];
                }
                
                RealTemp * tau_u = new RealTemp[18];
                memset(tau_u, 0, 18*sizeof(RealTemp));
                for(int i=0; i<6; i++)
                    for(int j=0; j<3; j++)
                        for(int k=0; k<3; k++)
                            tau_u[i+6*j] += tau[i][j*3+k]*uOnBoundary[i][k];
                
                const RealTemp divTau_u = muOnBoundary[0]*tau_u[0]-muOnBoundary[2]*tau_u[2]+muOnBoundary[1]*tau_u[7]-muOnBoundary[3]*tau_u[9]+muOnBoundary[4]*tau_u[16]-muOnBoundary[5]*tau_u[17];
                
                delete [] tau_u;
                
                //add to rhs of energy
                {
                    o(ix, iy, iz).dsdt.s = a*o(ix, iy, iz).dsdt.s + dtinvh/(32*h)*divTau_u;
                }
            }
}
