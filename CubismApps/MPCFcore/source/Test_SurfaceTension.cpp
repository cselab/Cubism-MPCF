/*
 *  Test_SurfaceTension.cpp
 *  MPCFcore
 *
 *  Created by Babak Hejazialhosseini on 2/20/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */
#include <cstdio>
#include <cstring>

#include "common.h"
#include "Test_SurfaceTension.h"

class SurfaceTension_CPP_M1
{
protected:	
	
    const Real dtinvh, sigma, G1, G2, h, a;
    
    Real project(const Real phi);
    void gradient(const Real c[3][2], Real gradc[3]);
    Real gradient_magnitude(const Real gradc[3]);
    void sten_tensor(const Real gradc[3], Real pi_tensor[6][9], const int i);
    Real delta_hs(const Real HSstencil[3][2]);  
    
public:	
	
    SurfaceTension_CPP_M1(const Real a, const Real _dtinvh, const Real G1, const Real G2, const Real _h, const Real _sigma=1): 
	a(a), dtinvh(_dtinvh), G1(G1), G2(G2), h(_h), sigma(_sigma) { }
    
    void compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
                 Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts);
};

void Test_SurfaceTension::_gold(TestLab_S2& lab, Block& block, const Real a, const Real _dtinvh, const Real G1, const Real G2, const Real _h, const Real _sigma)
{	
	SurfaceTension_CPP_M1 gold(a, _dtinvh, G1, G2, _h, _sigma);
	
	const int srcfloats  = sizeof(GP)/sizeof(Real);
	const int rowsrcs = _BLOCKSIZE_+2;
	const int slicesrcs =  (_BLOCKSIZE_+2)*(_BLOCKSIZE_+2);
	
	const int dstfloats =  sizeof(GP)/sizeof(Real);
	const int rowdsts =  _BLOCKSIZE_;
	const int slicedsts = _BLOCKSIZE_*_BLOCKSIZE_;
	
	const Real * const srcfirst = &lab(0,0, 0).s.r;
	Real * const dstfirst = &block(0,0,0).dsdt.r;
	
	gold.compute(srcfirst, srcfloats, rowsrcs, slicesrcs,
				 dstfirst, dstfloats, rowdsts, slicedsts);
}

Real SurfaceTension_CPP_M1::project(const Real phi)
{
	return reconstruct(1, 0, phi, 1/h);
}

void SurfaceTension_CPP_M1::gradient(const Real c[3][2], Real gradc[3])
{
    for(int iDir=0; iDir<3; iDir++)
        gradc[iDir] = c[iDir][1] - c[iDir][0];
}

Real SurfaceTension_CPP_M1::gradient_magnitude(const Real gradc[3])
{
    return sqrt(gradc[0]*gradc[0]+gradc[1]*gradc[1]+gradc[2]*gradc[2]);
}

void SurfaceTension_CPP_M1::sten_tensor(const Real gradc[3], Real pi_tensor[6][9], const int ipoint)
{
    const Real gradmag = gradient_magnitude(gradc);
    
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            pi_tensor[ipoint][j+3*i] = (gradmag==0? 0: ( (i==j ? (Real)1/3:0)*(gradmag*gradmag) - gradc[i]*gradc[j] )*sigma/gradmag );
}

Real SurfaceTension_CPP_M1::delta_hs(const Real HSstencil[3][2])
{
    return sqrt(pow(HSstencil[0][1]-HSstencil[0][0],2)+pow(HSstencil[1][1]-HSstencil[1][0],2)+pow(HSstencil[2][1]-HSstencil[2][0],2));
}

inline void assert_not_nan(const Real * const x, const int N, const int line= __LINE__)
{
	for(int i=0; i<N; ++i)
		assert(!isnan(x[i]));
}

struct AIGP_ST_M1 { Real r, u, v, w, s, l; }; 

void SurfaceTension_CPP_M1::compute(const Real * const srcfirst, const int srcfloats, const int rowsrcs, const int slicesrcs,
                                    Real * const dstfirst, const int dstfloats, const int rowdsts, const int slicedsts)
{
#define GETLS(ix,iy,iz) project((*(AIGP_ST_M1*)(srcfirst + srcfloats*((ix) + (iy)*rowsrcs + (iz)*slicesrcs))).l)
	//#define GETLS(ix,iy,iz) ((ix+1)*(ix+1))	
#define GETU(ix,iy,iz) (*(AIGP_ST_M1*)(srcfirst + srcfloats*((ix) + (iy)*rowsrcs + (iz)*slicesrcs))).u/(*(AIGP_ST_M1*)(srcfirst + srcfloats*((ix) + (iy)*rowsrcs + (iz)*slicesrcs))).r
#define GETV(ix,iy,iz) (*(AIGP_ST_M1*)(srcfirst + srcfloats*((ix) + (iy)*rowsrcs + (iz)*slicesrcs))).v/(*(AIGP_ST_M1*)(srcfirst + srcfloats*((ix) + (iy)*rowsrcs + (iz)*slicesrcs))).r
#define GETW(ix,iy,iz) (*(AIGP_ST_M1*)(srcfirst + srcfloats*((ix) + (iy)*rowsrcs + (iz)*slicesrcs))).w/(*(AIGP_ST_M1*)(srcfirst + srcfloats*((ix) + (iy)*rowsrcs + (iz)*slicesrcs))).r
    
#define OUTPUTU(ix,iy,iz) (*(AIGP_ST_M1*)(dstfirst + dstfloats*((ix) + (iy)*rowdsts + (iz)*slicedsts))).u
#define OUTPUTV(ix,iy,iz) (*(AIGP_ST_M1*)(dstfirst + dstfloats*((ix) + (iy)*rowdsts + (iz)*slicedsts))).v
#define OUTPUTW(ix,iy,iz) (*(AIGP_ST_M1*)(dstfirst + dstfloats*((ix) + (iy)*rowdsts + (iz)*slicedsts))).w
#define OUTPUTS(ix,iy,iz) (*(AIGP_ST_M1*)(dstfirst + dstfloats*((ix) + (iy)*rowdsts + (iz)*slicedsts))).s
    
    for(int iz=0; iz<rowdsts; iz++)
        for(int iy=0; iy<rowdsts; iy++)
            for(int ix=0; ix<rowdsts; ix++)
            {  
				const Real c000 = GETLS(ix,iy,iz);
                const Real c100 = GETLS(ix+1,iy,iz);
                const Real c110 = GETLS(ix+1,iy+1,iz);
                const Real c010 = GETLS(ix,iy+1,iz);
                const Real c101 = GETLS(ix+1,iy,iz+1);
                const Real c001 = GETLS(ix,iy,iz+1);
                const Real c111 = GETLS(ix+1,iy+1,iz+1);
                const Real c011 = GETLS(ix,iy+1,iz+1);
                const Real cm100 = GETLS(ix-1,iy,iz);
                const Real cm110 = GETLS(ix-1,iy+1,iz);
                const Real cm101 = GETLS(ix-1,iy,iz+1);
                const Real cm111 = GETLS(ix-1,iy+1,iz+1);
                const Real c0m10 = GETLS(ix,iy-1,iz);
                const Real cm1m10 = GETLS(ix-1,iy-1,iz);
                const Real c0m11 = GETLS(ix,iy-1,iz+1);
                const Real cm1m11 = GETLS(ix-1,iy-1,iz+1);
                const Real c1m10 = GETLS(ix+1,iy-1,iz);
                const Real c1m11 = GETLS(ix+1,iy-1,iz+1);
                const Real c10m1 = GETLS(ix+1,iy,iz-1);
                const Real c00m1 = GETLS(ix,iy,iz-1);
                const Real c11m1 = GETLS(ix+1,iy+1,iz-1);
                const Real c01m1 = GETLS(ix,iy+1,iz-1);
                const Real cm10m1 = GETLS(ix-1,iy,iz-1);
                const Real cm11m1 = GETLS(ix-1,iy+1,iz-1);
                const Real c0m1m1 = GETLS(ix,iy-1,iz-1);
                const Real cm1m1m1 = GETLS(ix-1,iy-1,iz-1);
                const Real c1m1m1 = GETLS(ix+1,iy-1,iz-1);  
				
				for(int dz=-1; dz<2; ++dz)
					for(int dy=-1; dy<2; ++dy)
						for(int dx=-1; dx<2; ++dx)
							assert(!isnan(GETLS(ix+dx,iy+dy,iz+dz)));
				
                const Real NormalVectOnBoundaryX[6] =
                {
                    (4*c100-4*c000+2*c1m10-2*c0m10+2*c101-2*c001+c1m11-c0m11) + (2*c110-2*c010+c111-c011) + (2*c10m1-2*c00m1+c1m1m1-c0m1m1) + (c11m1-c01m1),
                    (2*c100+2*c110+c101+c111) + (-2*cm100-2*cm110-cm101-cm111) + (c10m1+c11m1) + (-cm10m1-cm11m1),
                    (4*c000-4*cm100+2*c010-2*cm110+2*c001-2*cm101+c011-cm111) + (2*c0m10-2*cm1m10+c0m11-cm1m11) + (2*c00m1-2*cm10m1+c01m1-cm11m1) + (c0m1m1-cm1m1m1),
                    (-2*cm100-2*cm1m10-cm101-cm1m11) + (2*c100+2*c1m10+c101+c1m11) + (-cm10m1-cm1m1m1) + (c10m1+c1m1m1),
                    (2*c100+c110+2*c101+c111) + (-2*cm100-cm110-2*cm101-cm111) + (-cm1m10-cm1m11) + (c1m10+c1m11),                    
                    (2*c100+c110+2*c10m1+c11m1) + (-2*cm100-cm110-2*cm10m1-cm11m1) + (-cm1m10-cm1m1m1) + (c1m10+c1m1m1)					
                };
				
				assert_not_nan(NormalVectOnBoundaryX, sizeof(NormalVectOnBoundaryX)/sizeof(Real));
                
                const Real NormalVectOnBoundaryY[6] =
                {
                    (-2*c0m10-2*c1m10-c0m11-c1m11) + (2*c010+2*c110+c011+c111) + (-c0m1m1-c1m1m1) + (c01m1+c11m1),                    
                    (4*c010-4*c000+2*c110-2*c100+2*c011-2*c001+c111-c101) + (2*cm110-2*cm100+cm111-cm101) + (2*c01m1-2*c00m1+c11m1-c10m1) + (cm11m1-cm10m1),                    
                    (2*cm110+2*c010+cm111+c011) + (-2*cm1m10-2*c0m10-cm1m11-c0m11) + (cm11m1+c01m1) + (-cm1m1m1-c0m1m1),                    
                    (2*cm100-2*cm1m10+4*c000-4*c0m10+cm101-cm1m11+2*c001-2*c0m11) + (2*c100-2*c1m10+c101-c1m11) + (cm10m1-cm1m1m1+2*c00m1-2*c0m1m1) + (c10m1-c1m1m1),                                        
                    (2*c010+c110+2*c011+c111) + (cm110+cm111) + (-cm1m10-2*c0m10-cm1m11-2*c0m11) + (-c1m10-c1m11),                                        
                    (2*c010+c110+2*c01m1+c11m1) + (cm110+cm11m1) + (-cm1m10-2*c0m10-cm1m1m1-2*c0m1m1) + (-c1m10-c1m1m1)
                };
				
				assert_not_nan(NormalVectOnBoundaryY, sizeof(NormalVectOnBoundaryY)/sizeof(Real));
                
                const Real NormalVectZ[8] = 
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
				
				assert_not_nan(NormalVectZ, sizeof(NormalVectZ)/sizeof(Real));
                
                const Real NormalVectOnBoundaryZ[6] =
                {
                    NormalVectZ[3] + NormalVectZ[0] + NormalVectZ[7] + NormalVectZ[4],
                    NormalVectZ[0] + NormalVectZ[1] + NormalVectZ[4] + NormalVectZ[5],
                    NormalVectZ[1] + NormalVectZ[2] + NormalVectZ[5] + NormalVectZ[6],
                    NormalVectZ[2] + NormalVectZ[3] + NormalVectZ[6] + NormalVectZ[7],
                    NormalVectZ[0] + NormalVectZ[1] + NormalVectZ[2] + NormalVectZ[3],
                    NormalVectZ[4] + NormalVectZ[5] + NormalVectZ[6] + NormalVectZ[7]
                };
				
				assert_not_nan(NormalVectOnBoundaryZ, sizeof(NormalVectOnBoundaryZ)/sizeof(Real));
                
                Real pi_tensor[6][9];
                for(int i=0; i<6; i++)
                {
                    const Real gradc[3] = {NormalVectOnBoundaryX[i], NormalVectOnBoundaryY[i], NormalVectOnBoundaryZ[i]};
					assert(!isnan(gradc[0]));
					assert(!isnan(gradc[1]));
					assert(!isnan(gradc[2]));
                    sten_tensor(gradc, pi_tensor, i);
				}
                
                {                    
                    const Real divX = pi_tensor[0][0] - pi_tensor[2][0] + pi_tensor[1][3] - pi_tensor[3][3] + pi_tensor[4][6] - pi_tensor[5][6]; 
                    const Real divY = pi_tensor[0][1] - pi_tensor[2][1] + pi_tensor[1][4] - pi_tensor[3][4] + pi_tensor[4][7] - pi_tensor[5][7]; 
                    const Real divZ = pi_tensor[0][2] - pi_tensor[2][2] + pi_tensor[1][5] - pi_tensor[3][5] + pi_tensor[4][8] - pi_tensor[5][8];
                    
					const Real rhsu = dtinvh*divX/(16*h);
                    const Real rhsv = dtinvh*divY/(16*h);
                    const Real rhsw = dtinvh*divZ/(16*h);                        
					
                    OUTPUTU(ix,iy,iz) = a*OUTPUTU(ix,iy,iz) + rhsu;
                    OUTPUTV(ix,iy,iz) = a*OUTPUTV(ix,iy,iz) + rhsv;
                    OUTPUTW(ix,iy,iz) = a*OUTPUTW(ix,iy,iz) + rhsw;
					
					assert(!isnan(OUTPUTU(ix,iy,iz)));
					assert(!isnan(OUTPUTV(ix,iy,iz)));
					assert(!isnan(OUTPUTW(ix,iy,iz)));
                }  
                
                {                    
                    const Real vel[6][3] = 
                    {
                        GETU(ix+1,iy,iz)+GETU(ix,iy,iz)  , GETV(ix+1,iy,iz)+GETV(ix,iy,iz)  , GETW(ix+1,iy,iz)+GETW(ix,iy,iz),                        
                        GETU(ix,iy+1,iz)+GETU(ix,iy,iz)  , GETV(ix,iy+1,iz)+GETV(ix,iy,iz)  , GETW(ix,iy+1,iz)+GETW(ix,iy,iz),                        
                        GETU(ix,iy,iz)  +GETU(ix-1,iy,iz), GETV(ix,iy,iz)  +GETV(ix-1,iy,iz), GETW(ix,iy,iz)  +GETW(ix-1,iy,iz),
                        GETU(ix,iy,iz)  +GETU(ix,iy-1,iz), GETV(ix,iy,iz)  +GETV(ix,iy-1,iz), GETW(ix,iy,iz)  +GETW(ix,iy-1,iz),                        
                        GETU(ix,iy,iz+1)+GETU(ix,iy,iz)  , GETV(ix,iy,iz+1)+GETV(ix,iy,iz)  , GETW(ix,iy,iz+1)+GETW(ix,iy,iz),                        
                        GETU(ix,iy,iz)  +GETU(ix,iy,iz-1), GETV(ix,iy,iz)  +GETV(ix,iy,iz-1), GETW(ix,iy,iz)  +GETW(ix,iy,iz-1)
                    };
                    
                    Real * pi_tensor_u = new Real[18];
                    memset(pi_tensor_u, 0, 18*sizeof(Real));
                    for(int i=0; i<6; i++)
                        for(int j=0; j<3; j++)
                            for(int k=0; k<3; k++)
                                pi_tensor_u[i+6*j] += pi_tensor[i][j*3+k]*vel[i][k];
					
                    const Real divAll = pi_tensor_u[0]-pi_tensor_u[2]+pi_tensor_u[7]-pi_tensor_u[9]+pi_tensor_u[16]-pi_tensor_u[17]; 
                    delete [] pi_tensor_u;
                    
                    const Real rhss = dtinvh*divAll/(32*h);
					
                    OUTPUTS(ix,iy,iz) = a*OUTPUTS(ix,iy,iz) + rhss;
					
					assert(!isnan(OUTPUTS(ix,iy,iz)));
                }
            }
#undef GETLS
#undef GETU
#undef GETV
#undef GETW
    
#undef OUTPUTU
#undef OUTPUTV
#undef OUTPUTW
#undef OUTPUTS
}
