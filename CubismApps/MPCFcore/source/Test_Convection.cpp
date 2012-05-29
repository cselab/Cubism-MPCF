/*
 *  Test_Convection.cpp
 *  MPCFcore
 *
 *  Created by Diego Rossinelli on 5/17/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <cstring>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>

using namespace std;

#include "TestTypes.h"
#include "Test_Convection.h"

void Test_Convection::_initialize(TestLab& lab, Block& block)
{
	srand48(61651);
	
	for(int iz = -3; iz<_BLOCKSIZE_+3; iz++)
		for(int iy = -3; iy<_BLOCKSIZE_+3; iy++)
		{
			for(int ix = -3; ix<_BLOCKSIZE_+3; ix++)
			{
				memset(&lab(ix, iy, iz), 0, sizeof(GP));
				
				const int a = iy;
				const int b = iz; 
				const int c = ix;
				const double L = _BLOCKSIZE_;
				
				lab(ix, iy, iz).s.r = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L;
				lab(ix, iy, iz).s.u = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L;
				lab(ix, iy, iz).s.v = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L;
				lab(ix, iy, iz).s.w = (40+a*a+(c+3)*L + (a+3)*L*L)/(double)L/(double)L;
				lab(ix, iy, iz).s.s = 100+b/(double)L;
				lab(ix, iy, iz).s.G = (1+c+b*L + a*L*L)/(double)L;

#ifdef _LIQUID_
                lab(ix, iy, iz).s.P = (1+c+b*L + a*L*L)/(double)L;
#endif
			}
		}
	
	for(int iz = 0; iz<_BLOCKSIZE_; iz++)
		for(int iy = 0; iy<_BLOCKSIZE_; iy++)
			for(int ix = 0; ix<_BLOCKSIZE_; ix++)
			{
				memset(&block(ix, iy, iz), 0, sizeof(GP));
				
				block(ix, iy, iz).dsdt.r = 0.1+iz;
				block(ix, iy, iz).dsdt.u = ix;
				block(ix, iy, iz).dsdt.v = iy;
				block(ix, iy, iz).dsdt.w = (ix*iy/(double)_BLOCKSIZE_);
				block(ix, iy, iz).dsdt.s = 1+(iy*iz)/(double)_BLOCKSIZE_;
				block(ix, iy, iz).dsdt.G = -1 +  (iz+ix+iy)/(double)_BLOCKSIZE_;
 
#ifdef _LIQUID_
                block(ix, iy, iz).dsdt.P = -1 +  (iz+ix+iy)/(double)_BLOCKSIZE_;
#endif
			}
}

void Test_Convection::_print(Block& block)
{
	for(int iz = 0; iz<_BLOCKSIZE_; iz++)
		for(int iy = 0; iy<_BLOCKSIZE_; iy++)
			for(int ix = 0; ix<_BLOCKSIZE_; ix++)
			{
				StateVector v = block(ix, iy, iz).dsdt;
#ifndef _LIQUID_
				printf("%d %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e\n", ix+ iy*_BLOCKSIZE_ + iz*_BLOCKSIZE_*_BLOCKSIZE_, v.r, v.u, v.v, v.w, v.s, v.levelset);
#else
                printf("%d %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e %15.15e\n", ix+ iy*_BLOCKSIZE_ + iz*_BLOCKSIZE_*_BLOCKSIZE_, v.r, v.u, v.v, v.w, v.s, v.G, v.P);
#endif
			}
}

inline Real weno_minus(const Real a, const Real b, const Real c, const Real d, const Real e)
{
	const Real is0 = a*(a*(Real)(4./3.)  - b*(Real)(19./3.)  + c*(Real)(11./3.)) + b*(b*(Real)(25./3.)  - c*(Real)(31./3.)) + c*c*(Real)(10./3.);
	const Real is1 = b*(b*(Real)(4./3.)  - c*(Real)(13./3.)  + d*(Real)(5./3.)) + c*(c*(Real)(13./3.)  - d*(Real)(13./3.)) + d*d*(Real)(4./3.);
	const Real is2 = c*(c*(Real)(10./3.)  - d*(Real)(31./3.)  + e*(Real)(11./3.)) + d*(d*(Real)(25./3.)  - e*(Real)(19./3.)) + e*e*(Real)(4./3.);
	
	const Real is0plus = is0+WENOEPS;
	const Real is1plus = is1+WENOEPS;
	const Real is2plus = is2+WENOEPS;
	
	const Real alpha0 = (Real)(0.1)*((Real)1/(is0plus*is0plus));
	const Real alpha1 = (Real)(0.6)*((Real)1/(is1plus*is1plus));
	const Real alpha2 = (Real)(0.3)*((Real)1/(is2plus*is2plus));
	const Real alphasum = alpha0+alpha1+alpha2;
	
	const Real omega0=alpha0 * (((Real)1)/alphasum);
	const Real omega1=alpha1 * (((Real)1)/alphasum);
	const Real omega2= 1-omega0-omega1;
	
	return omega0*((Real)(1.0/3.)*a-(Real)(7./6.)*b+(Real)(11./6.)*c) + omega1*(-(Real)(1./6.)*b+(Real)(5./6.)*c+(Real)(1./3.)*d) + omega2*((Real)(1./3.)*c+(Real)(5./6.)*d-(Real)(1./6.)*e);
}

inline Real weno_plus(const Real b, const Real c, const Real d, const Real e, const Real f)
{	
	const Real is0 = d*(d*(Real)(10./3.) - e*(Real)(31./3.)  + f*(Real)(11./3.)) + e*(e*(Real)(25./3.)  - f*(Real)(19./3.)) +	f*f*(Real)(4./3.);
	const Real is1 = c*(c*(Real)(4./3.) - d*(Real)(13./3.)  + e*(Real)(5./3.)) + d*(d*(Real)(13./3.)  - e*(Real)(13./3.)) +	e*e*(Real)(4./3.);
	const Real is2 = b*(b*(Real)(4./3.) - c*(Real)(19./3.)  + d*(Real)(11./3.)) + c*(c*(Real)(25./3.)  - d*(Real)(31./3.)) +	d*d*(Real)(10./3.);
	
	const Real is0plus = is0+ ((Real)WENOEPS);
	const Real is1plus = is1+ ((Real)WENOEPS);
	const Real is2plus = is2+ ((Real)WENOEPS);
	
	const Real alpha0 = (Real)(0.1)*(((Real)1)/(is0plus*is0plus));
	const Real alpha1 = (Real)(0.6)*(((Real)1)/(is1plus*is1plus));
	const Real alpha2 = (Real)(0.3)*(((Real)1)/(is2plus*is2plus));
	const Real alphasum = alpha0+alpha1+alpha2;
	
	const Real omega0=alpha0 * (((Real)1)/alphasum);
	const Real omega1=alpha1 * (((Real)1)/alphasum);
	const Real omega2= 1-omega0-omega1;
	
	return omega0*((Real)(1./3.)*f-(Real)(7./6.)*e+(Real)(11./6.)*d) + omega1*(-(Real)(1./6.)*e+(Real)(5./6.)*d+(Real)(1./3.)*c) + omega2*((Real)(1./3.)*d+(Real)(5./6.)*c-(Real)(1./6.)*b);
}


void weno(const Real input[6], Real sides[2]) 
{
	for(int i=0; i<6; i++)
		assert(!isnan(input[i]));
	
	
	sides[0] = weno_minus(input[0],input[1],input[2],input[3],input[4]);
	sides[1] = weno_plus(input[1],input[2],input[3],input[4],input[5]);
	
	for(int i=0; i<2; i++)
		assert(!isnan(sides[i]));
}

struct Vec5
{
	Real  data[5];
	
	Vec5& operator=(Vec5 a)
	{
		for(int i=0; i<5; i++) data[i] = a.data[i];
		return *this;
	}
	
	Vec5& operator+=(Vec5& a)
	{
		for(int i=0; i<5; i++) data[i] += a.data[i];
		return *this;
	}
	
	Vec5()
	{
		data[0] = data[1] = data[2] = data[3] = data[4] = 0.;
	}
	
	Vec5(const Vec5& a)
	{
		for(int i=0; i<5; i++) data[i] = a.data[i];
	}
	
	Vec5(Real a, Real b, Real c, Real d, Real e)
	{
		data[0] = a;
		data[1] = b;
		data[2] = c;
		data[3] = d;
		data[4] = e;
	}
	
	const Real& operator[](const int i) const
	{
		return data[i];
	}
	
	Real& operator[](const int i)
	{
		return data[i];
	}
	
	inline void print() const
	{
		printf("Values of EU are : %e %e %e %e %e\n",data[0],data[1],data[2],data[3],data[4]);
	}
};


void setV5(Vec5& v, const Real a, const Real b, const Real c, const Real d, const Real e)
{
	v.data[0] = a;
	v.data[1] = b;
	v.data[2] = c;
	v.data[3] = d;
	v.data[4] = e;
}

void print(const Vec5& p ){p.print();}

Vec5 operator+(Vec5 a, Vec5 b)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a.data[i] + b.data[i];
	return r;
}

Vec5 operator+( Vec5 b, Real a)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a + b.data[i];
	return r;
}

Vec5 operator-(Vec5 a, Vec5 b)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a.data[i] - b.data[i];
	return r;
}

Vec5 operator-(Vec5 a, Real b)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a.data[i] - b;
	return r;
}

Vec5 operator*(Real a, Vec5 b)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a * b.data[i];
	return r;
}

Vec5 operator*(Vec5 b, Real a)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a * b.data[i];
	return r;
}

Vec5 operator*(Vec5 a, Vec5 b)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a.data[i] * b.data[i];
	return r;
}

Vec5 operator/(Vec5 a, Vec5 b)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a.data[i] / b.data[i];
	return r;
}

Vec5 operator/(Real a, Vec5 b)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = a/ b.data[i];
	return r;
}

Vec5 operator/(Vec5 b, Real a)
{
	Vec5 r;
	for(int i=0; i<5; i++) r.data[i] = b.data[i]/a;
	return r;
}

template <typename V5=Vec5> 
struct FluxBuilder_HLLE
{
	inline Real _CharacteristicVelocityLeft(const Real momentum[2], const Real density[2], const Real a[2]) const
	{
		return  min(momentum[0]/density[0] - a[0], momentum[1]/density[1] - a[1]);
	}
	
	inline Real _CharacteristicVelocityRight(const Real momentum[2], const Real density[2], const Real a[2]) const
	{
		return  max(momentum[0]/density[0] + a[0], momentum[1]/density[1] + a[1]);	
	}
	
	
	inline float _soundEOS(const float P, const float density, const float gamma, const bool bisIdealGas, const float corr_gamma) const
	{
		float result= (bisIdealGas)? sqrtf(gamma*max(P/density,(float)0.f)):
		sqrtf(gamma*max((P+gamma*corr_gamma)/density,(float)0.f));
		
		assert(result != 0);
		
		return result;
	}
	
	inline double _soundEOS(const double P, const double density, const double gamma, const bool bisIdealGas, const double corr_gamma) const
	{
		double result= (bisIdealGas)? sqrt(gamma*max(P/density,(double)0.0)):
		sqrt(gamma*max((P+gamma*corr_gamma)/density,(double)0.0));
		
		assert(result != 0);
		
		return result;
	}
	
	Real _heavyside(Real phi, Real h) const
	{
		const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/h)));
		const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
		const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
		return x<0 ? val_xneg : val_xpos;
	}
	
	inline Real _getGamma(Real phi, Real h)const
	{		
		const Real HS = _heavyside(phi, h);
		assert(HS>=0);
		assert(HS<=1);
		assert(gamma1>0);
		assert(gamma2>0);
		assert(gamma1*HS + gamma2*(((Real)1)-HS)>=1);
		return (gamma1*HS + gamma2*(((Real)1)-HS));
	}
	
	Real h, gamma1, gamma2;
	FluxBuilder_HLLE(Real h, Real gamma1, Real gamma2): h(h), gamma1(gamma1), gamma2(gamma2){}
	
	inline void HLLE(const Real fLEFT, const Real fRIGHT, const Real uLEFT, const Real uRIGHT, const Real aLEFT, const Real aRIGHT,  Real& r) const
	{
		
		Real candidates[3] = {
			fLEFT, 
			fRIGHT,
			(aRIGHT*fLEFT - aLEFT*fRIGHT + aLEFT*aRIGHT*(uRIGHT - uLEFT))*(1./(aRIGHT- aLEFT))
		};
		
		const bool flags[2] = {aLEFT>0.0, aRIGHT<0.0};
		
		r = flags[0] ? candidates[0] : (flags[1] ? candidates[1] : candidates[2]);
	}
	
	inline void HLLE5(const V5& fLEFT, const V5& fRIGHT, const V5& uLEFT, const V5& uRIGHT, const V5& aLEFT, const V5& aRIGHT, V5& r) const
	{
		V5 candidates[3] = {
			fLEFT, 
			fRIGHT,
			(aRIGHT*fLEFT - aLEFT*fRIGHT + aLEFT*aRIGHT*(uRIGHT - uLEFT))*(1./(aRIGHT- aLEFT))
		};
		
		bool flags[2][5];
		for(int i=0; i<5; i++) flags[0][i] = (int)(aLEFT[i]>0.0);
		for(int i=0; i<5; i++) flags[1][i] = (int)(aRIGHT[i]<0.0);
		
		for(int i=0; i<5; i++) r[i] = flags[0][i] ? candidates[0][i] : (flags[1][i] ? candidates[1][i] : candidates[2][i]);
	}
	
	template<int iDir>
	void buildUsingPrimitives(const Real U[6][2], Real FLUX[6]) 
	{
		V5 EU[2];		
		for(int j=0; j<5; j++)
			for(int i=0;i<2; i++)
				EU[i][j] = U[j][i];
		
		Real LS[2];
		for(int i=0;i<2; i++)
			LS[i] = U[5][i];
		
		V5 EUflux;
		Real LSflux;
		buildUsingPrimitives<iDir>(EU, LS, EUflux, LSflux);
		assert(!isnan(EUflux[0]));
		assert(!isnan(EUflux[1]));
		assert(!isnan(EUflux[2]));
		assert(!isnan(EUflux[3]));
		assert(!isnan(EUflux[4]));
		assert(!isnan(LSflux));
		
		for(int j=0; j<5; j++)
			FLUX[j] = EUflux[j];
		
		FLUX[5] = LSflux;
	}
	
	template<int iDir>
	void buildUsingPrimitives(const V5 EU[2], const Real LS[2], V5& EUflux, Real& LSflux) const
	{
		static const int r = 0;
		static const int u = 1;
		static const int v = 2;
		static const int w = 3;
		static const int p = 4;
		
		V5 uLEFT, uRIGHT, aLEFT, aRIGHT, fLEFT, fRIGHT;
		
		const Real Vel[3][2] = {
			EU[0][u] + 0, 
			EU[1][u] + 0, 			
			EU[0][v],  EU[1][v],
			EU[0][w],  EU[1][w]
		};
		
		const Real R[2] = { EU[0][r], EU[1][r] }; 
		
		const Real gamma[2] = {
		 	_getGamma(LS[0], h),
			_getGamma(LS[1], h)
		};
		
		const Real P[2] = {
			EU[0][p],
			EU[1][p]
		};
		
#ifndef _WATERBUBBLE_			
		const Real SOS[2] = { 
			_soundEOS(P[0], R[0], gamma[0], 1, (Real)0),  
			_soundEOS(P[1], R[1], gamma[1], 1, (Real)0) 
		};
#else
		const Real SOS[2] = { 
			_soundEOS(P[0], R[0], gamma[0], (LS[0]<0.0? 0:1), (Real)(LS[0]<0.0? _Pc_:0.0)),  
			_soundEOS(P[1], R[1], gamma[1], (LS[1]<0.0? 0:1), (Real)(LS[1]<0.0? _Pc_:0.0)) 
		};
#endif
		
		const Real mom[3][2] = {
			Vel[0][0]*R[0],	Vel[0][1]*R[1],
			Vel[1][0]*R[0],	Vel[1][1]*R[1],
			Vel[2][0]*R[0],	Vel[2][1]*R[1]
		};
		
		const Real CV[2] = { 
			_CharacteristicVelocityLeft(mom[iDir], R , SOS),
			_CharacteristicVelocityRight(mom[iDir], R , SOS)
		};
		
		const Real E[2] = {
			P[0]/(gamma[0]-(Real)1)+ ((Real)0.5)*R[0]*(Vel[0][0]*Vel[0][0] + Vel[1][0]*Vel[1][0] + Vel[2][0]*Vel[2][0] ),
			P[1]/(gamma[1]-(Real)1)+ ((Real)0.5)*R[1]*(Vel[0][1]*Vel[0][1] + Vel[1][1]*Vel[1][1] + Vel[2][1]*Vel[2][1] )	
		};
		
		setV5(fLEFT ,   Vel[iDir][0]*R[0], 
			  Vel[0][0]*Vel[iDir][0]*R[0] + (iDir==0 ? P[0] : 0.0), 
			  Vel[1][0]*Vel[iDir][0]*R[0] + (iDir==1 ? P[0] : 0.0),
			  Vel[2][0]*Vel[iDir][0]*R[0] + (iDir==2 ? P[0] : 0.0),
#ifndef _WATERBUBBLE_			  
			  Vel[iDir][0]*(P[0]+E[0]));
#else
		Vel[iDir][0]*((P[0]+gamma[0]*( LS[0]<0.0? _Pc_:0.0) )*(1.+1./(gamma[0]-1.) ) + 0.5*R[0]*(Vel[0][0]*Vel[0][0] + Vel[1][0]*Vel[1][0] + Vel[2][0]*Vel[2][0] ) ) );
#endif
		
		
		setV5(fRIGHT ,  Vel[iDir][1]*R[1], 
			  Vel[0][1]*Vel[iDir][1]*R[1] + (iDir==0 ? P[1] : 0.0), 
			  Vel[1][1]*Vel[iDir][1]*R[1] + (iDir==1 ? P[1] : 0.0),
			  Vel[2][1]*Vel[iDir][1]*R[1] + (iDir==2 ? P[1] : 0.0),
#ifndef _WATERBUBBLE_
			  Vel[iDir][1]*(P[1]+E[1]));
#else
		Vel[iDir][1]*((P[1]+gamma[1]*( LS[1]<0.0? _Pc_:0.0) )*(1.+1./(gamma[1]-1.) ) + 0.5*R[1]*(Vel[0][1]*Vel[0][1] + Vel[1][1]*Vel[1][1] + Vel[2][1]*Vel[2][1] ) ) );
#endif
		
#ifndef _WATERBUBBLE_
		setV5(uLEFT ,  R[0], R[0]*Vel[0][0], R[0]*Vel[1][0], R[0]*Vel[2][0], P[0]/(gamma[0]-1.) + 0.5*R[0]*(Vel[0][0]*Vel[0][0] + Vel[1][0]*Vel[1][0] + Vel[2][0]*Vel[2][0]) );
		setV5(uRIGHT , R[1], R[1]*Vel[0][1], R[1]*Vel[1][1], R[1]*Vel[2][1], P[1]/(gamma[1]-1.) + 0.5*R[1]*(Vel[0][1]*Vel[0][1] + Vel[1][1]*Vel[1][1] + Vel[2][1]*Vel[2][1]) );
#else
		setV5(uLEFT ,  R[0], R[0]*Vel[0][0], R[0]*Vel[1][0], R[0]*Vel[2][0], (P[0]+gamma[0]*( LS[0]<0.0? _Pc_:0.0) )/(gamma[0]-1.) + 0.5*R[0]*(Vel[0][0]*Vel[0][0] + Vel[1][0]*Vel[1][0] + Vel[2][0]*Vel[2][0]) );
		setV5(uRIGHT , R[1], R[1]*Vel[0][1], R[1]*Vel[1][1], R[1]*Vel[2][1], (P[1]+gamma[1]*( LS[1]<0.0? _Pc_:0.0) )/(gamma[1]-1.) + 0.5*R[1]*(Vel[0][1]*Vel[0][1] + Vel[1][1]*Vel[1][1] + Vel[2][1]*Vel[2][1]) );
#endif
		setV5(aLEFT ,  CV[0], CV[0], CV[0], CV[0], CV[0]);
		
		setV5(aRIGHT , CV[1], CV[1], CV[1], CV[1], CV[1]);
		HLLE5(fLEFT, fRIGHT, uLEFT, uRIGHT, aLEFT, aRIGHT, EUflux);
		
		assert(!isnan(EUflux[0]));
		assert(!isnan(EUflux[1]));
		assert(!isnan(EUflux[2]));
		assert(!isnan(EUflux[3]));
		assert(!isnan(EUflux[4]));
		
		HLLE(Vel[iDir][0]*LS[0], Vel[iDir][1]*LS[1], LS[0], LS[1], CV[0], CV[1], LSflux);		
		assert(!isnan(LSflux));
		
#if defined(_SINGLE_PHASE_)
		LSflux= 0.;
#endif
	}
};


void Test_Convection::_gold(TestLab& lab, Block& block, const Real h_gridpoint)
{
	for(int iz = 0; iz<_BLOCKSIZE_; iz++)
		for(int iy = 0; iy<_BLOCKSIZE_; iy++)
			for(int ix = 0; ix<_BLOCKSIZE_; ix++)
				assert(lab(ix, iy, iz).s.r>0);
	
	for(int iz = 0; iz<_BLOCKSIZE_; iz++)
		for(int iy = 0; iy<_BLOCKSIZE_; iy++)
			for(int ix = 0; ix<_BLOCKSIZE_; ix++)
			{
				const Real minuslambda = -1.0/h_gridpoint;
				
				Real out[6];
				
				Real WE_LEFT[6][6], WE_RIGHT[6][6];
				Real NS_DOWN[6][6], NS_UP[6][6];
				Real BF_BACK[6][6], BF_FRONT[6][6];
				
				for (int i=-3; i<3; i++)
				{
					_toPrimitive(lab(ix+i, iy, iz), out);
					
					WE_LEFT[0][i+3] = out[0];
					WE_LEFT[1][i+3] = out[1];
					WE_LEFT[2][i+3] = out[2];
					WE_LEFT[3][i+3] = out[3];
					WE_LEFT[4][i+3] = out[4];
					WE_LEFT[5][i+3] = out[5];
					
					_toPrimitive(lab(ix+i+1, iy, iz), out);
					
					WE_RIGHT[0][i+3] = out[0];
					WE_RIGHT[1][i+3] = out[1];
					WE_RIGHT[2][i+3] = out[2];
					WE_RIGHT[3][i+3] = out[3];
					WE_RIGHT[4][i+3] = out[4];
					WE_RIGHT[5][i+3] = out[5];
					
					_toPrimitive(lab(ix, iy+i, iz), out);
					
					NS_DOWN[0][i+3] = out[0];
					NS_DOWN[1][i+3] = out[1];
					NS_DOWN[2][i+3] = out[2];
					NS_DOWN[3][i+3] = out[3];
					NS_DOWN[4][i+3] = out[4];
					NS_DOWN[5][i+3] = out[5];
					
					_toPrimitive(lab(ix, iy+i+1, iz), out);
					
					NS_UP[0][i+3] = out[0];
					NS_UP[1][i+3] = out[1];
					NS_UP[2][i+3] = out[2];
					NS_UP[3][i+3] = out[3];
					NS_UP[4][i+3] = out[4];
					NS_UP[5][i+3] = out[5];
					
					_toPrimitive(lab(ix, iy, iz+i), out);
					
					BF_BACK[0][i+3] = out[0];
					BF_BACK[1][i+3] = out[1];
					BF_BACK[2][i+3] = out[2];
					BF_BACK[3][i+3] = out[3];
					BF_BACK[4][i+3] = out[4];
					BF_BACK[5][i+3] = out[5];
					
					_toPrimitive(lab(ix, iy, iz+i+1), out);
					
					BF_FRONT[0][i+3] = out[0];
					BF_FRONT[1][i+3] = out[1];
					BF_FRONT[2][i+3] = out[2];
					BF_FRONT[3][i+3] = out[3];
					BF_FRONT[4][i+3] = out[4];
					BF_FRONT[5][i+3] = out[5];
				}
				
				Real SIDES_WE_LEFT[6][2], SIDES_WE_RIGHT[6][2];
				Real SIDES_NS_DOWN[6][2], SIDES_NS_UP[6][2];
				Real SIDES_BF_BACK[6][2], SIDES_BF_FRONT[6][2];
				
				for(int i=0; i<6; i++)
				{
					weno(WE_LEFT[i] , SIDES_WE_LEFT[i] );
					weno(WE_RIGHT[i], SIDES_WE_RIGHT[i]);
					weno(NS_DOWN[i] , SIDES_NS_DOWN[i] );
					weno(NS_UP[i]   , SIDES_NS_UP[i]   );
					weno(BF_BACK[i] , SIDES_BF_BACK[i] );
					weno(BF_FRONT[i], SIDES_BF_FRONT[i]);
				}
				
				Real fluxN[6], fluxS[6], fluxW[6], fluxE[6], fluxB[6], fluxF[6];
				
				FluxBuilder_HLLE<> flux_builder(h_gridpoint, gamma1, gamma2);
				flux_builder.buildUsingPrimitives<0>(SIDES_WE_LEFT , fluxW);
				flux_builder.buildUsingPrimitives<0>(SIDES_WE_RIGHT, fluxE);
				flux_builder.buildUsingPrimitives<1>(SIDES_NS_DOWN , fluxS);
				flux_builder.buildUsingPrimitives<1>(SIDES_NS_UP   , fluxN);
				flux_builder.buildUsingPrimitives<2>(SIDES_BF_BACK , fluxB);
				flux_builder.buildUsingPrimitives<2>(SIDES_BF_FRONT, fluxF);
				
				for(int i=0; i<6; i++)
				{
					assert(!isnan(fluxE[i]));
					assert(!isnan(fluxW[i]));
					assert(!isnan(fluxN[i]));
					assert(!isnan(fluxS[i]));
					assert(!isnan(fluxF[i]));
					assert(!isnan(fluxB[i]));
				}
				
				Real RHS[6];
				
				for(int i=0; i<6; i++)
					RHS[i] = minuslambda*(fluxE[i] - fluxW[i] + fluxN[i] - fluxS[i] + fluxF[i] - fluxB[i]);
				
				//extra term for the levelset				
				RHS[5] -= minuslambda*(1./6)*
				(SIDES_WE_LEFT[5][1] + SIDES_WE_RIGHT[5][0] + SIDES_NS_UP[5][0] + SIDES_NS_DOWN[5][1] + SIDES_BF_FRONT[5][0] + SIDES_BF_BACK[5][1])*
				(SIDES_WE_RIGHT[1][0] - SIDES_WE_LEFT[1][1] + SIDES_NS_UP[2][0] - SIDES_NS_DOWN[2][1] + SIDES_BF_FRONT[3][0] - SIDES_BF_BACK[3][1]);
				
				
				for(int i=0; i<6; i++)
					assert(!isnan(RHS[i]));
				
				block(ix, iy, iz).dsdt.r = RHS[0];
				block(ix, iy, iz).dsdt.u = RHS[1];
				block(ix, iy, iz).dsdt.v = RHS[2];
				block(ix, iy, iz).dsdt.w = RHS[3];
				block(ix, iy, iz).dsdt.s = RHS[4];
				block(ix, iy, iz).dsdt.levelset = RHS[5];
			}
}
