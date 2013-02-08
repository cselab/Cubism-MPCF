/*
 *  Types.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#pragma once

#ifdef _FLOAT_PRECISION_
typedef float Real;
#else
typedef double Real;
#endif

#include <fstream>
#include "math.h"

using namespace std;

#include <Grid.h>
#include <GridMorton.h>
#include <BlockLab.h>
//#include <BlockProcessing.h>
#include <Profiler.h>
#include <ArgumentParser.h>
#include <SerializerIO.h>
#ifdef _USE_VTK_
#include <SerializerIO_ImageVTK.h>
#endif
#include <SerializerIO_VP.h>
#include <Timer.h>

#include "BoundaryConditions.h"

#define SETUP_MARKERS_IC \
const double mix_gamma = 1 + (G2*G1)/(G1*bubble+G2*(1-bubble)); \
const double mix_pinf  = (mix_gamma-1)/mix_gamma * (F1/G1*(1-bubble) + F2/G2*bubble); \
b(ix, iy, iz).G  = 1./(mix_gamma-1); \
b(ix, iy, iz).P = mix_gamma*mix_pinf/(mix_gamma-1); \
const double ke = 0.5*(pow(b(ix, iy, iz).u,2)+pow(b(ix, iy, iz).v,2)+pow(b(ix, iy, iz).w,2))/b(ix, iy, iz).rho; \
b(ix, iy, iz).energy   = pressure*b(ix, iy, iz).G + b(ix, iy, iz).P + ke;

struct sort_pred {
    bool operator()(const std::pair<Real,Real> &left, const std::pair<Real,Real> &right) {
        return abs(left.first) < abs(right.first);
    }
};

class Simulation_Environment
{
public:
	static Real EPSILON, GAMMA1, GAMMA2;	
	static Real shock_pos, mach;
    static Real reinit_tau_factor;
    static int reinit_freq, reinit_steps;
    static Real PC1, PC2;
	
    static Real heaviside(const Real phi)
    {
        return (phi>0? 0:1);
    }
    
    static Real heaviside_smooth(const Real phi)
    {
        const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/EPSILON)));
        const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
        const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
        return (x<0 ? val_xneg : val_xpos);
    }
    
    static void getPostShockRatio(const Real pre_shock[3], const Real mach, const Real gamma, const Real pc, Real postShock[3])
    {
        const double Mpost = sqrt( (pow(mach,(Real)2.)*(gamma-1.)+2.) / (2.*gamma*pow(mach,(Real)2.)-(gamma-1.)) );
        postShock[0] = (gamma+1.)*pow(mach,(Real)2.)/( (gamma-1.)*pow(mach,(Real)2.)+2.)*pre_shock[0] ;
        postShock[2] = 1./(gamma+1.) * ( 2.*gamma*pow(mach,(Real)2.)-(gamma-1.))*pre_shock[2];
        const double preShockU = mach*sqrt(gamma*(pc+pre_shock[2])/pre_shock[0]);
        const double postShockU = Mpost*sqrt(gamma*(pc+postShock[2])/postShock[0]);
        postShock[1] = preShockU - postShockU;
    }
    
    /**
     * Computes post shock values given the post shock pressure
     * 
     * @param array holding the pre shock values
     * @param pressure of the shock
     * @param array containing the computed post shock values
     *
     * @return void
     */
    static void getPostShockRatio(const Real pre_shock[3], const Real p_shock, Real postShock[3])
    {
        postShock[0] = pre_shock[0]*pow((double)p_shock/(3.31e4)+1,(1./7.));
        postShock[1] = pre_shock[1] + sqrt((postShock[0]-pre_shock[0])/(postShock[0]*pre_shock[0])*(p_shock-pre_shock[2]));
        postShock[2] = p_shock;
    }
};

class Simulation
{
public:
    
	virtual void run() = 0;
	virtual void paint() { }
	virtual void setup() { }
};

struct FluidElement
{
    Real rho, u, v, w, energy, G, P;

    void clear() { rho = u = v = w = energy = G = P = 0; }
};

struct StreamerGridPointASCII
{
	void operate(const FluidElement& input, ofstream& output) const 
	{
        output << input.rho << " " << input.u << " " << input.v << " " << 
        input.w << " " << input.energy << " " << input.G << " " << input.P << endl;
	}
	
	void operate(ifstream& input, FluidElement& output) const 
	{
		input >> output.rho;
		input >> output.u;
		input >> output.v;
		input >> output.w;
		input >> output.energy;
		input >> output.G;
		input >> output.P;
		
		cout << "reading: " <<  output.rho << " " 
		<<  output.u << " " 
		<<  output.v << " " 
		<<  output.w << " " 
		<<  output.energy << " " 
        <<  output.G << " " 
		<<  output.P << "\n" ;
	}
};

struct StreamerGridPoint //dummy
{
    static const int channels = 7;

	void operate(const FluidElement& input, Real output[7]) const
	{
		output[0] = input.rho;
		assert(input.rho >= 0);
		output[1] = input.u/input.rho;
		output[2] = input.v/input.rho;
		output[3] = input.w/input.rho;
        output[4] = (input.energy-0.5*(input.u*input.u+input.v*input.v+input.w*input.w)/input.rho - input.P)/input.G;
		output[5] = input.G;		
  		output[6] = input.P;		
	}
};

struct StreamerDensity
{
	static const int channels = 1;

	void operate(const FluidElement& input, Real output[1]) const
	{
		output[0] = input.rho;
	}
};

struct StreamerXVelocity
{
	static const int channels = 1;

	void operate(const FluidElement& input, Real output[1]) const
	{
		output[0] = input.u/input.rho;
	}
};

struct StreamerYVelocity
{
	static const int channels = 1;

	void operate(const FluidElement& input, Real output[1]) const
	{
		output[0] = input.v/input.rho;
	}
};

struct StreamerPressure
{
	static const int channels = 1;

	void operate(const FluidElement& input, Real output[1]) const
	{
        output[0] = (input.energy-0.5*(input.u*input.u+input.v*input.v+input.w*input.w)/input.rho - input.P)/input.G;
	}
};

struct StreamerG
{
	static const int channels = 1;

	void operate(const FluidElement& input, Real output[1]) const
	{
		output[0] = input.G;
	}
};

struct FluidBlock
{
	static const int sizeX = _BLOCKSIZE_;
	static const int sizeY = _BLOCKSIZE_;
	static const int sizeZ = _BLOCKSIZE_;

	static const int gptfloats = sizeof(FluidElement)/sizeof(Real);
	
	typedef FluidElement ElementType;
	typedef FluidElement element_type;
	
	FluidElement __attribute__((__aligned__(_ALIGNBYTES_))) data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
    
	Real __attribute__((__aligned__(_ALIGNBYTES_))) tmp[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][gptfloats];
    
	void clear_data()
	{
		const int N = sizeX*sizeY*sizeZ;
		FluidElement * const e = &data[0][0][0];
		for(int i=0; i<N; ++i) e[i].clear();
	}
    
	void clear_tmp()
	{    
        const int N = sizeX * sizeY * sizeZ * gptfloats;

        Real * const e = &tmp[0][0][0][0];
        for(int i=0; i<N; ++i) e[i] = 0;
	}
    
	void clear()
	{
		clear_data();
		clear_tmp();
	}
    
	inline FluidElement& operator()(int ix, int iy=0, int iz=0)
	{
		assert(ix>=0 && ix<sizeX);
		assert(iy>=0 && iy<sizeY);
		assert(iz>=0 && iz<sizeZ);
		
		return data[iz][iy][ix];
	}
	
	template <typename Streamer>
	inline void Write(ofstream& output, Streamer streamer) const
	{
		for(int iz=0; iz<sizeZ; iz++)
			for(int iy=0; iy<sizeY; iy++)
				for(int ix=0; ix<sizeX; ix++)
					streamer.operate(data[iz][iy][ix], output);
	}
	
	template <typename Streamer>
	inline void Read(ifstream& input, Streamer streamer)
	{
		for(int iz=0; iz<sizeZ; iz++)
			for(int iy=0; iy<sizeY; iy++)
				for(int ix=0; ix<sizeX; ix++)
					streamer.operate(input, data[iz][iy][ix]);
	}
};

template <> inline void FluidBlock::Write<StreamerGridPoint>(ofstream& output, StreamerGridPoint streamer) const
{
	output.write((const char *)&data[0][0][0], sizeof(FluidElement)*sizeX*sizeY*sizeZ);
}

template <> inline void FluidBlock::Read<StreamerGridPoint>(ifstream& input, StreamerGridPoint streamer)
{
	input.read((char *)&data[0][0][0], sizeof(FluidElement)*sizeX*sizeY*sizeZ);
}

struct StreamerDummy_HDF5 
{
	static const int NCHANNELS = 9;
	
	FluidBlock& ref;
	
	StreamerDummy_HDF5(FluidBlock& b): ref(b){}
	
	void operate(const int ix, const int iy, const int iz, Real output[9]) const
	{
		const FluidElement& input = ref.data[iz][iy][ix];
		
		output[0] = input.rho;
		assert(input.rho >= 0);
		output[1] = input.u/input.rho;
		output[2] = input.v/input.rho;
		output[3] = input.w/input.rho;        
        output[4] = (input.energy-0.5*(input.u*input.u+input.v*input.v+input.w*input.w)/input.rho - input.P)/input.G;
		output[5] = input.G;		
        output[6] = input.P;		
	}
	
	void operate(const Real output[9], const int ix, const int iy, const int iz) const
	{
		FluidElement& input = ref.data[iz][iy][ix];
		
		input.rho = output[0];
		assert(input.rho >= 0);
		input.u = output[1]*output[0];
		input.v = output[2]*output[0];
		input.w = output[3]*output[0];
		input.energy = output[4]*output[5] + 0.5*output[0]*(output[1]*output[1]+output[2]*output[2]+output[3]*output[3])+output[6];
		input.G = output[5];		        
  		input.P = output[6];		
	}
	
	static const char * getAttributeName() { return "Tensor"; } 
};

struct StreamerFromTemp_HDF5
{
	static const int NCHANNELS = 3;
	
	FluidBlock& ref;
	
	StreamerFromTemp_HDF5(FluidBlock& b): ref(b) {}
	
	void operate(const int ix, const int iy, const int iz, Real output[3]) const
	{
		output[0] = ref.tmp[iz][iy][ix][0];
		output[1] = ref.tmp[iz][iy][ix][1];
		output[2] = ref.tmp[iz][iy][ix][2];
	}
	
	static const char * getAttributeName() { return "Vector"; } 
};

//typedef Grid<FluidBlock, tbb::scalable_allocator> FluidGrid;
typedef Grid<FluidBlock, std::allocator> FluidGrid;
//typedef BlockProcessing_TBB<FluidBlock> BlockProcessing;
