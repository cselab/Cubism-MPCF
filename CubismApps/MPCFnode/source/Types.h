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

#include <tbb/scalable_allocator.h>
#include <fstream>
#include "math.h"

using namespace std;

#include <Grid.h>
#include <GridMorton.h>
#include <BlockLab.h>
#include <BlockProcessing.h>
#include <Profiler.h>
#include <ArgumentParser.h>
#include <SerializerIO.h>
#ifdef _USE_VTK_
#include <SerializerIO_ImageVTK.h>
#endif
#include <SerializerIO_VP.h>
#include <Timer.h>

#include "BoundaryConditions.h"

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
		if(EPSILON!=0)
		{
			const Real x = min((Real)1, max((Real)-1, phi*(((Real)1)/EPSILON)));
			const Real val_xneg = (((Real)-0.5)*x - ((Real)1))*x + ((Real)0.5);
			const Real val_xpos = (((Real)+0.5)*x - ((Real)1))*x + ((Real)0.5);
			return (x<0 ? val_xneg : val_xpos);
		}
		else
			return (phi>0 ? 0:1);
	}
	
	static Real getGamma(const Real x)
	{
		const Real hs = heaviside(x);
		return GAMMA1*hs + GAMMA2*(((Real)1)-hs);
	}
    
    static Real getPressureCorrection(const Real x)
	{
		const Real hs = heaviside(x);
		return PC1*hs + PC2*(((Real)1)-hs);
	}
    
	static void getPostShockRatio(const Real mach, const Real gamma, Real postShock[])
	{
		const double Mpost = sqrt( (pow(mach,(Real)2.)*(gamma-1.)+2.) / (2.*gamma*pow(mach,(Real)2.)-(gamma-1.)) );
		postShock[0] = (gamma+1.)*pow(mach,(Real)2.)/( (gamma-1.)*pow(mach,(Real)2.)+2. ) ;
		postShock[2] = 1./(gamma+1.) * ( 2.*gamma*pow(mach,(Real)2.)-(gamma-1.) );
		const double preShockU = mach*sqrt(gamma);
		const double postShockU = Mpost*sqrt(gamma*postShock[2]/postShock[0]);
		postShock[1] = preShockU - postShockU;
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
	Real rho, u, v, w, energy, levelset;
	
	void clear() { rho = u = v = w = energy = levelset = 0; }
};

struct StreamerGridPointASCII
{
	void operate(const FluidElement& input, ofstream& output) const 
	{
		output << input.rho << " " << input.u << " " << input.v << " " << 
        input.w << " " << input.energy << " " << input.levelset << endl;
	}
	
	void operate(ifstream& input, FluidElement& output) const 
	{
		input >> output.rho;
		input >> output.u;
		input >> output.v;
		input >> output.w;
		input >> output.energy;
		input >> output.levelset;
		
		cout << "reading: " <<  output.rho << " " 
		<<  output.u << " " 
		<<  output.v << " " 
		<<  output.w << " " 
		<<  output.energy << " " 
		<<  output.levelset << "\n" ;
	}
};

struct StreamerGridPoint //dummy
{
	static const int channels = 6;
	void operate(const FluidElement& input, Real output[6]) const
	{
		output[0] = input.rho;
		assert(input.rho >= 0);
		output[1] = input.u/input.rho;
		output[2] = input.v/input.rho;
		output[3] = input.w/input.rho;
		output[4] = (input.energy-0.5*(input.u*input.u+input.v*input.v+input.w*input.w)/input.rho)*(Simulation_Environment::getGamma(input.levelset)-1);
		output[5] = input.levelset;		
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
		output[0] = (input.energy-0.5*(input.u*input.u+input.v*input.v+input.w*input.w)/input.rho)*(Simulation_Environment::getGamma(input.levelset)-1);
	}
};

struct StreamerLevelset
{
	static const int channels = 1;
	void operate(const FluidElement& input, Real output[1]) const
	{
		output[0] = input.levelset;
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
	
	FluidElement __attribute__((aligned(_ALIGNBYTES_))) data[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_];
	Real __attribute__((aligned(_ALIGNBYTES_))) tmp[_BLOCKSIZE_][_BLOCKSIZE_][_BLOCKSIZE_][6];
    
	void clear_data()
	{
		const int N = sizeX*sizeY*sizeZ;
		FluidElement * const e = &data[0][0][0];
		for(int i=0; i<N; ++i) e[i].clear();
	}
    
	void clear_tmp()
	{
		const int N = sizeX*sizeY*sizeZ*6;
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
		output[1] = input.u;
		output[2] = input.v;
		output[3] = input.w;
		output[4] = input.energy;
		output[5] = input.levelset;		
	}
	
	void operate(const Real output[9], const int ix, const int iy, const int iz) const
	{
		FluidElement& input = ref.data[iz][iy][ix];
		
		input.rho = output[0];
		assert(input.rho >= 0);
		input.u = output[1];
		input.v = output[2];
		input.w = output[3];
		input.energy = output[4];
		input.levelset = output[5];		
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

typedef Grid<FluidBlock, tbb::scalable_allocator> FluidGrid;
typedef BlockProcessing_TBB<FluidBlock> BlockProcessing;

