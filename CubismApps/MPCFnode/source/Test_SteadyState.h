/*
 *  Test_SteadyState.h
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "math.h"
#include <limits>

#include "Types.h"
#include "FlowStep_LSRK3.h"

#include <Profiler.h>

class Test_SteadyState: public Simulation
{
protected:
	//"constants" of the sim
    int BPDX, BPDY, BPDZ, DUMPPERIOD, SAVEPERIOD, RAMP, VERBOSITY, REPORT_FREQ, MOLLFACTOR, ANALYSISPERIOD;
	Real CFL, TEND;
	bool bRESTART, bASCIIFILES, bVP;
	
	// parameters required for testing/benchmarking
	int NSTEPS;
	bool bAWK;
	
	//state of the sim
	double t;
	int step_id;
    
    FluidGrid * grid;
	FlowStep_LSRK3 * stepper;
    
	//helpers
	ArgumentParser parser;
	Profiler profiler;
	
	void _dump(string filename);
	
	void _setup_constants();
	void _ic(FluidGrid& grid);
    
	virtual void _restart();
	virtual void _save();
	virtual void _vp(FluidGrid& grid);
	
	Real _initial_dt(int nsteps);
	void _tnext(double &tnext, double& tnext_dump, double& tend);
	
	void _refine(bool bUseIC);
	void _compress(bool bUseIC);    
public:
	
	void run();
	void paint();
	void setup();
	void verbosity(int v) { VERBOSITY = v; }
	Real get_time() const { return t; }
	int get_stepid() const { return step_id; }
    
	Test_SteadyState(const int argc, const char ** argv);
	virtual ~Test_SteadyState() { }
};

