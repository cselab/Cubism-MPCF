/*
 *  main.cpp
 *  MPCFnode
 *
 *  Created by Diego Rossinelli on 6/14/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include <iostream>

#ifdef _SSE_
#include <xmmintrin.h>
#endif

#ifdef _USE_NUMA_
#include <numa.h>
#include <omp.h>
#endif

#include "Test_SteadyState.h"
#include "Test_ShockBubble.h"
#include "Test_SIC.h"

using namespace std;

Simulation * sim = NULL;

#ifdef _GLUT_VIZ 
struct VisualSupport
{	
	static void display(){}
	
	static void idle(void)
	{
		glClear(GL_COLOR_BUFFER_BIT);	
		sim->run();
		glutSwapBuffers();
	}
	
	static void run(int argc, const char ** argv)
	{
		static bool bSetup = false;
		
		if (!bSetup)
		{
			setup(argc, argv);
			bSetup = true;
		}
		
		glutDisplayFunc(display);
		glutIdleFunc(idle);
		
		glutMainLoop();
	}
	
	static void setup(int argc,  const char ** argv)
	{
		glutInit(&argc, const_cast<char **>(argv));
		glutInitWindowSize(800,800);
		glutInitWindowPosition(0, 0);
		glutInitDisplayMode(GLUT_DEPTH | GLUT_STENCIL |GLUT_RGBA | GLUT_DOUBLE );
		
		glutCreateWindow("MPCF-node");
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		
		glOrtho(0.0, 1.0, 0.0, 1.0, -1, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glEnable(GL_TEXTURE_2D);
		
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
		
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	}
};
#endif

int main (int argc, const char ** argv) 
{
#ifdef _USE_NUMA_
	if (numa_available() < 0)
		printf("WARNING: The system does not support NUMA API!\n");
	else
		printf("NUMA API supported!\n");
#endif
    
	ArgumentParser parser(argc, argv);	
	const bool bFlush2Zero = parser("-f2z").asBool(true);
	
#ifdef _SSE_
	if (bFlush2Zero)
#pragma omp parallel
	{
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	}
#endif
	
	
	if( parser("-sim").asString() == "steady" )
		sim = new Test_SteadyState(argc, argv);
	else if( parser("-sim").asString() == "sb" )
		sim = new Test_ShockBubble(argc, argv);
    else if( parser("-sim").asString() == "sic" )
        sim = new Test_SIC(argc, argv);
    else
	{
		printf("Study case not defined!\n"); 
		abort();
	}
	
	sim->setup();
	
	double wallclock;
	
	{
		Timer timer;
        
		timer.start();		
#ifdef _MRAG_GLUT_VIZ 
		VisualSupport::run(argc, argv);
#else
		sim->run();
#endif
		
		wallclock = timer.stop();
	}
	
	printf("we spent: %2.2f \n", wallclock);
	
    return 0;
}
