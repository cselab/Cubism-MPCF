SHELL := /bin/bash

CC ?= gcc
LD ?= gcc

bs ?= 16
ap ?= float
weps ?= 1e-6
precdiv ?= 0
sse ?= 1
avx ?= 0
align ?= 16
nthreads ?= 24
omp ?= 1
hdf ?= 1
vtk ?= 1
numa ?= 0
cvt ?= 0

ifeq "$(ap)" "float"
        CPPFLAGS += -D_FLOAT_PRECISION_
endif

ifeq "$(config)" "release"
	ifeq "$(CC)" "icc"
		OPTFLAGS+= -DNDEBUG -O3 -xHOST  -ip -ansi-alias -fno-fnalias -inline-level=1
		ifeq "$(ap)" "float"
			CPPFLAGS += -fp-model precise 
		else
			CPPFLAGS += -fast
		endif
	else
		OPTFLAGS += -O3 -msse -msse2 -fstrict-aliasing -fschedule-insns2 -falign-functions=16 -DNDEBUG -ftree-vectorize
	endif
else
	ifeq "$(CC)" "icc"
		CPPFLAGS += -g 
	else
		CPPFLAGS += -g -msse -msse2
	endif
endif

ifeq "$(CC)" "icc"
	ifeq "$(avx)" "1"
		CPPFLAGS += -xAVX -D_AVX_
	endif
	CPPFLAGS += -Wbrief -wd18000 -diag-disable remark
else
	ifeq "$(avx)" "1"
		CPPFLAGS += -mavx -D_AVX_
	endif	
endif

CPPFLAGS += -Wno-deprecated 

ifneq "$(findstring rosa,$(shell hostname))" ""
	ifeq "$(CC)" "cc"
		LIBS += -dynamic
	endif
endif

ifneq "$(findstring todi,$(shell hostname))" ""
	ifeq "$(CC)" "cc"
		LIBS += -dynamic
	endif
endif

ifeq "$(cvt)" "1"
	CPPFLAGS += -D_USE_CVT_
endif

ifeq "$(precdiv)" "1"
CPPFLAGS += -D_PREC_DIV_
endif 

ifeq "$(precdiv)" "-1"
	CPPFLAGS += -D_PREC_DIV_NONE_
	ifeq "$(CC)" "icc"
		#OPTFLAGS += -ipo
	endif
endif 

ifeq "$(ap)" "float"
	CPPFLAGS += -D_SP_COMP_
endif 

ifeq "$(sse)" "1"
	CPPFLAGS += -D_SSE_
endif

ifeq "$(omp)" "1"
	ifeq "$(CC)" "icc"
		CPPFLAGS += -openmp	
		OPTFLAGS += -openmp
	else
		CPPFLAGS += -fopenmp
		OPTFLAGS += -fopenmp
	endif	
endif

vtk-inc ?=/usr/local/include/vtk-5.2/ 
vtk-lib ?=/usr/local/lib/vtk-5.2/ 
tbb-inc ?=/Users/paralleluser/tbb22_012oss/include/
tbb-lib ?=/Users/paralleluser/tbb22_012oss/ia32/cc4.0.1_os10.5.4/lib/
numa-inc ?=/cluster/work/infk/diegor/numactl-2.0.8-rc4 
numa-lib ?=/cluster/work/infk/diegor/numactl-2.0.8-rc4
hdf-inc ?=.
hdf-lib ?=.

ifeq "$(shell hostname)" "reynolds.inf.ethz.ch"
	vtk-inc =/usr/local/include/vtk-5.6/
	vtk-lib =/usr/local/lib/vtk-5.6/
	tbb-inc =/Users/babak/tbb30_196oss/include/
	tbb-lib =/Users/babak/tbb30_196oss/lib/
endif

ifeq "$(shell hostname)" "silvio"
	vtk-inc = /usr/include/vtk-5.6/
	vtk-lib = /usr/lib/vtk-5.6/
	tbb-inc = /home/hbabak/my_packages/tbb30_196oss/include/
	tbb-lib = /home/hbabak/my_packages/tbb30_196oss/build/linux_intel64_icc_cc4.4.5_libc2.11.3_kernel2.6.36_release/
	hdf-inc = /home/hbabak/hdf5-1.8.8_gcc_serial/include/
	hdf-lib = /home/hbabak/hdf5-1.8.8_gcc_serial/lib/
	numa-inc = /home/hbabak/numactl-2.0.7/include/
	numa-lib = /home/hbabak/numactl-2.0.7/
endif

ifneq "$(findstring brutus,$(shell hostname))" ""
	vtk-inc = /cluster/work/infk/cconti/VTK5.8_gcc/include/vtk-5.8/
	vtk-lib = /cluster/work/infk/cconti/VTK5.8_gcc/lib/vtk-5.8/
	tbb-inc = /cluster/apps/intel/Compiler/11.1/069/tbb/include/
	tbb-lib = /cluster/apps/intel/Compiler/11.1/069/tbb/lib/intel64/
	hdf-inc = /cluster/work/infk/hbabak/apps/hdf5-1.8.8_gcc_serial/include/
	hdf-lib = /cluster/work/infk/hbabak/apps/hdf5-1.8.8_gcc_serial/lib/
	numa-inc = /cluster/work/infk/hbabak/numactl-2.0.7/include/
	numa-lib = /cluster/work/infk/hbabak/numactl-2.0.7/
endif

ifneq "$(findstring rosa,$(shell hostname))" ""
        vtk-inc = /users/petrosk/VTK/include/vtk-5.8/
        vtk-lib = /users/petrosk/VTK/lib/vtk-5.8/
        tbb-inc = /apps/rosa/intel/composerxe-2011.3.174/tbb/include/
        tbb-lib = /apps/rosa/intel/composerxe-2011.3.174/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/
	ifeq "$(CC)" "icc"
		hdf-inc = /users/petrosk/hdf5-1.8.8_icc_serial/include/
		hdf-lib = /users/petrosk/hdf5-1.8.8_icc_serial/lib/
	else
		hdf-inc = /users/petrosk/hdf5-1.8.8_gcc_serial/include/
		hdf-lib = /users/petrosk/hdf5-1.8.8_gcc_serial/lib/
	endif
	numa-inc = /users/petrosk/numactl-2.0.7/
	numa-lib = /users/petrosk/numactl-2.0.7/lib64/
endif

ifneq "$(findstring todi,$(shell hostname))" ""
        vtk-inc = /users/petrosk/VTK/include/vtk-5.8/
        vtk-lib = /users/petrosk/VTK/lib/vtk-5.8/
        tbb-inc = /apps/rosa/intel/composerxe-2011.3.174/tbb/include/
        tbb-lib = /apps/rosa/intel/composerxe-2011.3.174/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/
	ifeq "$(CC)" "icc"
		hdf-inc = /users/petrosk/hdf5-1.8.8_icc_serial/include/
		hdf-lib = /users/petrosk/hdf5-1.8.8_icc_serial/lib/
	else
		hdf-inc = /users/petrosk/hdf5-1.8.8_gcc_serial/include/
		hdf-lib = /users/petrosk/hdf5-1.8.8_gcc_serial/lib/
	endif
	numa-inc = /users/petrosk/numactl-2.0.7/
	numa-lib = /users/petrosk/numactl-2.0.7/lib64/
endif

CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -D_BLOCKSIZEX_=$(bs) -D_BLOCKSIZEY_=$(bs) -D_BLOCKSIZEZ_=$(bs) -DWENOEPS=$(weps) -DNTHREADS=$(nthreads)
CPPFLAGS += -I../../../Cubism/source/ -I../source -I../../MPCFcore/source/
CPPFLAGS += -I$(vtk-inc) -I$(tbb-inc)


ifeq "$(hdf)" "1"
	CPPFLAGS += -I$(hdf-inc) -D_USE_HDF_
	LIBS += -L$(hdf-lib) -lhdf5
endif

CPPFLAGS += -I$(tbb-inc) -I$(hdf-inc)

LIBS += -lstdc++

LIBS += -L$(tbb-lib) \
	-ltbb \
	-ltbbmalloc

ifeq "$(vtk)" "1"
	CPPFLAGS += -I$(vtk-inc) -D_USE_VTK_

	LIBS += -L$(vtk-lib) \
		-lvtkVolumeRendering \
		-lvtkRendering \
		-lvtkIO \
		-lvtkGenericFiltering \
		-lvtkGraphics \
		-lvtkImaging \
		-lvtkFiltering \
		-lvtkCommon \
		-lvtkftgl \
		-lvtkDICOMParser \
		-lvtksys \
		-lvtkexoIIc \
		-lvtkNetCDF \
		-lvtkNetCDF_cxx \
		-lvtkmetaio \
		-lvtksqlite \
		-lvtkverdict \
		-lvtkexpat \
		-lvtktiff \
		-lvtkpng \
		-lvtkjpeg \
		-lvtkzlib \
		-lvtkHybrid \
		-lvtkfreetype

	ifneq "$(findstring brutus,$(shell hostname))" ""
		LIBS+= -lvtkParallel -lVPIC -lCosmo
	endif
endif

ifeq "$(numa)" "1"
	CPPFLAGS += -D_USE_NUMA_ -I$(numa-inc)
	LIBS += -L$(numa-lib) -lnuma
endif

VPATH := ../source/ ../../MPCFcore/source/ ../../../Cubism/source/
.DEFAULT_GOAL := mpcf-node

OBJECTS = main.o FlowStep_LSRK3.o Test_SteadyState.o Test_ShockTube.o Test_ShockBubble.o Test_CVT.o Convection_CPP.o Update.o MaxSpeedOfSound.o Profiler.o Types.o SurfaceTension_CPP.o DivTensor_CPP.o Diffusion_CPP.o Test_TG.o

ifeq "$(sse)" "1"
	ifeq "$(ap)" "double"
		OBJECTS += Convection_SSEd.o
	else
		OBJECTS += Convection_SSE.o
	endif
endif

ifeq "$(avx)" "1"
	OBJECTS += Convection_AVX.o
endif

#on brutus flowstep performs better without native arch
ifeq "$(findstring brutus,$(shell hostname))" ""
	FSDIEGOEXTRA=$(extra)
endif


mpcf-node: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

FlowStep_SSE_diego.o: FlowStep_SSE_diego.cpp
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $(FSDIEGOEXTRA) -c $^ -o $@

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra) -c $^ -o $@

clean:
	rm -f *.o mpcf-node
