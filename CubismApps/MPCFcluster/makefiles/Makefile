SHELL := /bin/bash

CC ?= mpicxx
LD ?= mpicxx

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
fftw ?= 1
vtk ?= 1
numa ?= 0

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
		OPTFLAGS += -O3  -fstrict-aliasing -fschedule-insns2 -falign-functions=16 -DNDEBUG -ftree-vectorize
	endif
else
	ifeq "$(CC)" "icc"
		CPPFLAGS += -g 
	else
		CPPFLAGS += -g
	endif
endif

ifeq "$(CC)" "icc"
	CPPFLAGS += -Wbrief -wd18000 -diag-disable remark	
endif

ifeq "$(avx)" "1"
	CPPFLAGS += -mavx -D_AVX_
endif

CPPFLAGS += -Wno-deprecated 

ifneq "$(findstring rosa,$(shell hostname))" ""
	ifeq "$(CC)" "cc"
		LIBS += -dynamic
	endif
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
	CPPFLAGS += -D_SSE_  -msse -msse2
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
mpi-inc ?=.
mpi-lib ?=.
hdf-inc ?=.
hdf-lib ?=.
fftw-inc ?=.
fftw-lib ?=.
hdf-inc =.
hdf-lib =.

ifeq "$(shell hostname)" "reynolds.inf.ethz.ch"
	vtk-inc =/usr/local/include/vtk-5.6/
	vtk-lib =/usr/local/lib/vtk-5.6/
	tbb-inc =/Users/babak/tbb30_196oss/include/
	tbb-lib =/Users/babak/tbb30_196oss/lib/
endif

ifneq "$(findstring silvio,$(shell hostname))" ""
        vtk-inc = /usr/include/vtk-5.6/
        vtk-lib = /usr/lib/vtk-5.6/
        tbb-inc = /home/hbabak/my_packages/tbb30_196oss/include/
	ifeq "$(CC)" "icc"
		tbb-lib = /home/hbabak/tbb30_196oss/build/linux_intel64_icc_cc4.6.1_libc2.13_kernel3.0.0_release/
	else
		tbb-lib = /home/hbabak/tbb30_196oss/build/linux_intel64_gcc_cc4.6.1_libc2.13_kernel3.0.0_release/
	endif
	mpi-inc = /home/menahel/openmpi-1.5.5/include/
	mpi-lib = /home/menahel/openmpi-1.5.5/lib/
	ifeq "$(CC)" "icc"
		hdf-inc = /home/hbabak/hdf5-1.8.8/hdf5/include/
		hdf-lib = /home/hbabak/hdf5-1.8.8/hdf5/lib/
		fftw-inc = /home/hbabak/fftw-2.1.5_icc/include/
		fftw-lib = /home/hbabak/fftw-2.1.5_icc/lib/
	else
                hdf-inc = /home/hbabak/hdf5-1.8.8_gcc/include/
                hdf-lib = /home/hbabak/hdf5-1.8.8_gcc/lib/
                fftw-inc = /home/hbabak/fftw-2.1.5_gcc/include/
                fftw-lib = /home/hbabak/fftw-2.1.5_gcc/lib/
	endif
	numa-inc = /home/hbabak/numactl-2.0.7/include/
	numa-lib = /home/hbabak/numactl-2.0.7/
endif

ifneq "$(findstring brutus,$(shell hostname))" ""
        vtk-inc = /cluster/work/infk/cconti/VTK5.8_gcc/include/vtk-5.8/
        vtk-lib = /cluster/work/infk/cconti/VTK5.8_gcc/lib/vtk-5.8/
        tbb-inc = /cluster/apps/intel/Compiler/11.1/069/tbb/include/
        tbb-lib = /cluster/apps/intel/Compiler/11.1/069/tbb/lib/intel64/
        mpi-inc = ${MPI_ROOT}/include/
        mpi-lib = ${MPI_ROOT}/lib/
        ifeq "$(CC)" "icc"
                hdf-inc = /cluster/work/infk/hbabak/apps/hdf5-1.8.8/include/
                hdf-lib = /cluster/work/infk/hbabak/apps/hdf5-1.8.8/lib/
        else
                hdf-inc = /cluster/work/infk/hbabak/apps/hdf5-1.8.8_gcc/include/
                hdf-lib = /cluster/work/infk/hbabak/apps/hdf5-1.8.8_gcc/lib/
        endif
	fftw-inc = /cluster/work/infk/hbabak/apps/fftw-2.1.5/include/
	fftw-lib = /cluster/work/infk/hbabak/apps/fftw-2.1.5/lib/
	numa-inc = /cluster/work/infk/hbabak/numactl-2.0.7/include/
	numa-lib = /cluster/work/infk/hbabak/numactl-2.0.7/
endif

ifneq "$(findstring rosa,$(shell hostname))" ""
	vtk-inc = /users/petrosk/VTK/include/vtk-5.8/
        vtk-lib = /users/petrosk/VTK/lib/vtk-5.8/
        tbb-inc = /apps/rosa/intel/composerxe-2011.3.174/tbb/include/
        tbb-lib = /apps/rosa/intel/composerxe-2011.3.174/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/
	mpi-inc = $(CRAY_MPICH2_DIR)/include/
	mpi-lib = $(CRAY_MPICH2_DIR)/lib/
	ifeq "$(CC)" "icc"	
		hdf-inc = /users/petrosk/hdf5-1.8.8/include/
		hdf-lib = /users/petrosk/hdf5-1.8.8/lib/
	else
                hdf-inc = /users/petrosk/hdf5-1.8.8_gcc/include/
                hdf-lib = /users/petrosk/hdf5-1.8.8_gcc/lib/
	endif
	fftw-inc = /users/petrosk/fftw-2.1.5_gcc/include/
	fftw-lib = /users/petrosk/fftw-2.1.5_gcc/lib/
	numa-inc = /users/petrosk/numactl-2.0.7/
	numa-lib = /users/petrosk/numactl-2.0.7/lib64/
endif

CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -D_BLOCKSIZEX_=$(bs) -D_BLOCKSIZEY_=$(bs) -D_BLOCKSIZEZ_=$(bs) -DWENOEPS=$(weps) -DNTHREADS=$(nthreads) $(extra)
CPPFLAGS += -I../../../Cubism/source/ -I../source -I../../MPCFcore/source/ -I../../MPCFnode/source/
CPPFLAGS += -I$(tbb-inc) -I$(mpi-inc)

LIBS += -lstdc++ -lz -lm -lrt

LIBS += -L$(tbb-lib) \
	-ltbb \
	-ltbbmalloc

ifneq "$(findstring rosa,$(shell hostname))" ""
	ifeq "$(CC)" "icc"
		LIBS += -L$(mpi-lib) -lmpich_intel -lmpichcxx_intel
	else
		LIBS += -L$(mpi-lib) -lmpich -lmpichcxx
	endif
else
	LIBS += -L$(mpi-lib) -lmpi -lmpi_cxx
endif

ifeq "$(hdf)"  "1"
	CPPFLAGS += -I$(hdf-inc) -D_USE_HDF_
	LIBS += -L$(hdf-lib) -lhdf5
endif

ifeq "$(fftw)"  "1"
	#FFTW always in double precision
	CPPFLAGS += -I$(fftw-inc) -D_USE_FFTW_ 
	LIBS += -L$(fftw-lib) -ldfftw -ldfftw_mpi -lfftw -lfftw_mpi
endif

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
        -lvtksys
	
	ifneq "$(findstring rosa,$(shell hostname))" ""
        	LIBS += -lvtkexoIIc -lvtkNetCDF -lvtkNetCDF_cxx -lvtkmetaio -lvtksqlite -lvtkverdict -lvtkexpat -lvtktiff -lvtkpng -lvtkjpeg -lvtkzlib -lvtkHybrid -lvtkfreetype
	endif
endif

ifeq "$(numa)" "1"
	CPPFLAGS += -D_USE_NUMA_ -I$(numa-inc)
	LIBS += -L$(numa-lib) -lnuma
endif

VPATH := ../source/ ../../MPCFcore/source/ ../../MPCFnode/source/ ../../../Cubism/source/
.DEFAULT_GOAL := mpcf-cluster

OBJECTS = main.o FlowStep_LSRK3.o Test_SteadyState.o Convection_CPP.o Update.o MaxSpeedOfSound.o Profiler.o Types.o SurfaceTension_CPP.o DivTensor_CPP.o Diffusion_CPP.o Histogram.o
OBJECTS += Test_ShockBubble.o Test_SIC.o

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

mpcf-cluster: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o mpcf-cluster
