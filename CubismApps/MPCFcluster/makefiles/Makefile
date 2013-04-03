SHELL := /bin/bash

CC ?= mpicxx
LD ?= mpicxx

include ../../Makefile.config

mpi-inc ?=.
mpi-lib ?=.
fftw-inc ?=.
fftw-lib ?=.

ifneq "$(findstring silvio,$(shell hostname))" ""
	mpi-inc = /home/menahel/openmpi-1.5.5/include/
	mpi-lib = /home/menahel/openmpi-1.5.5/lib/
	ifeq "$(CC)" "icc"
		fftw-inc = /home/hbabak/fftw-2.1.5_icc/include/
		fftw-lib = /home/hbabak/fftw-2.1.5_icc/lib/
	else
                fftw-inc = /home/hbabak/fftw-2.1.5_gcc/include/
                fftw-lib = /home/hbabak/fftw-2.1.5_gcc/lib/
	endif
endif

ifneq "$(findstring brutus,$(shell hostname))" ""
	mpi-inc = ${MPI_ROOT}/include/
	mpi-lib = ${MPI_ROOT}/lib/
	fftw-inc = /cluster/work/infk/hbabak/apps/fftw-2.1.5/include/
	fftw-lib = /cluster/work/infk/hbabak/apps/fftw-2.1.5/lib/
endif

ifneq "$(findstring rosa,$(shell hostname))" ""
	mpi-inc = $(CRAY_MPICH2_DIR)/include/
	mpi-lib = $(CRAY_MPICH2_DIR)/lib/
	fftw-inc = /users/petrosk/fftw-2.1.5_gcc/include/
	fftw-lib = /users/petrosk/fftw-2.1.5_gcc/lib/
endif

ifneq "$(findstring daint,$(shell hostname))" ""
        mpi-inc = $(CRAY_MPICH2_DIR)/include/
        mpi-lib = $(CRAY_MPICH2_DIR)/lib/
        fftw-inc = /users/petrosk/fftw-2.1.5_gcc/include/
        fftw-lib = /users/petrosk/fftw-2.1.5_gcc/lib/
endif

#CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -D_BLOCKSIZEX_=$(bs) -D_BLOCKSIZEY_=$(bs) -D_BLOCKSIZEZ_=$(bs) -DWENOEPS=$(weps) -DNTHREADS=$(nthreads) $(extra)
CPPFLAGS += -I../../../Cubism/source/ -I../source -I../../MPCFcore/source/ -I../../MPCFnode/source/
CPPFLAGS +=  -I$(mpi-inc)

ifneq "$(findstring rosa,$(shell hostname))" ""
	ifeq "$(CC)" "icc"
		LIBS += -L$(mpi-lib) -lmpich_intel -lmpichcxx_intel
	else
		LIBS += -L$(mpi-lib) -lmpich -lmpichcxx
	endif
else ifneq "$(findstring daint,$(shell hostname))" ""
        ifeq "$(CC)" "icc"
                LIBS += -L$(mpi-lib) -lmpich_intel -lmpichcxx_intel
        else
                LIBS += -L$(mpi-lib) -lmpich -lmpichcxx
        endif
else
	ifeq "$(findstring mpi,$(CC))" ""
		LIBS += -L$(mpi-lib) -lmpi -lmpi_cxx
    	endif

	#LIBS += -L$(mpi-lib) -lmpi -lmpi_cxx
endif

ifeq "$(fftw)"  "1"
	#FFTW always in double precision
	CPPFLAGS += -I$(fftw-inc) -D_USE_FFTW_ 
	LIBS += -L$(fftw-lib) -ldfftw -ldfftw_mpi -lfftw -lfftw_mpi
endif

##################
VPATH := ../source/ ../../MPCFcore/source/ ../../MPCFnode/source/ ../../../Cubism/source/
.DEFAULT_GOAL := mpcf-cluster

OBJECTS = main.o 
OBJECTS +=  ../../MPCFnode/makefiles/FlowStep_LSRK3.o  ../../MPCFnode/makefiles/Test_SteadyState.o ../../MPCFnode/makefiles/Test_ShockBubble.o ../../MPCFnode/makefiles/Test_SIC.o ../../MPCFnode/makefiles/Types.o ../../MPCFnode/makefiles/Test_Cloud.o ../../MPCFnode/makefiles/WaveletCompressor.o
OBJECTS += ../../MPCFcore/makefiles/Convection_CPP.o ../../MPCFcore/makefiles/Update.o ../../MPCFcore/makefiles/MaxSpeedOfSound.o 

#../../Cubism/makefiles/
OBJECTS += Profiler.o Histogram.o

ifeq "$(qpx)" "1"
          OBJECTS += ../../MPCFcore/makefiles/WenoSOA2D_QPX.o
          OBJECTS += ../../MPCFcore/makefiles/HLLESOA2D_QPX.o
          OBJECTS += ../../MPCFcore/makefiles/DivSOA2D_QPX.o
endif

ifeq "$(sse)" "1"
	ifeq "$(ap)" "double"
		OBJECTS += ../../MPCFcore/makefiles/Convection_SSEd.o
	else
		OBJECTS += ../../MPCFcore/makefiles/Convection_SSE.o
	endif
endif

ifeq "$(avx)" "1"
	OBJECTS += ../../MPCFcore/makefiles/Convection_AVX.o
endif

all: mpcf-cluster

mpcf-cluster: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

../../MPCFcore/makefiles/%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra) -c $^ -o $@

../../MPCFnode/makefiles/%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra) -c $^ -o $@
	
../../MPCFnode/makefiles/WaveletCompressor.o: WaveletCompressor.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) $(extra)  -c $^ -o $@

clean:
	rm -f *.o mpcf-cluster

cleanall:
	(cd ../../MPCFnode/makefiles; make cleanall)
	make clean;

