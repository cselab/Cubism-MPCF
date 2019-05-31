SHELL := /bin/bash

CC = mpic++

bs ?= 32
align ?= 16
bgq ?= 0
omp ?= 1
ap ?= float
hdf ?= 1
config ?= release
alphaeps ?= 0.0

ifeq "$(omp)" "1"
	ifeq "$(CC)" "icc"
		CPPFLAGS += -openmp
		OPTFLAGS += -openmp
	else
		ifeq "$(bgq)" "1"
			CPPFLAGS += -qsmp=speculative
			OPTFLAGS += -qsmp=speculative
		else
			CPPFLAGS += -fopenmp
			OPTFLAGS += -fopenmp
		endif
	endif
endif

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
		ifeq "$(bgq)" "1"
			OPTFLAGS += -O3 -qsmp=speculative -qhot=level=1  -qarch=qp -qtune=qp -qprefetch=aggressive
		else
			OPTFLAGS += -O3 -fno-expensive-optimizations -DNDEBUG
		endif
	endif
endif

hdf-inc ?=.
hdf-lib ?=.
fftw-inc ?=.
fftw-lib ?=.

ifneq "$(findstring brutus,$(shell hostname))" ""
	mpi-inc = ${MPI_ROOT}/include/
	mpi-lib = ${MPI_ROOT}/lib/
        hdf-inc = $(HDF5_ROOT)/include/
        hdf-lib = $(HDF5_ROOT)/lib/
endif

ifneq "$(findstring fen,$(shell hostname))" ""
       hdf-inc = ${HDF5_INC}
       hdf-lib = ${HDF5_LIB}
       fftw-inc = ${FFTW3_INC}
       fftw-lib = ${FFTW3_LIB}
endif

ifneq "$(findstring cetus,$(shell hostname))" ""
       hdf-inc = /soft/libraries/hdf5/current/cnk-xl/current/include
       hdf-lib = /soft/libraries/hdf5/current/cnk-xl/current/lib
       fftw-inc = /soft/libraries/alcf/current/xl/FFTW3/include/
       fftw-lib = /soft/libraries/alcf/current/xl/FFTW3/lib/
endif

ifneq "$(findstring mira,$(shell hostname))" ""
       hdf-inc = /soft/libraries/hdf5/current/cnk-xl/current/include
       hdf-lib = /soft/libraries/hdf5/current/cnk-xl/current/lib
       fftw-inc = /soft/libraries/alcf/current/xl/FFTW3/include/
       fftw-lib = /soft/libraries/alcf/current/xl/FFTW3/lib/
endif

ifneq "$(findstring juq,$(shell hostname))" ""
         hdf-inc = ${HDF5_INCLUDE}
         hdf-lib = ${HDF5_LIB}
         fftw-inc = ${FFTW_INCLUDE}
         fftw-lib = ${FFTW_LIB}
endif

ifneq "$(findstring panda,$(shell hostname))" ""
	hdf-inc = ${HDF5_ROOT}/include
	hdf-lib = ${HDF5_ROOT}/lib
endif
ifneq "$(findstring falcon,$(shell hostname))" ""
	hdf-inc = ${HDF5_ROOT}/include
	hdf-lib = ${HDF5_ROOT}/lib
endif

ifeq "$(hdf)" "1"
	CPPFLAGS += -I$(hdf-inc) -D_USE_HDF_
	LIBS += -L$(hdf-lib) -lhdf5

	ifneq "$(findstring juq,$(shell hostname))" ""
	      LIBS += -L/bgsys/local/szip/lib/ -lsz -L/bgsys/local/zlib/lib/ -lz
	endif

        ifneq "$(findstring fen,$(shell hostname))" ""
              #CPPFLAGS += -DH5_USE_16_API
              #LIBS += -L$SZIP_LIB -lsz -L$ZLIB_LIB -lz
              LIBS += -L/cineca/prod/libraries/szip/2.1/bgq-xl--1.0/lib -lsz -L/cineca/prod/libraries/zlib/1.2.7/bgq-gnu--4.4.6/lib -lz
              CPPFLAGS += -D_ON_FERMI_
        endif

        ifneq "$(findstring cetus,$(shell hostname))" ""
              LIBS += -L/soft/libraries/alcf/current/xl/SZIP/lib/ -L/soft/libraries/alcf/current/xl/ZLIB/lib/ -lz
        endif

        ifneq "$(findstring mira,$(shell hostname))" ""
              LIBS += -L/soft/libraries/alcf/current/xl/SZIP/lib/ -L/soft/libraries/alcf/current/xl/ZLIB/lib/ -lz
        endif
endif

CPPFLAGS += -g
CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -D_BLOCKSIZEX_=$(bs) -D_BLOCKSIZEY_=$(bs) -D_BLOCKSIZEZ_=$(bs) -DALPHAEPS=$(alphaeps)
CPPFLAGS += -I.
CPPFLAGS += -I../../../Cubism/include
CPPFLAGS += -I../../../CubismApps/MPCFcore/source
CPPFLAGS += -I../../../CubismApps/MPCFnode/source
CPPFLAGS += -I../../../CubismApps/MPCFcluster/source
CPPFLAGS += -I../../../CubismApps/tools/CubismZ/CubismApps/tools/fpzip-1.1.0/inc
LIBS += -L../../../CubismApps/tools/CubismZ/CubismApps/tools/fpzip-1.1.0/lib -lfpzip -lz -ldl

.DEFAULT_GOAL := gridmanip

OBJECTS = main.o ../../../CubismApps/MPCFnode/source/Types.o \
		  ../../../CubismApps/MPCFnode/source/Streamer.o

all: gridmanip

gridmanip: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o gridmanip *~
