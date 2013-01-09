SHELL := /bin/bash

CC ?= gcc
LD ?= gcc
bs ?= 16
ap ?= float
weps ?= 1e-6
precdiv ?= 0
sse ?= 1
avx ?= 0
omp ?= 0
align ?= 16
bgq ?= 0

CPPFLAGS+= $(extra)

ifeq "$(ap)" "float"
       CPPFLAGS += -D_FLOAT_PRECISION_
endif

ifeq "$(config)" "release"
	ifeq "$(CC)" "icc"
		CPPFLAGS+= -DNDEBUG -O3 -xHOST  -ip -ansi-alias -fno-fnalias -inline-level=1
		ifeq "$(ap)" "float"
			CPPFLAGS += -fp-model precise 
		else
			CPPFLAGS += -fast
		endif
	else
		ifeq "$(bgq)" "1"
			OPTFLAGS += -O3 -qarch=qp -qtune=qp -qhot -qstrict -qsource -qlist -qreport -DNDEBUG
		else
			OPTFLAGS += -O3  -fstrict-aliasing -fschedule-insns2 -falign-functions=16 -DNDEBUG -ftree-vectorize
		endif
		#CPPFLAGS += -O3 -fstrict-aliasing -fschedule-insns2 -falign-functions=16 -DNDEBUG -ftree-vectorize 
	endif
else
	CPPFLAGS += -g
endif

ifeq "$(sse)" "1"
	CPPFLAGS += -D_SSE_
	ifneq "$(CC)" "icc"
		CPPFLAGS += -msse -msse2
	endif
endif

ifeq "$(avx)" "1"
	CPPFLAGS += -D_AVX_
	ifneq "$(CC)" "icc"
		CPPFLAGS += -mavx
	endif
endif

ifeq "$(precdiv)" "1"
CPPFLAGS += -D_PREC_DIV_
endif 

ifeq "$(precdiv)" "-1"
CPPFLAGS += -D_PREC_DIV_NONE_
endif 

ifeq "$(ap)" "float"
CPPFLAGS += -D_SP_COMP_
endif 

ifeq "$(omp)" "1"
	ifneq "$(CC)" "icc"
		CPPFLAGS += -fopenmp
		LIBS += -fopenmp
	else
		ifeq "$(bgq)" "1"
			CPPFLAGS += -qsmp=omp
			OPTFLAGS += -qsmp=omp
		else
			CPPFLAGS += -fopenmp
			OPTFLAGS += -fopenmp
		endif
	endif
endif

ifneq "$(CC)" "icc"
	ifneq "$(config)" "release"
		DIEGO_AVXGCC_ISSUE = -O1
	endif
endif

CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -DWENOEPS=$(weps)
CPPFLAGS += -I../../../Cubism/source

LIBS += -lm -lstdc++

ifeq "$(bgq)" "1"
	LIBS += /gpfs/DDNgpfs1/walkup/mpi_trace/bgq/lib/libmpihpm.a \
			/bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a \
			-L/bgsys/drivers/ppcfloor/spi/lib -lSPI_upci_cnk
endif

OBJECTS = main.o Convection_CPP.o Test_Convection.o Update.o MaxSpeedOfSound.o 

ifeq "$(sse)" "1"
	ifeq "$(ap)" "double"
		OBJECTS+= Convection_SSEd.o
	else
		OBJECTS+= Convection_SSE.o
	endif
endif
ifeq "$(avx)" "1"
	OBJECTS+= Convection_AVX.o
endif

VPATH := ../source/
.DEFAULT_GOAL := mpcf-core

mpcf-core: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC) $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o mpcf-core
