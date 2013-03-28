hdf = 0
vtk = 0
numa = 0
cvt = 0
##nthreads ?= 24

include ../../Makefile.config
#include ../Makefile.config.peh


#ifneq "$(CC)" "icc"
#	ifneq "$(config)" "release"
#		DIEGO_AVXGCC_ISSUE = -O1
#	endif
#endif
#CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -DWENOEPS=$(weps)
#CPPFLAGS += -I../../../Cubism/source

#CPPFLAGS += -D_L1PAR
CPPFLAGS += -D_L2PAR

###############

OBJECTS = main.o Convection_CPP.o Test_Convection.o Update.o MaxSpeedOfSound.o

ifeq "$(qpx)" "1"
OBJECTS += WenoSOA2D_QPX.o 
OBJECTS += HLLESOA2D_QPX.o
OBJECTS += DivSOA2D_QPX.o
NASTYFLAGS = -qlist -qsource -qreport -O5 -qarch=qp -qtune=qp  -qhot=level=2 -qunroll
endif

OBJECTS += Convection_CPP_omp.o

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
header_files = $(wildcard ../source/*.h)

all: mpcf-core

mpcf-core: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $^ -o $@ $(LIBS)

%.o: %.cpp $(header_files)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) -c -o $@ $<

WenoSOA2D_QPX.o: WenoSOA2D_QPX.cpp WenoSOA2D_QPX.h common.h ../../MPCFthread/source/Weno_CPP.h ../../MPCFthread/source/Weno_QPX.h 
	$(CC) $(NASTYFLAGS) $(CPPFLAGS) -c -o $@ $<

HLLESOA2D_QPX.o: HLLESOA2D_QPX.cpp HLLESOA2D_QPX.h common.h
	$(CC) $(NASTYFLAGS) $(CPPFLAGS) -c -o $@ $<

DivSOA2D_QPX.o: DivSOA2D_QPX.cpp DivSOA2D_QPX.h common.h
	$(CC) $(NASTYFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o mpcf-core
