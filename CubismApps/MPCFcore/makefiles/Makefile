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

ifeq "$(bgq)" "1"
OBJECTS += WenoSOA2D_QPX.o 
OBJECTS += HLLESOA2D_QPX.o
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

WenoSOA2D_QPX.o:	WenoSOA2D_QPX.cpp $(header_files)
	$(CC) -qreport -O5 -qarch=qp -qtune=qp  -qhot=level=2 -qunroll $(CPPFLAGS) -c -o $@ $<

HLLESOA2D_QPX.o:        HLLESOA2D_QPX.cpp $(header_files)
	$(CC) -qreport -O5 -qarch=qp -qtune=qp  -qhot=level=2 -qunroll $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f *.o mpcf-core
