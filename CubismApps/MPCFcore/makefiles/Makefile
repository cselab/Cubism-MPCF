hdf = 0
vtk = 0
numa = 0
cvt = 0
##nthreads ?= 24

include ../../Makefile.config


#ifneq "$(CC)" "icc"
#	ifneq "$(config)" "release"
#		DIEGO_AVXGCC_ISSUE = -O1
#	endif
#endif
#CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -DWENOEPS=$(weps)
#CPPFLAGS += -I../../../Cubism/source

###############

OBJECTS = main.o Convection_CPP.o Test_Convection.o Update.o MaxSpeedOfSound.o DivTensor_CPP.o Diffusion_CPP.o

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

all: mpcf-core

mpcf-core: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC) $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o mpcf-core
