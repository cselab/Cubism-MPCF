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
		CPPFLAGS += -O3 -fstrict-aliasing -fschedule-insns2 -falign-functions=16 -DNDEBUG -ftree-vectorize 
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
		CPPFLAGS += -openmp
		LIBS += -openmp
	endif
endif

ifneq "$(CC)" "icc"
	ifneq "$(config)" "release"
		DIEGO_AVXGCC_ISSUE = -O1
	endif
endif

CPPFLAGS += -D_ALIGNBYTES_=$(align) -D_BLOCKSIZE_=$(bs) -DWENOEPS=$(weps)
CPPFLAGS += -I../../../Cubism/v2/source

LIBS += -lm -lstdc++

OBJECTS = main.o FlowStep_CPP.o MaxSpeedOfSound.o FlowStep_Test.o Update.o SurfaceTension_CPP.o DivTensor_CPP.o Diffusion_CPP.o SurfaceTension_Test.o Diffusion_Test.o

ifeq "$(sse)" "1"
	OBJECTS+= FlowStep_SSE_diego.o
endif
ifeq "$(avx)" "1"
	OBJECTS+= FlowStep_AVX_diego.o
endif

VPATH := ../source/
.DEFAULT_GOAL := mpcfkernels

mpcfkernels: $(OBJECTS)
	$(CC) $^ -o $@ $(LIBS)

FlowStep_AVX_diego.o: FlowStep_AVX_diego.cpp
	$(CC)  $(CPPFLAGS) $(DIEGO_AVXGCC_ISSUE) -c $^ -o $@

%.o: %.cpp
	$(CC)  $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o mpcfkernels
