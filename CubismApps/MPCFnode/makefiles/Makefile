include ../../Makefile.config

#########

VPATH := ../source/ ../../MPCFcore/source/ ../../../Cubism/source/
.DEFAULT_GOAL := mpcf-node

OBJECTS = main.o FlowStep_LSRK3.o Test_SteadyState.o Test_ShockBubble.o  Types.o Test_SIC.o WaveletCompressor.o
OBJECTS += ../../MPCFcore/makefiles/Convection_CPP.o ../../MPCFcore/makefiles/Update.o ../../MPCFcore/makefiles/MaxSpeedOfSound.o 

ifeq "$(qpx)" "1"
	OBJECTS += ../../MPCFcore/makefiles/WenoSOA2D_QPX.o
	OBJECTS += ../../MPCFcore/makefiles/HLLESOA2D_QPX.o
	OBJECTS += ../../MPCFcore/makefiles/DivSOA2D_QPX.o
endif

OBJECTS += Profiler.o 

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

#on brutus flowstep performs better without native arch
ifeq "$(findstring brutus,$(shell hostname))" ""
	FSDIEGOEXTRA=$(extra)
endif

all: mpcf-node

mpcf-node: $(OBJECTS)
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $^ -o $@ $(LIBS)

FlowStep_SSE_diego.o: FlowStep_SSE_diego.cpp
	$(CC) $(OPTFLAGS) $(CPPFLAGS) $(FSDIEGOEXTRA) -c $^ -o $@

WaveletCompressor.o: WaveletCompressor.cpp
	$(CC) $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@
#	$(CC) $(OPTFLAGS) $(CPPFLAGS) --std=c++11 -c $^ -o $@

%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

../../MPCFcore/makefiles/%.o: %.cpp
	$(CC)  $(OPTFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o mpcf-node

cleanall:
	(cd ../../MPCFcore/makefiles; make clean)
	make clean;
