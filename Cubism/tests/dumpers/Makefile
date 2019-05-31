SHELL := /bin/bash

CC = mpic++

precision ?= single
hdf5precision ?= single
nonuniform ?= no

extra ?=

CFLAGS = -fopenmp -g -O2 -I../../include
CFLAGS+= $(extra)
# CFLAGS += -I${HDF5_ROOT}/include
CFLAGS+= -Wall -Wextra -Wfloat-equal -Wundef -Wcast-align
CFLAGS+= -Wwrite-strings -Wmissing-declarations -Wredundant-decls
# CFLAGS+= -Wlogical-op
CFLAGS+= -Wshadow -Woverloaded-virtual -Wuninitialized
CFLAGS+= -Wpedantic -Wno-unused-parameter # does not make much sense with function overloading


LIBS =
ifneq (${HDF5_ROOT},)
	LIBS += -L${HDF5_ROOT}/lib
endif

.PHONY: all clean

$(warning using precision=$(precision))
$(warning using hdf5precision=$(hdf5precision))
$(warning using nonuniform=$(nonuniform))
ifeq ("$(precision)","double")
	CFLAGS += -DCUBISM_TEST_DOUBLE
endif

ifeq ("$(hdf5precision)","double")
	CFLAGS += -DCUBISM_TEST_HDF5_DOUBLE_PRECISION
endif

ifeq ("$(nonuniform)","yes")
	CFLAGS += -DCUBISM_TEST__NONUNIFORM
endif

all: testDumpsMPI.cpp
	$(CC) $(CFLAGS) -o testDumpsMPI testDumpsMPI.cpp $(LIBS) -lhdf5

clean:
	rm -f testDumpsMPI *.o *.h5 *.xmf *.bin *.zbin
