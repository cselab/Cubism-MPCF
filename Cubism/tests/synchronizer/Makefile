SHELL := /bin/bash

CC = mpic++

precision ?= single
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


.PHONY: all clean

$(warning using precision=$(precision))
$(warning using nonuniform=$(nonuniform))
ifeq ("$(precision)","double")
	CFLAGS += -DCUBISM_TEST_DOUBLE
endif

ifeq ("$(nonuniform)","yes")
	CFLAGS += -DCUBISM_TEST_NONUNIFORM
endif

all: testSynchronizerMPI.cpp
	$(CC) $(CFLAGS) -o testSynchronizerMPI testSynchronizerMPI.cpp

clean:
	rm -f testSynchronizerMPI *.o
