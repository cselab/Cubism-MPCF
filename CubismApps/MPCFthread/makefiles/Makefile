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

LIBS += /gpfs/DDNgpfs1/walkup/mpi_trace/bgq/lib/libmpihpm_smp.a \
			/bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a \
			-L/bgsys/drivers/ppcfloor/spi/lib -lSPI_upci_cnk

OBJECTS = main-weno.o 

VPATH := ../source/
.DEFAULT_GOAL := weno

all: weno

weno: $(OBJECTS)
	$(CC)  $(extra) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CC) $(extra) -c $^ -o $@

clean:
	rm -f *.o weno