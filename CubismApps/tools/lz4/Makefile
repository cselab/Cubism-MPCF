include ../../Makefile.config

#/bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc64-bgq-linux-ar

AR = ar
ARFLAGS = rcs

all: liblz4.a

liblz4.a: lz4.o
	$(AR) $(ARFLAGS) liblz4.a lz4.o 

lz4.o: lz4.c lz4.h
	$(CC) $(CPPFLAGS) -c lz4.c

clean:
	rm -f *.o *.a
