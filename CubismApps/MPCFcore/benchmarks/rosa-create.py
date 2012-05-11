#this is my shit
import create

commonflags = " omp=1 -j32 config=release  avx=1 align=32 "
tgcc = create.create_executables("", " CC=gcc" + commonflags, "gcc461", 
	["", " extra=\"-march=bdver1 -mtune=bdver1\" ", " precdiv=-1 "], "rosa-exec-gcc461")

ticc = 0
#create.create_executables(" source /opt/intel/bin/compilervars.sh intel64; ", " CC=icc " + commonflags, 
#	"icc121", ["", "precdiv=-1"], "rosa-exec-icc121")

print 'Compilation with gcc.4.6.1 took %d min %d sec.' % (int(tgcc/60), tgcc - 60*int(tgcc/60))
print 'Compilation with icc.12.1 took %d min %d sec.' % (int(ticc/60), ticc - 60*int(ticc/60))
tcompilation = tgcc + ticc
print 'Overall, compilation took %d min %d sec.' % (int(tcompilation/60), tcompilation - 60*int(tcompilation/60))
