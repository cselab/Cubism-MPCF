#this is my shit
import create

commonflags = " omp=1 -j8 config=release avx=1 align=32"
tgcc = create.create_executables("", " CC=gcc" + commonflags, "gcc462", 
	["", " extra=\"-march=native -mtune=native\" ", " precdiv=-1 "], "bulldozer-exec-gcc462")

ticc = create.create_executables(" source /opt/intel/composerxe/bin/compilervars.sh intel64; ", " CC=icc" + commonflags, 
	"icc121", ["", " precdiv=-1 "], "bulldozer-exec-icc121")

print 'Compilation with gcc.4.6.2 took %d min %d sec.' % (int(tgcc/60), tgcc - 60*int(tgcc/60))
print 'Compilation with icc.12.1 took %d min %d sec.' % (int(ticc/60), ticc - 60*int(ticc/60))
tcompilation = tgcc + ticc
print 'Overall, compilation took %d min %d sec.' % (int(tcompilation/60), tcompilation - 60*int(tcompilation/60))
