#this is my shit
import create

commonflags = " omp=1 -j8 config=release precdiv=-1 vtk=0 "
tgcc = create.create_executables("source /etc/profile.d/modules.sh; module purge ; module load gcc/4.6.1; ", 
		" CC=gcc" + commonflags, "gcc461", [" extra=\"-march=native -mtune=native\" ", " numa=1 "], "brutus-exec-gcc461")
ticc = create.create_executables("source /etc/profile.d/modules.sh; module purge ; module load intel/12; ", 
		" CC=icc" + commonflags, "icc12", ["", " numa=1  "], "brutus-exec-icc12")

print 'Compilation with gcc.4.6.1 took %d min %d sec.' % (int(tgcc/60), tgcc - 60*int(tgcc/60))
print 'Compilation with icc/12 took %d min %d sec.' % (int(ticc/60), ticc - 60*int(ticc/60))
tcompilation = tgcc + ticc
print 'Overall, compilation took %d min %d sec.' % (int(tcompilation/60), tcompilation - 60*int(tcompilation/60))