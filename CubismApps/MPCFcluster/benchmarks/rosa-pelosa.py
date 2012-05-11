import os

print "Hello rosa ...so hairy!"

#create the benchmark folder
getnodedata = os.popen('hostname ; date \'+-%y-%m-%d-%H_%M_%S\'')
folderpath = ("pelo-" + getnodedata.readline() + getnodedata.readline()).replace('\n', '')
print "folderpath is %s" % folderpath
if not os.path.exists(folderpath): os.makedirs(folderpath)

ncores = 32
numa = input("Do you want to map one MPI-process per NUMA NODE? [False, True]")

dispatcher = raw_input("Which dispatcher do you want? [overlap, omp]")
assert dispatcher == "overlap" or dispatcher == "omp"

taskspernode = 1
if numa: taskspernode = 4
threads = ncores/taskspernode

#compile the program, put it in the benchmark folder
try:
	build_dir = "../makefiles"
	cwd = "/users/petrosk/projects/CubismApps/MPCFcluster/benchmarks/" #os.getcwd()
	print "current working directory is %s" % cwd
	os.chdir(build_dir)
	
	command1 = "	make clean ; make hdf=0 vtk=0 extra=\" -march=bdver1 -mtune=bdver1 -mprefer-avx128 -ftree-vectorize \" \
	config=release numa=0  precdiv=-1  bs=32 nthreads="+ str(threads) +" fftw=0  CC=gcc -j8 "
	print 'COMMAND1 is %s' % (command1)
	os.system(command1)
	
	command2 = "mv mpcfcluster " + cwd + "/" + folderpath + "/mpcfcluster"
	print 'COMMAND2 is %s' % (command2)
	os.system(command2)
finally:
	os.chdir(cwd)
	
print "done till there"
nodesamount = [27,45,90,180,360,735,1482]
nodes2pesize = {27:(3,3,3) , 45:(5,3,3) , 90:(6,5,3) , 180:(6,6,5), 360:(9,8,5), 735:(15,7,7), 1482:(19,13,6)}

#prepare the submission script. NUMA!
bpd = 16
for nodes in nodesamount:
	pesize = nodes2pesize[nodes]

	pesizetotal = pesize[0]*pesize[1]*pesize[2]*taskspernode
	assert nodes*taskspernode == pesizetotal
	
	print "nodes is %d --> %d %d %d --> %d" % (nodes, pesize[0],pesize[1],pesize[2], pesize[0]*pesize[1]*pesize[2])
	
	submissionfile = open(folderpath + "/submit-"+str(nodes), "w")
	submissionfile.writelines("#!/bin/bash -l\n")
	submissionfile.writelines("#SBATCH --ntasks=%d\n" % pesizetotal)
	submissionfile.writelines("#SBATCH --ntasks-per-node=%d\n" % taskspernode)
	submissionfile.writelines("#SBATCH --cpus-per-task=%d\n" % threads)
	submissionfile.writelines("#SBATCH --time=00:20:00\n")
	submissionfile.writelines("#SBATCH --job-name=%s-nodes%d-numa%d-dispatcher%s\n" % (folderpath , nodes, numa, dispatcher))
	submissionfile.writelines("#SBATCH --account=s70\n")
	submissionfile.writelines("\n")
	submissionfile.writelines("mkdir -p results-%dnodes-%dbpd-numa%d\n" % (nodes, bpd, numa))
	submissionfile.writelines("cd results-%dnodes-%dbpd-numa%d\n" % (nodes, bpd, numa))
	submissionfile.writelines("export OMP_NUM_THREADS=%d\n" % threads)
	submissionfile.writelines("\n")

	mybpd = [ bpd for i in range(3)]
	mypesize = [ x for x in pesize]

	myexec = "aprun -n "+ str(pesizetotal) + " -N " + str(taskspernode) + " -d " + str(threads) 
	
	if numa: 
		myexec += " -ss -S 1 "
		mybpd[1] /= 2
		mybpd[2] /= 2
		mypesize[1] *= 2
		mypesize[2] *= 2
		
	myexec += " ../mpcfcluster -ascii 0 -awk 1 " + \
	" -bpdx " + str(mybpd[0]) + " -bpdy " + str(mybpd[1]) + " -bpdz " + str(mybpd[2]) + \
	" -bubx 0.2 -buby " + str(pesize[1]*0.5/pesize[0]) + " -bubz " + str(pesize[2]*0.5/pesize[0]) + " -cfl 0.6 -dumpperiod 10 -g1 1.4 -g2 1.667 " + \
	"  -mach 1.2 -mollfactor 2 -nsteps 101 -pb 60 -pc1 0 -pc2 0 -pp 17.2 -rad 0.05 -restart 0 -saveperiod 50000 -shockpos 0.1 -sim sb -tend 0.2 -verb 1 -vp 0 " \
	" -dispatcher " + dispatcher + " -kernels sse  -nthreads " + str(threads) + \
	" -xpesize " + str(mypesize[0]) + " -ypesize " + str(mypesize[1]) + " -zpesize " + str(mypesize[2]) + \
	" -gsync 1 -report 10 -sten 1e-4 " 
	
	submissionfile.writelines("\necho myexec is \'%s\'\n\n" % myexec)
	
	submissionfile.writelines("\n" + myexec + "\n")
	submissionfile.close()
	
#shoot the mofos
try:
	os.chdir(folderpath)
	command3 = " for S in $(ls submit-*) ; do sbatch $S ; done"
	os.system(command3)
finally: os.chdir(cwd)

