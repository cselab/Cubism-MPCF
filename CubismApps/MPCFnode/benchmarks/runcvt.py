import sys
import os
import fnmatch
import datetime
import math

benches = []
#all the kernels names
#export OMP_NUM_THREADS=$NTHREADS ; ./mpcf-node  -dumpfreq 100 -savefreq 100 -ascii 0 -awk  
#-bpdx $2 -bpdy $3 -bpdz $3 -bubx 0.2 -buby 0.5 -bubz 0 -cfl 0.6 
#-dispatcher $DISPATCHER -dumpperiod 100 -g1 1.4 -g2 1.667 -gr 0 -kernels sse 
#-mach 1.2 -mollfactor 2 -nsteps 51 -nthreads $NTHREADS -pc1 0 -pc2 0 -rad 0.05 
#-reinit 0 -reinitfreq 100 -report 10 -restart 0 -saveperiod 100 -shockpos 0.1 
#-sim sb -sten 0 -sumrho 3 -tend 2.0 -verb 1 -vp 0 -peano 1
dispatchers = [" -dispatcher omp", " -dispatcher tbblight"]
peanos  = [ " -peano 0 ", " -peano 1 "]
stdkernels = [" -kernels sse ", " -kernels cpp "]

def _runall(niters, args, prerun, execfolder, outputfolder, nthreads, threads_per_node):
	start = datetime.datetime.now()
	
	for bench in benches:
		kernels = stdkernels
		
		if bench[3]:
			kernels += [ " -kernels avx " ]
			
		for disp in dispatchers:
			for pea in peanos:
				for ker in kernels:
					bs = bench[1];
					bpd = str(512/bs)
					
					newargs = disp + pea + " -nthreads " + str(nthreads) + ker + " -bpdx " + bpd + " -bpdy " + bpd +  " -bpdz " + bpd + " -nsteps " + str(niters)
					allargs =  args + newargs
					outputname = ("output bs" + str(bench[1]) + " e " + str(bench[2]) + " avx " + str(bench[3]) + " " + newargs).replace(' ', '-')
					print 'Evaluating ' + bench[0] + ": " + allargs + '...'
					maxnode = (nthreads-1)/threads_per_node;
					command1 = prerun + " export OMP_NUM_THREADS=" + str(nthreads) + " ; numactl --cpunodebind=0-"+str(maxnode) + " " + execfolder + bench[0] + allargs +  " > " + outputfolder + "/" + outputname + " 2> errors.log"
					print "COMMAND IS %s" % command1
					os.system(command1)
	
	finish = datetime.datetime.now()
	
	return (finish - start).seconds
	

def run(runargs, tdesired_seconds, core_pp_gfs, core_pb_gbs, prerun, execfolder, outputfolder, execprefix, nthreads_list, threads_per_node):
	print '======================================================================================='
	print '==========================  RUN.PY: BENCHMARKING MPCFNODE  ============================'
	print '======================================================================================='
	
	#create finalfolder if necessary
	if not os.path.exists(outputfolder):
		os.makedirs(outputfolder)
		
	print 'Fetching all executables from %s...' % (execfolder)

	files=os.listdir(execfolder)
	for fn in files:
		if fnmatch.fnmatch(fn, '*bench*'):
			info = fn.split('-')
			print "info is ", info
			bs = int(info[1].split('b')[1])
			extra = int(info[2].split('e')[1])
			
			avx = False
			if len(info) >= 4 and info[3] == "avx": 
				avx = True
				
			benches.append((fn, bs, extra, avx))

	print 'Printing all found executables...'
	for X in benches:
		print('(name, blocksize, extra flags, avx)', X)

	print 'Evaluating all benchmarks...'	

	args = runargs + " -pp " + core_pp_gfs + " -pb " + core_pb_gbs #+ " -nblocks " + nblocks 
	tdesired_seconds_per_bunch = tdesired_seconds/len(nthreads_list)
	
	for nthreads in nthreads_list:
	
		time1 = _runall(1, args, prerun, execfolder, outputfolder, nthreads, threads_per_node)
		print "With one iteration it took " + str(time1) + " seconds."

		if time1 > tdesired_seconds_per_bunch/4.:
			print "Alright we have spent enough time. Benchmarking done"
		else:
			time2 = _runall(10, args, prerun, execfolder, outputfolder, nthreads, threads_per_node)
			print "With two iteration it took " + str(time2) + " seconds."

			#estimate the number of iterations to cover the desired time
			a = max(1e-6, (time2 - time1)/9)
			b = max(0, time1 - a)
			td = tdesired_seconds_per_bunch - time1 - time2 - 30 #safety margin
			niters = int(max(1, min(1e3, math.floor((td - b)/a))))
			print "a=%.2f [s/iteration], b=%.2f [s], td=%.2f [s], niters=%d [iteration]\n" % (a, b, td, niters)
			
			if niters > 10: 
				timefinal = _runall(niters, args, prerun, execfolder, outputfolder, nthreads, threads_per_node)
				print "In total we have spent " + str(time1 + time2 + timefinal) + " seconds."
			else:
				print "Alright we have spent enough time. Benchmarking done"
