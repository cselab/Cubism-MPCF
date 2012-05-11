import sys
import os
import fnmatch
import datetime
import math

benches = []
#all the kernels names
ssekernels = [ "FS_SSE_diego", "Update_SSE", "MaxSOS_SSE", "SurfaceTension_SSE", "Diffusion_SSE", "Diffusion_CPP", "SurfaceTension_CPP", "FS_CPP"]
avxkernels = [ "FS_AVX_diego", "Diffusion_AVX", "SurfaceTension_AVX"] 
	
def _runall(niters, args, prerun, execfolder, outputfolder, footprint_megabytes):
	start = datetime.datetime.now()
		
	for bench in benches:
		kernels = ssekernels
		if bench[3]:
			kernels += avxkernels
		
		bs = bench[1];
		footprint_per_block_megabytes = bs**3 * 12. * 4 / 1024. / 1024.
		print "footprint_per_block_megabytes is %.2f MB" % footprint_per_block_megabytes
		nblocks = max(1, math.floor(footprint_megabytes / footprint_per_block_megabytes))
			
		for k in kernels:
			allargs =  args + " -n " + str(niters) + " -kernel " + k +  " -nblocks " + str(nblocks) 
			outputname = ("output bs" + str(bench[1]) + " e " + str(bench[2]) + " avx " + str(bench[3]) + " " + allargs).replace(' ', '-')
			print 'Evaluating ' + bench[0] + ": " + allargs + '...'
			command1 = prerun + execfolder + bench[0] + allargs +  " > " + outputfolder + "/" + outputname + " 2>" + outputfolder + "/errors.log"
			print "COMMAND IS %s" % command1
			os.system(command1)
	
	finish = datetime.datetime.now()
	
	return (finish - start).seconds
	

def run(runargs, tdesired_seconds, core_pp_gfs, core_pb_gbs, footprint_megabytes, prerun, execfolder, outputfolder, execprefix):
	print '======================================================================================='
	print '==========================  RUN.PY: BENCHMARKING MPCFKERNELS =========================='
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

	time1 = _runall(10, args, prerun, execfolder, outputfolder, footprint_megabytes)
	print "With one iteration it took " + str(time1) + " seconds."

	if time1 > tdesired_seconds/4.:
		print "Alright we have spent enough time. Benchmarking done"
	else:
		time2 = _runall(100, args, prerun, execfolder, outputfolder, footprint_megabytes)
		print "With two iteration it took " + str(time2) + " seconds."

		#estimate the number of iterations to cover the desired time
		a = max(1e-6, (time2 - time1)/90)
		b = max(0, time1 - 10*a)
		td = tdesired_seconds - time1 - time2 - 30 #safety margin
		niters = int(max(1, min(1e3, math.floor((td - b)/a))))
		print "a=%.2f [s/iteration], b=%.2f [s], td=%.2f [s], niters=%d [iteration]\n" % (a, b, td, niters)
		
		if niters > 100: 
			timefinal = _runall(niters, args, prerun, execfolder, outputfolder, footprint_megabytes)
			print "In total we have spent " + str(time1 + time2 + timefinal) + " seconds."
		else:
			print "Alright we have spent enough time. Benchmarking done"
