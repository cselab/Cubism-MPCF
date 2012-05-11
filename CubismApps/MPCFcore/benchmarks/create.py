import sys;
import os;
import re;
import datetime;

def mymerge(extra, n):
	return "".join([extra[s] for s in range(0,n+1)])

def create_executables(commands_precompilation, compileopt, prefix_exe, extra, finalfolder):
	print '======================================================================================='
	print '===========  CREATE.PY: CREATING EXECUTABLES FOR BENCHMARKING MPCFKERNELS ============='
	print '======================================================================================='
	start = datetime.datetime.now()

	#create finalfolder if necessary
	if not os.path.exists(finalfolder):
		os.makedirs(finalfolder)	 
		
	if compileopt != "": print "Compilation additional flags are %s" % compileopt

	useavx=""
	if re.search("avx=1", compileopt) != None:
		useavx="-avx"
			
	#generate executables
	blocksizes = [16, 24, 32];#[16, 24, 28, 32, 36, 40]

	build_dir = "../makefiles"
	cwd = os.getcwd()
	
	try:
		os.chdir(build_dir);
		for ie in range(0, len(extra)): 
			for ibs in range(0, len(blocksizes)):
				command1 = commands_precompilation + " make clean ; " + "make " + compileopt + " bs=" + str(blocksizes[ibs]) + mymerge(extra, ie);
				print 'COMMAND1 is %s' % (command1)
				os.system(command1)
				newname = prefix_exe + "bench-" + "b" + str(blocksizes[ibs]) + "-e" + str(ie) + useavx
				
				command2 = "mv mpcfkernels " + cwd + "/" + finalfolder + "/" + newname
				print 'COMMAND2 is %s' % (command2)
				os.system(command2)
	finally:
		os.chdir(cwd)
	
	finish = datetime.datetime.now()	
	return (finish - start).seconds
	