import sys
import re
import os

#this is my shit
import run

print 'usage: python benchmark-brutus.py {gcc|icc} <desired-time-in-minutes> <footprint-per-thread-in-megabytes>'

if len(sys.argv) > 4:
	print('Error: passing more arguments than allowed. Aborting now.')
	sys.exit("Aborting because of errors!")
	
print "Args are:"
for x in sys.argv: print x
	
stage = "gcc icc"
if len(sys.argv) > 1:
	if re.search("gcc|icc", sys.argv[1]) == None:
		print('No match for the first argument')
		sys.exit("Aborting because of errors!")
	else:
		stage = sys.argv[1]

tdesired_seconds = 5.*60
if len(sys.argv) >= 3:
	tdesired_seconds = float(sys.argv[2])*60;
	
footprint_megabytes = 128
if len(sys.argv) >= 4:
	footprint_megabytes = int(sys.argv[3]);
	
getnodedata = os.popen('hostname ; date \'+-%y-%m-%d-%H_%M_%S\'')
nodedata = ("result-" + getnodedata.readline() + getnodedata.readline()).replace('\n', '')
common = " source /etc/profile.d/modules.sh; export OMP_NUM_THREADS=48; "

if re.search("gcc", stage) != None:
	run.run("", tdesired_seconds, "17.6", "2.", footprint_megabytes, common + " module load gcc/4.6.1; ", "./brutus-exec-gcc461/", nodedata+"-gcc461", "")
	
if re.search("icc", stage) != None:
	run.run("", tdesired_seconds, "17.6", "2.", footprint_megabytes, common + " module load intel/12; ", "./brutus-exec-icc12/", nodedata+"-icc12", "")
