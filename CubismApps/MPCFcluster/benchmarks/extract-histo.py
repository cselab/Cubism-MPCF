import sys
import operator
import math
import os

print "usage: python dotheshit.py <weak-scaling-folder> <reportfreq>"

folderpath = "."
if len(sys.argv) > 1: folderpath = sys.argv[1].rstrip("/");

print "do the shit on folder %s !" % folderpath

nssteps = 10
if len(sys.argv) > 2: nssteps = int(sys.argv[2]);

tokens = folderpath.split('-');
#nprocs = int(tokens[2])*int(tokens[3])*int(tokens[4])
#print tokens[6].split('n')[0];
#sys.exit(0)
nprocs = int(tokens[6].split('n')[0]) #int(tokens[2])*int(tokens[3])*int(tokens[4])

if nprocs < 0 :
	print "Oops nprocs is %d! Aborting now. " % nprocs
	sys.exit(-1)
else:
	print "Number of processes: %d" % nprocs
	
filenames = [ "hist_STEPID" ]
#, "hist_FLOWSTEP", "hist_UPDATE", "hist_DIFFUSION" , "hist_NSYNCH", "hist_SURFACETENSION" ]

for fname in filenames:
	print "working with %s " % fname
	
	data = {(-1,-1,-1): -1.0}
	counter = 0
	maxval = 0
	minval = 1e8
	myfile = open(folderpath + "/" + fname,"r")
	for line in myfile:
		currentry = float(line.strip("\n"))
		maxval = max(maxval, currentry)
		minval = min(minval, currentry)
		
		ssid = counter % nssteps
		pid = (counter/nssteps) % nprocs
		tid = counter / (nssteps*nprocs)

		#print "line: %s" % currentry
		#print "tid %d pid %d ssid %d" %(tid, pid, ssid)
		
		counter += 1
		
		data[(tid, ssid, pid)] = currentry
		
	myfile.close();
	
	NT = counter / (nssteps*nprocs)
	
	h = (maxval-minval)/nprocs*4
	actualmaxval = maxval
	maxval *= 1.1
	minval = minval - (maxval - actualmaxval)
	
	N =  int(math.ceil(maxval/h))

	for tid in range(0, NT):
		for ssid in range(0, nssteps):
			histo = [0] * N	
			scatterplot = []
			for pid in range(0, nprocs):
				val = data[(tid, ssid, pid)]
				histo[int(math.floor(val/h))] += 1
				scatterplot.append((pid, val))
				
			path2file = folderpath + "/" + fname + "-timestep-" + str(tid*nssteps + ssid) + ".txt"
			myoutput = open(path2file, "w")
			for i in range(0, N): myoutput.writelines("%f %f\n" % ((i+0.5)*h, histo[i]))
			#for entry in scatterplot: myoutput.writelines("%f %f\n" % entry)
			myoutput.close()
			
			gnuplottext = "set terminal png ;\
			set output \""+ path2file.rstrip(".txt") + ".png\" ;\
			set xrange["+ str(minval) +":"+ str(maxval) +"]; \
			set yrange[-0.5:100];\
			plot \""+path2file+"\" using 1:2 w lines ; "
			
			#gnuplottext = "set terminal png ; \
			#set output \"" + path2file.rstrip(".txt") + ".png\" ;\
			#set xrange[-0.5:"+ str(nprocs) +"]; \
			#set yrange["+ str(0) +":"+ str(maxval) +"];\
			#plot \""+path2file+"\" using 1:2 "
			
			gnuplotfile = open("diego.gnuplot", "w");
			gnuplotfile.write(gnuplottext);
			gnuplotfile.close();
			
			os.system("gnuplot diego.gnuplot; ");
	
	os.system("ffmpeg -y -r 30 -b 3000k -i " + folderpath + "/" + fname + "-timestep-%d.png " + 
	folderpath + "-" + str(nprocs) + "-" + fname + ".avi ; ");
	#sys.exit(0)
	
			
print "Histograms done. Ciao!"
