#MYSRCDIR=/tmp/CUBISM-DATA
#MYDSTDIR=/tmp/CUBISM-DATA/HDFCACHE
MYSRCDIR=/home/chatzidp/LIQ_BGQ/10042013A/CUBISM-MPCF/CubismApps/tools/vp2hdf/tmpin
MYDSTDIR=/home/chatzidp/LIQ_BGQ/10042013A/CUBISM-MPCF/CubismApps/tools/vp2hdf/tmp
QUALITYCHECK=0

#get the number of threads
{
	printf "inquiring number of cores..."
	
	if [ $(uname) == "Darwin" ]
	then
		THREADS=$(sysctl hw.ncpu | cut -d" " -f2) 
	else
		{ THREADS=$(egrep -c "core id" /proc/cpuinfo 2> /dev/null) ; }  || \
		{ printf "not able to recover number of cores\n" ; THREADS=1; }
	fi
	
	printf "setting it to %d.\n" $THREADS
}

function generatechannel()
{
	(( $# == 3 )) || { printf "I was expecting 3 arguments, but i received %b\n" "$#"  ;  exit 2; }
	
	local MYCHANNEL=$1
	local MYMAX=$2
	local MYMIN=$3
	
	for F in ${MYSRCDIR}/*channel${MYCHANNEL} 
	do
		#ok, lets do it!
		mpirun -mca btl sm,self -n $THREADS ./vp2hdf -simdata $F -h5file ${MYDSTDIR}/$(basename $F)
		#echo $F
	done
}

T_CHECK_MSEC=0
TSTART="$(date +%s)"

generatechannel 0 0.6 1e3

#it does not make sense to process velocity, so i skip this
#generatechannel 1 -0.1 0.1
#generatechannel 2 -0.1 0.1
#generatechannel 3 -0.1 0.1

generatechannel 4 0.02 100.0
generatechannel 5 0.17 2.5; 
#generatechannel 6 3.47 4773.44

TOTALTIME="$(($(date +%s)-TSTART))"
T_CHECK=0
echo Total time: $TOTALTIME seconds. Time spent in quality checks: $T_CHECK seconds. Ciao!
