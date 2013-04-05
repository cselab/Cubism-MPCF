MYSRCDIR=/tmp/CUBISM-DATA
MYDSTDIR=/tmp/CUBISM-DATA/VPCACHE
QUALITYCHECK=1

function generatechannel()
{
	(( $# == 3 )) || { printf "I was expecting 3 arguments, but i received %b\n" "$#"  ;  exit 2; }
	
	local MYCHANNEL=$1
	local MYMAX=$2
	local MYMIN=$3
	
	for F in ${MYSRCDIR}/*${MYCHANNEL} 
	do
		#ok, lets do it!
		mpirun -n 8 ./vpcachier -simdata $F -vp ${MYDSTDIR}/$(basename $F).vpcache -min $MYMIN -max $MYMAX
		
		echo -n verifying...
		
		local MYTSTART="$(date +%s%N)"
		
		if (( $QUALITYCHECK ))
		then
			./vpcachier -simdata $F -vp ${MYDSTDIR}/$(basename $F).vpcache -min $MYMIN -max $MYMAX  -read
			
			(( $? )) && { printf "Ooops! quality check of the output failed.\n" ; exit 2; } 
		fi
		
		local MYDELTAT_NSEC="$(($(date +%s%N)-MYTSTART))" 
		(( T_CHECK_MSEC = T_CHECK_MSEC + MYDELTAT_NSEC / 1000000 ))
		
		echo quality check passed!
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
generatechannel 6 3.47 4773.44

TOTALTIME="$(($(date +%s)-TSTART))"
T_CHECK=$(( T_CHECK_MSEC / 1000 ))
echo Total time: $TOTALTIME seconds. Time spent in quality checks: $T_CHECK seconds. Ciao!
