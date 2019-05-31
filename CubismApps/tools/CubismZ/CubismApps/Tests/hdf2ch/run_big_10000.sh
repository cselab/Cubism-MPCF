#!/bin/bash -l
#

#SBATCH --job-name="big64"
#SBATCH --output=cubismz-%j.txt
#SBATCH --error=cubismz-%j.txt
#SBATCH --time=12:00:00
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=gpu
#SBATCH --account=s659

#!/bin/bash
set -x #echo on

source setup_daint.sh

##################
#h5file=../../../../fabdata/big/data_010000-p-big.h5
#./runone_big.sh $h5file

##################
h5file=../../../../fabdata/big/data_010000-rho-big.h5
./runone_big.sh $h5file

##################
h5file=../../../../fabdata/big/data_010000-E-big.h5
./runone_big.sh $h5file

##################
h5file=../../../../fabdata/big/data_010000-a2-big.h5
./runone_big.sh $h5file

##################
#h5file=../../../../fabdata/data_010000-divU-big.h5
#./runone_big.sh $h5file

##################
#h5file=../../../../fabdata/data_010000-Omegax-big.h5
#./runone_big.sh $h5file

##################
#h5file=../../../../fabdata/data_010000-Ux-big.h5
#./runone_big.sh $h5file

