#!/bin/bash

set -x #echo on
ulimit -s 1600000000
export OMP_NUM_THREADS=1
rm tmp00000.StreamerGridPointIterative.channel0;

wt=3
#thrval=0.0001
thrval=$1
inputfile=../../../../fabdata/data_010000-p.h5
nb=16

#inputfile=../../../../fabdata/data_010000-p-med.h5
#nb=32

export CUBISMZ_NOIO=1

export OMP_PROC_BIND=TRUE
export OMP_NUM_THREADS=1;  srun --ntasks=1 -c 1   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt
export OMP_NUM_THREADS=2;  srun --ntasks=1 -c 2   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt
export OMP_NUM_THREADS=4;  srun --ntasks=1 -c 4   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt
export OMP_NUM_THREADS=6;  srun --ntasks=1 -c 6   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt
export OMP_NUM_THREADS=8;  srun --ntasks=1 -c 8   ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt
export OMP_NUM_THREADS=10; srun --ntasks=1 -c 10  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt
export OMP_NUM_THREADS=12; srun --ntasks=1 -c 12  ./hdf2ch -xpesize 1 -bpdx $nb -bpdy $nb -bpdz $nb -sim io -simdata $inputfile -outdata tmp -threshold $thrval -wtype_write $wt
