#!/bin/bash
if [ -z $MYKERNEL ]; then
	MYKERNEL=all
fi

if [ -z $MYNBLOCKS ]; then
        MYNBLOCKS=1
fi

if [ -z $MYNTIMES ]; then
        MYNTIMES=10
fi

if [ -z $MYPROFILE ]; then
        MYPROFILE=0
fi

   runjob --block $4  \
  -n $1 \
  -p $2 \
  --cwd /gpfs/DDNgpfs1/bekas/BGQ/PETROS/diegor/CUBISM-MPCF/CubismApps/MPCFcore/makefiles \
  --exe /gpfs/DDNgpfs1/bekas/BGQ/PETROS/diegor/CUBISM-MPCF/CubismApps/MPCFcore/makefiles/mpcf-core \
  --envs PAMI_DEVICE=B \
  --envs PAMID_COLLECTIVES=1\
  --envs PAMI_MEMORY_OPTIMIZED=1\
  --envs BG_SHAREDMEMSIZE=128 \
  --envs BG_MEMSIZE=16384 \
  --envs BG_THREADLAYOUT=2 \
  --envs XLSMPOPTS=parthds=$3 \
  --envs OMP_PROC_BIND=TRUE\
  --envs OMP_STACKSIZE=1M \
  --envs OMP_WAIT_POLICY=ACTIVE\
  --args -n $MYNTIMES \
  --args -nblocks $MYNBLOCKS\
  --args -pp 3.2 \
  --args -pb 0.66 \
  --args -profile $MYPROFILE \
  --args -kernel $MYKERNEL \
  --args -accuracy 5e-5

#  --args -kernel Convection_QPX


