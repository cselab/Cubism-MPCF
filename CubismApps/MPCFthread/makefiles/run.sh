#!/bin/bash
   runjob --block $4  \
  -n $1 \
  -p $2 \
  --cwd /gpfs/DDNgpfs1/bekas/BGQ/PETROS/CUBISM-MPCF/CubismApps/MPCFthread/makefiles \
  --exe /gpfs/DDNgpfs1/bekas/BGQ/PETROS/CUBISM-MPCF/CubismApps/MPCFthread/makefiles/weno \
  --envs PAMI_DEVICE=B \
  --envs PAMID_COLLECTIVES=1\
  --envs PAMI_MEMORY_OPTIMIZED=1\
  --envs BG_SHAREDMEMSIZE=128 \
  --envs BG_MEMSIZE=16384 \
  --envs BG_THREADLAYOUT=2 \
  --envs XLSMPOPTS=parthds=$3 \
  --envs OMP_SCHEDULE=$5\
  --envs OMP_WAIT_POLICY=ACTIVE\
  --envs L1P_POLICY=std\
  --envs USE_FAST_WAKEUP=TRUE
