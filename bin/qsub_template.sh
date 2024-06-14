#!/bin/bash
#PBS -S /bin/bash
#PBS -N DTT_DOME_BEV_all_eff_cuda
#PBS -o DTT_DOME_BEV_all_eff_cuda.out
#PBS -e DTT_DOME_BEV_all_eff_cuda.err
#PBS -q dtt
#PBS -l select=1:ncpus=1:ngpus=1

module load autoload petsc/3.15_dp_hypre_omp--hpc-sdk--20.11--binary

export KMP_AFFINITY=compact
# export KMP_AFFINITY=verbose,compact
# export KMP_AFFINITY="verbose,granularity=core,compact"
# export KMP_STACKSIZE=64m

export OMP_STACKSIZE=1G
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export OMP_PLACES=cores
# export OMP_PROC_BIND=true

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${NVHPC_HOME}/Linux_x86_64/24.3/cuda/lib64/

${HTHETA}/htheta.exe