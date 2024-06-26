#!/bin/bash

#SBATCH --nodes=1                    # nodes
#SBATCH --ntasks-per-node=1          # tasks per node
#SBATCH --cpus-per-task=1            # cores per task
#SBATCH --gres=gpu:1                 # GPUs per node
#SBATCH --mem=494000                 # mem per node (MB)
#SBATCH --time=00:20:00              # time limit (d-hh:mm:ss)
#SBATCH --account=tra24_openhack     # account ($ saldo -b)
#SBATCH --partition=boost_usr_prod   # partition name
#SBATCH --reservation=s_tra_hack
#SBATCH --qos=normal                 # quality of service

module load nvhpc/24.3

export KMP_AFFINITY=compact
# export KMP_AFFINITY=verbose,compact
# export KMP_AFFINITY="verbose,granularity=core,compact"
# export KMP_STACKSIZE=64m

export OMP_STACKSIZE=1G
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export OMP_PLACES=cores
# export OMP_PROC_BIND=true

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${NVHPC_HOME}/Linux_x86_64/24.3/cuda/lib64/

#srun ./htheta.exe
#srun nvprof -o prof_out ./htheta.exe
#ncu --set=full -o myreport --import-source=yes --launch-skip=2*5 --launch-count=5 ./htheta.exe
nsys profile -t cuda,nvtx,openacc ./htheta.exe