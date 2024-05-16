#!/bin/bash

#SBATCH --nodes=1                    # nodes
#SBATCH --ntasks-per-node=1          # tasks per node
#SBATCH --cpus-per-task=1            # cores per task
#SBATCH --gres=gpu:1                 # GPUs per node
#SBATCH --mem=494000                 # mem per node (MB)
#SBATCH --time=24:00:00              # time limit (d-hh:mm:ss)
#SBATCH --account=tra24_openhack_0   # account ($ saldo -b)
#SBATCH --partition=gll_all_serial   # partition name
#SBATCH --qos=normal                 # quality of service

module purge
module load autoload nvhpc/22.3

export KMP_AFFINITY=compact
# export KMP_AFFINITY=verbose,compact
# export KMP_AFFINITY="verbose,granularity=core,compact"
# export KMP_STACKSIZE=64m

export OMP_STACKSIZE=1G
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# export OMP_PLACES=cores
# export OMP_PROC_BIND=true

srun ./htheta.exe