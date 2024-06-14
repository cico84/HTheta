# Load required modules for compilation
module purge
module load autoload petsc/3.15_dp_hypre_omp--hpc-sdk--20.11--binary

# Compile code, if necessary
make clean -f "${HTHETA}/src/makefile_openacc_eni.txt"
make HTHETA -f "${HTHETA}/src/makefile_openacc_eni.txt"

# Launch sbatch executable (a file named "sbatch_sim_name.sh" must be available in the bin folder)
qsub "qsub.sh"