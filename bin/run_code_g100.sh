# Load required modules for compilation
module purge
module load nvhpc/24.3

# Compile code, if necessary
make clean -f "${HTHETA}/src/makefile.txt"
make HTHETA -f "${HTHETA}/src/makefile.txt"

# Launch sbatch executable (a file named "sbatch_sim_name.sh" must be available in the bin folder)
sbatch "sbatch.sh"