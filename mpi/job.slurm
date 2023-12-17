#!/bin/bash
#SBATCH --job-name=f90-mpi       # create a short name for your job
#SBATCH --nodes=2                # node count
#SBATCH --ntasks=32              # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1G         # memory per cpu-core (4G is default)
#SBATCH --time=00:00:50          # total run time limit (HH:MM:SS)


module purge
module load intel-mpi/gcc/2021.3.1
srun ./block_under_tension_mpi
