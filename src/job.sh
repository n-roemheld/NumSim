#!/bin/bash
#
#SBATCH --job-name=submission
#SBATCH --output=result.txt
#
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=4
#SBATCH --time=10:00

module use /usr/local.nfs/sgs/modulefiles
module load vtk/8.2
module load gcc/8.2
module load cmake/cmake-3.12.3

srun -n 2 ../build/src/numsim_parallel lid_driven_cavity.txt

