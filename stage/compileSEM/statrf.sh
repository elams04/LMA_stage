#!/bin/sh
#SBATCH -J RFMat
#SBATCH -p skylake
#SBATCH -n 1
#SBATCH -A b401
#SBATCH -t 0:30:00
#SBATCH -o ./%N.%x.out
#SBATCH -e ./%N.%x.err
# #SBATCH --mail-type=BEGIN,END
# #SBATCH --mail-user=XXX@XXX

# chargement des modules
module purge
module load userspace/all
module load cmake/3.11.2
module load blas/gcc72/3.7.1
module load hdf5/gcc72/openmpi/1.10.1
module load fftw3/gcc72/openmpi/3.3.6-pl2

# echo of commands
set -x

# To compute in the submission directory
cd ${SLURM_SUBMIT_DIR}

# execution with 'ntasks' MPI processes
srun -n $SLURM_NTASKS /home/${SLURM_JOB_USER}/SEM/buildRF/statistics.exe > outputmatStat.log
