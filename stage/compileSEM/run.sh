#!/bin/sh
#SBATCH -J TradWM
#SBATCH -p skylake
#SBATCH -n 128
#SBATCH -N 4
#SBATCH -A b401
#SBATCH -t 05:00:00
#SBATCH -o ./%N.%x.out
#SBATCH -e ./%N.%x.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=de-abreu-correa@lma.cnrs-mrs.fr
#SBATCH --mem=160G

# chargement des modules
module purge
module load userspace/all
module load cmake/3.20.2
module load intel-compiler/64/2018.0.128
module load intel-mpi/64/2018.0.128
module load intel-mkl/64/2018.0.128
module load hdf5/icc18/impi/1.10.1

### export I_MPI_LINK=opt_mt
export CC=icc CXX=icpc FC=ifort
export I_MPI_CC=icc I_MPI_CXX=icpc I_MPI_FC=ifort I_MPI_F90=ifort
export EXTRAMPI=/trinity/shared/apps/tr17.10/x86_64/intel-2018.0.128/compilers_and_libraries_2018.0.128/linux/mpi/intel64
export HDF5_PATH=/trinity/shared/apps/tr17.10/x86_64/hdf5-icc18-impi-1.10.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/trinity/shared/apps/tr17.10/x86_64/hdf5-icc18-impi-1.10.1/lib
export LIBRARY_PATH=$LIBRARY_PATH:/trinity/shared/apps/tr17.10/x86_64/hdf5-icc18-impi-1.10.1/lib
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so


# echo of commands
set -x

# To compute in the submission directory
cd ${SLURM_SUBMIT_DIR}

# execution with 'ntasks' MPI processes
srun -n $SLURM_NTASKS /home/${SLURM_JOB_USER}/SEM/buildSEMRea/SEM3D/sem3d.exe > output.log
