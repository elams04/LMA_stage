#!/bin/bash
module purge
module load userspace/all
module load cmake/3.11.2
module load blas/gcc72/3.7.1
module load hdf5/gcc72/openmpi/1.10.1
module load fftw3/gcc72/openmpi/3.3.6-pl2

make -f makefileRFamu clean
make -f makefileRFamu all
