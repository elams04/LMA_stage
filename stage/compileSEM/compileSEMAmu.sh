#!/bin/bash
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

## MPI_C_COMPILER_INCLUDE_DIRS:STRING=/trinity/shared/apps/tr17.10/x86_64/intel-2018.0.128/compilers_and_libraries_2018.0.128/linux/mpi/intel64/include/gfortran/5.1.0;/trinity/shared/apps/tr17.10/x86_64/intel-2018.0.128/compilers_and_libraries_2018.0.128/linux/mpi/intel64/include

make clean

cmake -DCMAKE_BUILD_TYPE=Debug -DDEBUG=ON -DCMAKE_Fortran_FLAGS_DEBUG="-g -rdynamic" -DCMAKE_C_FLAGS_DEBUG="-g -rdynamic" -DCMAKE_CXX_FLAGS_DEBUG="-g -rdynamic" -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_C_FLAGS="-std=gnu99" -DCMAKE_Fortran_FLAGS="-O0" -DCMAKE_EXE_LINKER_FLAGS=-lhdf5 -DOPT_CPML=0 -DOPT_VEC=1 -DCMAKE_VERBOSE_MAKEFILE=ON -DMPI_CXX_ADDITIONAL_INCLUDE_DIRS=$EXTRAMPI/include -DMPI_C_ADDITIONAL_INCLUDE_DIRS=$EXTRAMPI/include -DHDF5_hdf5_fortran_LIBRARY_RELEASE=$HDF5_PATH/lib/libhdf5_fortran.a -DHDF5_hdf5_LIBRARY_RELEASE=$HDF5_PATH/lib/libhdf5.a -DMPI_Fortran_ADDITIONAL_INCLUDE_DIRS=$EXTRAMPI/include -DMPI_C_COMPILER_INCLUDE_DIRS=$EXTRAMPI/include -DMPI_CXX_COMPILER_INCLUDE_DIRS=$EXTRAMPI/include ../

make -j 1
