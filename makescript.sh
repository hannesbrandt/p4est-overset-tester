#!/bin/bash

# compilers
CXX=mpicxx
CC=mpicc

# path to compilers
CC_PATH="`which $CC`"
CXX_PATH="`which $CXX`"

# ========== #
# make p4est #
# ========== #
cd 3PL/p4est
git checkout feature-multi
git submodule init && git submodule update
./bootstrap
./configure --enable-shared  \
	    --disable-static \
            --enable-mpi     \
            --without-blas   \
	    "CC=$CC" "CXX=$CXX" \
            | tee config.out
make; make install
cd ../..

# ============ #
# make pdriver #
# ============ #
MY_PATH="$( cd "$( dirname "$0" )" && pwd )"
P4EST_DIRECTORY=${MY_PATH}/3PL/p4est

mkdir -p build
cd build
cmake -D CMAKE_C_COMPILER=${CC_PATH}        \
      -D CMAKE_CXX_COMPILER=${CXX_PATH}     \
      -D p4est_dir=${P4EST_DIRECTORY}/local \
      -D p4est_home=${P4EST_DIRECTORY} \
      -G "Unix Makefiles" ..

make
cp -r ../example/* bin/.
cd ..
ln -s -f build/bin
