#!/bin/bash

# ========== #
# make p4est #
# ========== #
#cd 3PL/p4est
#git checkout feature-multi
#git submodule init && git submodule update
#./bootstrap
#./configure --enable-shared \
#            --enable-mpi    \
#            --without-blas  \
#            | tee config.out
#make; make install
#cd ../..

# ============ #
# make pdriver #
# ============ #
MY_PATH="$( cd "$( dirname "$0" )" && pwd )"
P4EST_DIRECTORY=${MY_PATH}/3PL/p4est

mkdir build
cd build
cmake -D p4est_dir=${P4EST_DIRECTORY}/local \
      -D p4est_home=${P4EST_DIRECTORY} \
	 -G "Unix Makefiles" ..
make
cp ../example/input.pdriver bin/.
mkdir bin/lib
cd ..
ln -s build/bin
