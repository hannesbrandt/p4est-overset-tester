#!/bin/bash

#================#
# config options #
#================#
CLEAN_DIST=0
CONF_P4EST=0
BUILD_P4EST=0
CONF_SOLVER=0
BUILD_SOLVER=0
DEBUG_ON=1
BUILD_TYPE="Debug"

#===========#
# compilers #
#===========#
CC=mpicc
CXX=mpicxx

#CC=mpiicc
#CXX=mpiicpc

# compiler paths #
CC_PATH="`which $CC`"
CXX_PATH="`which $CXX`"
LD_PATH="`which ld`"

help() {
    echo " +--------- p4est-overset-tester ---------+"
    echo "  Configure + Build: ./makescript.sh -go   "
    echo "  ======================================== "
    echo "    Configure p4est: ./makescript.sh -cp   "
    echo "        Build p4est: ./makescript.sh -bp   "
    echo "   Configure solver: ./makescript.sh -cs   "
    echo "       Build solver: ./makescript.sh -bs   "
    echo " +----------------------------------------+"
}


if [[ $# -lt 1 ]]; then
  help
  exit 1
fi

for var in "$@"
do
  if [ "$var" == "--help" -o "$var" == "-help" -o "$var" == "-h" ]; then
    help
    exit 0

  elif [ "$var" == "-distclean" -o "$var" == "-dc" ]; then
    echo "Found known argument: $var"
    echo "Cleaning the distribution..."
    CLEAN_DIST=1

  elif [ "$var" == "-go" ]; then
    echo -e "Found known argument: ${gC}$var${eC}"
    CONF_P4EST=1
    CONF_SOLVER=1
    BUILD_P4EST=1
    BUILD_SOLVER=1

  elif [ "$var" == "-cp" ]; then
    echo "Found known argument: $var"
    echo "Configuring p4est..."
    CONF_P4EST=1

  elif [ "$var" == "-bp" ]; then
    echo "Found known argument: $var"
    echo "Building p4est..."
    BUILD_P4EST=1

  elif [ "$var" == "-cs" ]; then
    echo "Found known argument: $var"
    echo "Configuring solver..."
    CONF_SOLVER=1

  elif [ "$var" == "-bs" ]; then
    echo -e "Found known argument: $var"
    BUILD_SOLVER=1

  else
    echo "Unknown option: $var"
    echo "See available options: ./makescript.sh -help"
    exit
  fi
done

# ========================= #
# display command line args
# ========================= #
echo "$0 $@"
cmd_args="${@:1}"

# =================================================================== #
if [ $CLEAN_DIST == 1 ]; then
  rm -rf build
  exit 0
fi
# =================================================================== #

# ========== #
# make p4est #
# ========== #

cd 3PL/p4est
if [ $CONF_P4EST == 1 ]; then
    git submodule init
    git submodule update
    git checkout feature-multi
    git submodule init && git submodule update
    ./bootstrap
fi

if [ $BUILD_P4EST == 1 ]; then
    if [ $DEBUG_ON == 1 ]; then
        ./configure --enable-debug      \
                    --enable-mpi        \
                    --enable-shared     \
                    --disable-static    \
                    "CC=$CC" "CXX=$CXX" \
                    "CFLAGS=-Wall -O0 -g -Wextra -Wno-unused-parameter" \
                    | tee config.out
    else 
        ./configure --enable-shared  \
                    --disable-static \
                    --enable-mpi     \
                    --without-blas   \
	                "CC=$CC" "CXX=$CXX" \
                    | tee config.out
    fi
    make;
    make install
fi
cd ../..

# ============ #
# make pdriver #
# ============ #
MY_PATH="$( cd "$( dirname "$0" )" && pwd )"
P4EST_DIRECTORY=${MY_PATH}/3PL/p4est
P4EST_DIR_INSTALL=${P4EST_DIRECTORY}/local

if [ $CONF_SOLVER == 1 ]; then

    if [ ! -d "${P4EST_DIR_INSTALL}" ]; then
        echo "Error:"
        echo "${P4EST_DIR_INSTALL} does not exist."
        echo " Please build p4est first: ./makescript.sh -bp"
        exit
    fi

    mkdir -p build
    cd build
    cmake -D CMAKE_C_COMPILER=${CC_PATH}        \
          -D CMAKE_CXX_COMPILER=${CXX_PATH}     \
          -D CMAKE_LINKER=${LD_PATH}            \
          -D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
          -D p4est_dir=${P4EST_DIRECTORY}/local \
          -D p4est_home=${P4EST_DIRECTORY}      \
          -G "Unix Makefiles" .. | tee cmake_config.out
    cd ..
fi

if [ $BUILD_SOLVER == 1 ]; then

    if [ ! -d "${P4EST_DIR_INSTALL}" ]; then
        echo "Error:"
        echo "${P4EST_DIR_INSTALL} does not exist."
        echo " Please build p4est first: ./makescript.sh -p4est"
        exit
    fi

    cd build
    make
    cp -r ../example/* bin/.
    cd ..
    ln -s -f build/bin
fi
