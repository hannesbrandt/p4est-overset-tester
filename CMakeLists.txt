CMAKE_MINIMUM_REQUIRED(VERSION 2.17)
PROJECT(pdriver)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CC_COMPILER mpicc)

# set C++17 standard for OCCA compatibility
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# --------------------------------------------------------------------------- #
# force GNU linker to inlude all shared libraries
if (NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed")
else ()
  set(CMAKE_MACOSX_RPATH ON)
endif ()

message(${CMAKE_SOURCE_DIR})
set(p4est_home ${CMAKE_SOURCE_DIR}/3PL/p4est CACHE PATH "p4est home directory location")
set(p4est_dir  ${CMAKE_SOURCE_DIR}/3PL/p4est/local CACHE PATH "p4est install directory location")
# ======================================= #
# Library header files to be included
include_directories(${p4est_dir}/include
                    ${p4est_home}/example/multi)

# Libraries to be linked
find_library(p4est_lb     NAMES p4est           PATHS ${p4est_dir}/lib)
find_library(sc_lb        NAMES sc              PATHS ${p4est_dir}/lib)
# ======================================= #

# ======================================= #
# Source directories for this project
set(src_path src)
add_subdirectory(${src_path})
# ======================================= #

