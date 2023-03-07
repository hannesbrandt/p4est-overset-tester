set(p4est_dir    "" CACHE PATH "p4est install directory location")

# ============================================= #
# Solver Mesh Ghost Connectivity Type:
#  CONN_TYPE:
#    1-Connect faces only
#    2-Connect faces and edges only
#    3-Connect faces, edges, and corners [full]
# ============================================= #
set(CONN_TYPE 1)
if(CONN_TYPE)
  add_definitions(-DCONN_TYPE=${CONN_TYPE})
endif(CONN_TYPE)

# ======================================= #
# Set number of fields per solution point #
set(NFIELDS2D 1 CACHE STRING "2D Simulation number of fields")
set(NFIELDS3D 1 CACHE STRING "3D Simulation number of fields")

if(NFIELDS2D)
  add_definitions(-DNFIELDS2D=${NFIELDS2D})
endif(NFIELDS2D)

if(NFIELDS3D)
  add_definitions(-DNFIELDS3D=${NFIELDS3D})
endif(NFIELDS3D)
# ======================================= #

# ======================================= #
# Library header files to be included
include_directories(${p4est_dir}/include
                    ./include)

# Libraries to be linked
find_library(p4est_lb     NAMES p4est           PATHS ${p4est_dir}/lib)
find_library(sc_lb        NAMES sc              PATHS ${p4est_dir}/lib)
# ======================================= #

# ======================================= #
# Source directories for this project
set(src_path src)
add_subdirectory(${src_path})
# ======================================= #
