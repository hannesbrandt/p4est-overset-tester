include_directories(${p4est_dir}/include)
include_directories(${p4est_home}/example/multi)

find_library(p4est_lb NAMES p4est PATHS ${p4est_dir}/lib)
find_library(sc_lb    NAMES sc    PATHS ${p4est_dir}/lib)

file(GLOB PDRIVER_H  *.h *.hxx)

set(SRC
    input.cxx
    flow_solver.cxx
    p4est_overset_interface.cxx
    driver.cxx
)

add_library(pdriver SHARED ${SRC})
target_link_libraries(pdriver ${sc_lb} ${p4est_lb} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

add_executable(pdriver.mpi main.cxx)
target_link_libraries(pdriver.mpi pdriver)

install(TARGETS pdriver pdriver.mpi
        RUNTIME DESTINATION bin/
        LIBRARY DESTINATION lib/
        ARCHIVE DESTINATION lib/
)

install(FILES ${PDRIVER_H} DESTINATION include)
