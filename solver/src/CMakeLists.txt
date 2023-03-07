include_directories(${p4est_dir}/include)

find_library(p4est_lb NAMES p4est PATHS ${p4est_dir}/lib)
find_library(sc_lb    NAMES sc    PATHS ${p4est_dir}/lib)

file(GLOB AMR_H  *.h *.hxx)

set(AMR_SRC
    amr_var_data.cxx
    amr_utilities.cxx
    amr_p4est_utilities.cxx
    amr_tag.cxx
    amr_mpi.cxx
    amr_vtk.cxx
    amr_postprocess.cxx
    amr_external_func.cxx
    amr_initialize.cxx
    amr_transfer_data.cxx
    amr_regrid.cxx
    amr_regrid_solution.cxx
    amr_gridfile.cxx
    amr_containers.cxx
    amr_geometry.cxx
    wake3d_driver_interface.cxx
)

add_library(ovtest2d SHARED ${AMR_SRC} ${AMR_H} ${GIT_FILE})
add_library(ovtest3d SHARED ${AMR_SRC} ${AMR_H} ${GIT_FILE})
target_compile_definitions(ovtest2d PUBLIC -D_2D_)

target_link_libraries(ovtest2d ${sc_lb} ${p4est_lb} m)
target_link_libraries(ovtest3d ${sc_lb} ${p4est_lb} m)

install(TARGETS ovtest2d ovtest3d
        RUNTIME DESTINATION bin/
        LIBRARY DESTINATION lib/
        ARCHIVE DESTINATION lib/
)

install(FILES ${AMR_H} DESTINATION include)
