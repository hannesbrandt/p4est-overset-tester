/**
 * \file    wake3d_driver_interface.h
 * \author  akirby
 * \ingroup amr_group
 */

#ifndef WAKE3D_DRIVER_INTERFACE_H
#define WAKE3D_DRIVER_INTERFACE_H

/* header files */
#include "precision_types.h"
#include "amr_var_data.h"
#include "amr_var_defines.h"

#include "amr_mpi.h"
#include "amr_vtk.h"
#include "amr_utilities.h"
#include "amr_external_func.h"
#include "amr_external_func.h"
#include "amr_p4est_utilities.h"
#include "amr_initialize.h"

#include "dg4est_solver.hxx"

#ifdef __cplusplus
extern "C" {
#endif

ctx_t * driver_interface_get_ctx();
p4est_t * driver_interface_get_p4est();
p4est_connectivity_t * driver_interface_get_p4est_conn();

/** Driver interface function:
 *      set the logging threshold value
 *      *Values:
 *      * 0     log everything
 *      * 1     prefix file and line number
 *      * 2     information on the internal state
 *      * 3     information on conditions, decisions
 *      * 4     main information a function is doing
 *      * 5     important consistency/performance information
 *      * 6     few lines for a major API function
 *      * 7     logs a few lines max per program
 *      * 8     logs errors only
 *      * 9     never logs anything
 *
 * @param [in] log_threshold    value to set logging threshold
 */
void driver_interface_set_logging(int log_threshold);

/** Driver interface callback function:
 *      create an MPI communicator using group_comm
 *
 * @param [in] group_comm    Fortran MPI integer to create a C-MPI communicator
 */
void driver_interface_initialize_group_mpi(MPI_Fint *group_comm);

/** Driver interface callback function:
 *      Initialize code and additional driver interface data
 */
void driver_interface_initialize(void);

/** Driver interface callback function:
 *      Deallocate driver interface data
 */
int driver_interface_shutdown(void);

/** Driver interface callback function:
 *      get data from code to wake3d driver
 *
 * @param [out] nnode_out       number of nodes in mesh
 * @param [out] nwbc_out        number of wall bc nodes
 * @param [out] nobc_out        number of outer bc nodes
 * @param [out] ntri_out        number of triangles (always 0)
 * @param [out] nquad_out       number of quads: # 2d quads
 * @param [out] ntet_out        number of tetrahedra (always 0)
 * @param [out] npyramid_out    number of pyramids (always 0)
 * @param [out] nprism_out      number of prisms (always 0)
 * @param [out] nhex_out        number of hexahedron; number of quadrants
 */
void driver_interface_get_data(int *body_tag,int *nnode_out,
                               int *nwbc_out,int *nobc_out,
                               int *ntri,int *nquad,
                               int *ntet_out,int *npyramid_out,
                               int *nprism_out,int *nhex_out);

/** Driver interface callback function:
 *      access to pointers such as solution vector and node locations
 *
 * @param [out] ug_ptr          pointer to solution data
 * @param [out] xgeom_ptr       pointer to node locations 3*nnodes
 * @param [out] iblank_ptr      pointer to node iblank, Tioga uses this and we hold the storage
 * @param [out] iwbcnode_ptr    pointer to array of indices containing walls boundaries
 * @param [out] iobcnode_ptr    pointer to array of indices containing outer boundaries
 * @param [out] ndc4_ptr        pointer to tetrahedra nodes (NULL)
 * @param [out] ndc5_ptr        pointer to pyramid nodes (NULL)
 * @param [out] ndc6_ptr        pointer to prism nodes (NULL)
 * @param [out] ndc8_ptr        pointer to hexahedron nodes 8*num_quadrants
 * @param [out] iblank_cell_ptr pointer to cell iblank array; -1=fringe, 0=hole, 1=field
 * @param [out] corn_ptr        function pointer to count_receptor_nodes routine
 * @param [out] crrn_ptr        function pointer to create_receptor_nodes routine
 * @param [out] dit_ptr         function pointer to donor_inclusion_test routine
 * @param [out] cdf_ptr         function pointer to create_donor_frac routine
 * @param [out] crc_ptr         function pointer to convert_to_receptor_coefficients routine
 */
void driver_interface_set_pointers(Real **ug_ptr,
                                   Real **xgeom_ptr,
                                   int **iblank_ptr,
                                   int **iwbcnode_ptr,
                                   int **iobcnode_ptr,
                                   int **ndcTri_ptr,
                                   int **ndcQuad_ptr,
                                   int **ndcTet_ptr,
                                   int **ndcPyr_ptr,
                                   int **ndcPrism_ptr,
                                   int **ndcHex_ptr,
                                   int **iblank_cell_ptr,
                                   void (**corn_ptr)(int*,int*),
                                   void (**crrn_ptr)(int*,int*,Real*),
                                   void (**dit_ptr)(int*,Real*,int*,Real*),
                                   void (**cdf_ptr)(int*,Real*,int*,int*,Real*,Real*,int*),
                                   void (**crc_ptr)(int*,int*,Real*,int*,int*,Real*));

/** Driver interface callback function:
 *      output solution
 *
 * @param [in] index    used to index the output name
 */
void driver_interface_output_solution(int *index);

/** Driver interface callback function:
 *      point inclusion test
 *
 * @param [in]  npts        number of points to test inclusion
 * @param [in]  pts         points to test inclusion
 * @param [out] cell_id     cell ids of points found
 */
void driver_interface_point_inclusion(int *npts,Real *pts,int *cell_id);

/** Driver interface callback function:
 *      bounding box intersection check
 *
 * @param [in]  nbb         number of bounding boxes
 * @param [in]  bb_xlo      lo coordinates of the bounding boxes
 * @param [in]  bb_xhi      hi coordinates of the bounding boxes
 * @param [out] bb_flag     flag indicator if intersecting bounding box
 */
void driver_interface_bounding_box_intersection(int *nbb,
                                                Real *bb_xlo,
                                                Real *bb_xhi,
                                                int *bb_flag);

/** Driver interface callback function:
 *      basis and its derivative evaluation
 *
 * @param [in]  pdegree     basis degree
 * @param [in]  npts        number of points to evaluate basis and derivative
 * @param [in]  xi_pts      points to evaluate basis and derivative
 * @param [out] phi         basis evaluated at the points
 * @param [out] dphid       derivative of basis evaluated at the points
 */
void driver_interface_basis_and_derivative(int *pdegree,int *npts,Real *xi_pts,
                                           Real *phi,Real *dphi);

/** Driver interface callback function:
 *      print stuff
 *
 * @param [in] value    integer to be displayed to terminal
 */
void driver_interface_print_stuff(int *value);

/** Driver interface callback function:
 *      callback function for Tioga that counts receptor nodes
 *
 * @param [in]  cell_index          cell index local to processor
 * @param [out] num_receptor_nodes  number of receptor nodes for cell with index \a cell_index
 */
void wake3d_count_receptor_nodes(int *cell_index,int *num_receptor_nodes);

/** Driver interface callback function:
 *      callback function for Tioga that creates receptor nodes
 *
 * @param [in]  cell_index          cell index local to processor
 * @param [in]  num_receptor_nodes  number of receptor nodes for cell with index \a cell_index can be used to allocate
 * @param [out] receptor_nodes      location of receptor nodes 3*num_receptor_nodes
 */
void wake3d_create_receptor_nodes(int *cell_index,int *num_receptor_nodes,Real *receptor_nodes);

/** Driver interface callback function:
 *      callback function for Tioga that tests high order inclusion of a point with a cell
 *
 * @param [in]  cell_index  index of the cell testing for inclusion
 * @param [in]  pt_xyz      physical location of point used for inclusion (length 3)
 * @param [out] pass_flag   set to 1 if \a xyz is inside cell and 0 if it is not
 * @param [out] rst_out     if \a pass_flag is 1 then put the natural coordinates here (-1<r,s,t<1, length 3)
 */
void wake3d_donor_inclusion_test(int *cell_index,Real *pt_xyz,int *pass_flag,Real *rst_out);

/** Driver interface callback function:
 *      callback function for Tioga that creates interpolation weights
 *
 * @param [in]  cell_index  index of cell creating interpolation weights for
 * @param [in]  xyz         physical location (not used anymore since \a rst_in is always passed in)
 * @param [out] nfrac       number of modes in the basis used for interpolation
 * @param [out] index       index to access solution array index is base 1
 * @param [out] frac        basis evaluated at point \a rst_in length \a nfrac
 * @param [in]  rst_in      natural coordinates input to evaluate the basis at  (-1<r,s,t<1, length 3)
 * @param [in]  dimfrac     large integer (set inside Tioga to be 11^3=1331 for a p=10 hex) used for allocation frac if necessary
 */
void wake3d_create_donor_frac(int *cell_index,Real *xyz,int *nfrac,int *index,
                              Real *frac,Real *rst_in,int *dimfrac);

/** Driver interface callback function:
 *      callback function for Tioga that converts nodal points into solution coefficients
 *      * nodal DG code this is simply a copy
 *
 * @param [in]  cell_index  index of cell
 * @param [in]  npts        number of points to convert
 * @param [in]  f           value of nodal function to be converted
 * @param [out] tm          number of nodes/modes being passed out
 * @param [out] index_out   index to modify solution array base 1
 * @param [out] a_out       solution coefficients
 */
void wake3d_convert_to_receptor_coefficients(int *cell_index,int *npts,Real *f,
                                             int *tm,int *index_out,Real *a_out);

#ifdef __cplusplus
}
#endif
#endif /* WAKE3D_DRIVER_INTERFACE_H */
