/**
 * File:   p4est_overset_interface.cxx
 * Author: akirby
 *
 * Created on February 25, 2023, 1:00 PM
 */

/* header files */
#include "p4est_overset_interface.h"

/** Execute multi-mesh overset algorithm.
 * \param [in] glocomm          Global communicator over all meshes.
 * \param [in] leadercomm       If global rank is first for a mesh, a
 *                              communicator over all first mesh ranks.
 * \param [in] rolecomm         Separate communicator over each mesh.
 * \param [in] mesh_index       Index of mesh: 0 for background mesh,
 *                              starting from 1 for overset meshes.
 * \param [in] num_meshes       Number of meshes including background.
 * \param [in] mesh_offsets     Array of ascending global ranks,
 *                              one for the first of each mesh, and then
 *                              one more for the end (exclusive of the last).
 * \param [in] bgp4est          For \a myrole zero, the background forest.
 *                              NULL otherwise.
 * \param [in] qpoints          Query points: 4-double-tuples (x, y, z, v).
 *                              \b a is -1 for the overset boundary, -2 for
 *                              the wall boundary, and a non-negative
 *                              representative volume of the point otherwise.
 *                              The points are provided process-local over
 *                              all overset meshes present on this process.
 *                              NULL for \a myrole zero.
 * \param [in] callback         Placeholder for callback to be added.
 * \param [in] user             Context data passed to callback.
 */
void p4est_overset_interface::p4est_overset_init(MPI_Comm glocomm,
                                                 MPI_Comm leadercomm,
                                                 MPI_Comm rolecomm,
                                                 int mesh_index,
                                                 int num_meshes,
                                                 const int *mesh_offsets,
                                                 p4est_t *bgp4est,
                                                 sc_array_t *qpoints){
  p4est_multi_overset(glocomm,leadercomm,rolecomm,
                      mesh_index,num_meshes,mesh_offsets,
                      bgp4est,qpoints,
                      NULL,NULL);
}