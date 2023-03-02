/**
 * \file   amr_gridfile.h
 * \author mugolott
 */

#ifndef AMR_GRIDFILE_H
#define AMR_GRIDFILE_H

/* header files */
#include "amr_var_defines.h"
#include "precision_types.h"
#include "dg4est_solver.hxx"
#include "amr_utilities.h"
#include "amr_containers.h"
#include "amr_geometry.h"
#include "amr_external_func.h"

#ifndef DOXYGEN_IGNORE
#  include <math.h>
#  include <stdio.h>
#  include <stdlib.h>
#  include <stdint.h>
#  include <string.h>
#  include <unistd.h>
#  include <sys/types.h>
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <err.h>
#endif

#ifndef P4_TO_P8
#  include <p4est_extended.h>
#else
#  include <p8est_extended.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Read line
 *
 * @param [in]  stream    File stream
 */
char *get_line(FILE *stream,int ln);

/** Check if periodic inputs are provided in input file
 */
char check_periodic_inputs(ctx_t *ctx);

/** Grid file reader
 *
 * @param [in]  filename    Name of the grid file
 */
p4est_connectivity_t *gridfile_reader(ctx_t *ctx,const char *filename);

/** Grid file reader: reads all the quadrant's coordinates
 *
 * @param [in]  which_tree  quadrant ID
 * @param [out] xyz         quadrant coordinates
 * @return q-degree         quadrant geometry degree
 */
int gridfile_read_highorder(p4est_topidx_t which_tree,Real *xyz);

/** Read the GMSH gridfile data
 *
 * @param [in]  stream          file stream of grid file
 * @param [in]  ctx             simulation context data
 * @param [out] periodic_v2v    periodic vertex-to-vertex map
 * @return connectivity         p4est_connectivity data structure
 */
p4est_connectivity_t* gridfile_reader_msh_stream(FILE *stream,ctx_t *ctx,
                                                 p4est_topidx_t **periodic_v2v);

/** Allocate gridfile data to bookkeep line number for each element
 *
 * @param [in] gfile           grid file data
 * @param [in] num_trees       number of trees
 */
void gridfile_allocate_data(gridfile_t *gfile,p4est_topidx_t num_trees);

/** Deallocate gridfile data
 *
 * @param [in]     ctx             Context data
 */
void gridfile_deallocate_data(ctx_t *ctx);

/** Broadcasts gridfile data to all mpi ranks
 *
 * @param [in] ctx              context data
 * @param [in] num_trees        number of trees in mesh
 */
void gridfile_broadcast_data(ctx_t *ctx,p4est_topidx_t num_trees);

/** Constructs periodic boundaries for unstructured meshes
 *
 * @param [in]  ctx     context data
 * @param [in]  conn    p4est connectivity
 * @param [out] pv2v    periodic vertex-to-vertex map
 */
void gridfile_periodic_boundary_construction(ctx_t *ctx,
                                             p4est_connectivity_t *conn,
                                             p4est_topidx_t *pv2v);

#ifdef __cplusplus
}
#endif
#endif /* AMR_GRIDFILE_H */