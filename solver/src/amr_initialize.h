/**
 * \file   amr_initialize.h
 * \author akirby
 */

#ifndef AMR_INITIALIZE_H
#define AMR_INITIALIZE_H

/* header files */
#include "precision_types.h"
#include "dg4est_solver.hxx"
#include "amr_var_tag.h"
#include "amr_var_quad.h"
#include "amr_var_defines.h"
#include "amr_vtk.h"
#include "main.h"
#include "amr_utilities.h"
#include "amr_regrid_solution.h"
#include "amr_p4est_utilities.h"
#include "amr_gridfile.h"
#include "amr_containers.h"
#include "amr_geometry.h"


#ifndef DOXYGEN_IGNORE
#  include <math.h>
#  include <stdio.h>
#  include <stdlib.h>
#  include <string.h>
#endif

#ifndef P4_TO_P8
#  include <p4est_extended.h>
#else
#  include <p8est_extended.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Initializes the AMR code module from input files
 *
 * Builds output directories, reads input files, builds the mesh, and
 * initializes the solution data.
 *
 * @param [in]  filename    input file name
 * @param [out] dg4est      simulation data
 * @param [out] p4est       p4est tree structure
 * @param [out] conn        p4est tree connectivity
 */
void initialize_amr_file(char *filename,dg4est_t *dg4est,
                         p4est_t **p4est,p4est_connectivity_t **conn);

/** Reads the input file from the command line.
 *
 * @param [in]  filename        input file name
 * @param [out] dg4est          simulation data
 * @param [in]  noinput         flag if input file was found
 * @param [in]  reread_inputs   flag to reread dynamic inputs
 */
void initialize_inputs_from_file(char *filename,dg4est_t *dg4est,
                                 char noinput,char reread_inputs);

/** Defaults the context variables of the simulation.
 *
 * @param [out] ctx simulation context data
 */
void initialize_default_inputs(ctx_t *ctx);

/** Builds new grid or reads in restart file. Initializes external solver data.
 *
 * @param [in]  dg4est          simulation data
 * @param [out] p4est           p4est tree structure
 * @param [out] conn            p4est tree connectivity
 * @param [in]  regrid_now      flag to regrid
 * @param [in]  partition_now   flag to repartition
 * @param [in]  visualize_now   flag to visualize solution
 */
void initialize_grid_and_solver(dg4est_t *dg4est,
                                p4est_t **p4est,
                                p4est_connectivity_t **conn,
                                char regrid_now,char partition_now,char visualize_now);

/** Build grid constants from input file.
 *
 * @param [in,out] ctx   simulation context data
 */
void initialize_grid_constants(ctx_t *ctx);

/** Builds new Cartesian grid given user dimensions.
 *
 * @param [in,out] ctx   simulation context data
 * @param [out] conn    p4est tree connectivity
 */
void initialize_grid(ctx_t *ctx,
                     p4est_connectivity_t **conn);

/** Initialize p4est forset of trees.
 *
 * @param [in]     dg4est          simulation data
 * @param [in,out] ctx             simulation context data
 * @param [out]    conn            p4est tree connectivity
 * @return         p4est           p4est datatype
 */
p4est_t *initialize_p4est(dg4est_t *dg4est,ctx_t *ctx,p4est_connectivity_t *conn);

/** Calculates the total grid volume.
 *
 * @param [in,out] ctx      simulation context data
 * @param [in]     p4est    p4est tree structure
 */
void initialize_grid_volume(ctx_t *ctx,p4est_t *p4est);

/** Callback function to initialize solution data in each amr quadrant/octant.
 *
 * @param [in]  p4est       p4est tree structure
 * @param [in]  which_tree  tree in the forest of octrees
 * @param [out] quad        quadrant solution data
 */
void initialize_quadrant_data(p4est_t *p4est,
                              p4est_topidx_t which_tree,
                              p4est_quadrant_t *quad);

/** Callback function to initialize boundary conditions in each amr quadrant/octant.
 *
 * @param [in]  p4est       p4est tree structure
 * @param [in]  which_tree  tree in the forest of octrees
 * @param [out] quad        quadrant solution data
 */
void initialize_quadrant_bc_data(p4est_iter_face_info_t *info,void *user_data);


/** Callback function to calculate quadrant's volume
 *
 * @param [in]    info          p4est volume iteration information
 * @param [inout] user_data     local domain grid volume accumulator
 */
void domain_volume_callback(p4est_iter_volume_info_t *info,void *user_data);

#ifdef __cplusplus
}
#endif
#endif /* AMR_INITIALIZE_H */
