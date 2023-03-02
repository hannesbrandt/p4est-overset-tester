/**
 * \file   amr_regrid.h
 * \author akirby
 */

#ifndef AMR_REGRID_H
#define AMR_REGRID_H

/* header files */
#include "amr_tag.h"
#include "amr_utilities.h"
#include "amr_initialize.h"
#include "amr_p4est_utilities.h"
#include "dg4est_solver.hxx"
#include "precision_types.h"
#include "amr_transfer_data.h"

#ifndef P4_TO_P8
#  include "p4est_search.h"
#else
#  include "p8est_search.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Regrid coarsening function
 *
 * @param [in,out] p4est        p4est forest structure
 * @param [in]  initial         flag if initial time step
 * @return coarsen wall-clock time
 */
Real regrid_coarsen(p4est_t *p4est,int initial);

/** Mesh 2:1 quadrant balancing
 *
 * @param [in,out] p4est        p4est forest structure
 * @return time to balance mesh
 */
Real regrid_balance(p4est_t *p4est);

/** Regridding function for refining features
 *
 * @param [in] ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @return time to regrid features
 */
Real regrid_feature(ctx_t *ctx,p4est_t *p4est,int initial,int ref_max_level);

/** Regridding function for refining from level 0 to mim AMR level
 *
 * @param [in] ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @return time to regrid features
 */
Real regrid_level(ctx_t *ctx,p4est_t *p4est);

/** Mesh 2:1 quadrant balancing and mpi partitioning function
 *
 * @param [in,out] p4est            p4est forest structure
 * @return time to partition mesh
 */
Real regrid_partition(p4est_t *p4est);

/** MPI partition weighting function per quadrant (external)
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     p4est quadrant that needs to be weighted
 * @return partition weight for this quadrant
 */
int  regrid_load(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant);

/** Solution replacement function when coarsening or refining (external)
 *
 * @param [in] p4est            p4est forest structure
 * @param [in] which_tree       tree id in p4est forest
 * @param [in] num_outgoing     number of quadrants being replaced
 * @param [out] outgoing        quadrants being replaced
 * @param [in] num_incoming     number of quadrants being introduced
 * @param [in] incoming         quadrant being introduced
 */
void regrid_replace_quads(p4est_t *p4est,p4est_topidx_t which_tree,
                          int num_outgoing,p4est_quadrant_t *outgoing[],
                          int num_incoming,p4est_quadrant_t *incoming[]);

#ifdef __cplusplus
}
#endif
#endif /* AMR_REGRID_H */