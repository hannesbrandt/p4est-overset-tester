/**
 * \file   amr_regrid_solution.h
 * \author akirby
 */

#ifndef AMR_REGRID_SOLUTION_H
#define AMR_REGRID_SOLUTION_H

/* header files */
#include "amr_regrid.h"
#include "amr_utilities.h"
#include "amr_transfer_data.h"
#include "dg4est_solver.hxx"
#include "precision_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Top-level regridding coordination function
 *
 * @param [in] ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in] initial          flag if initial time step
 * @return regrid wall-clock time
 */
Real regrid_solution(ctx_t *ctx,p4est_t *p4est,int initial);

/** Top-level regridding initial function:
 *      This regridding strategy is used to refine the grid
 *      up to the minimum AMR level when min_level>0.
 *      This is needed since the initialization of p4est
 *      must be carried out at level 0 for high-order
 *      unstructured grids.
 *
 * @param [in] ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in] initial          flag if initial time step
 * @return regrid wall-clock time
 */
Real regrid_2minlevel(ctx_t *ctx,p4est_t *p4est,int initial);

/** Top-level regridding strategy function #2
 *      This strategy does not allow features to be tagged above
 *      the maximum level allowed for the user pdegree_max
 *
 * @param [in] ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in] initial          flag if initial time step
 * @param [in] ref_max_level    maximum level to allow feature refinement
 * @return regrid wall-clock time
 */
Real regrid_callback(ctx_t *ctx,p4est_t *p4est,int initial,int ref_max_level);

#ifdef __cplusplus
}
#endif
#endif /* AMR_REGRID_SOLUTION_H */