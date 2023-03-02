/**
 * \file   amr_external_func.h
 * \author akirby
 */

#ifndef AMR_EXTERNAL_FUNC_H
#define AMR_EXTERNAL_FUNC_H

/* header files */
#include "precision_types.h"
#include "amr_utilities.h"
#include "amr_p4est_utilities.h"
#include "amr_geometry.h"

/* system header files */
#ifndef DOXYGEN_IGNORE
#  include <mpi.h>
#  include <vector>
#endif

#ifndef P4_TO_P8
#  include <p4est_extended.h>
#  include <p4est_nodes.h>
#  include <p4est_lnodes.h>
#else
#  include <p8est_extended.h>
#  include <p8est_nodes.h>
#  include <p8est_lnodes.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Initializes external wake3d storage and sets up Tioga data
 *
 * @param [in,out] p4est    p4est tree data structure
 * @param [in,out] ctx      simulation context data
 */
void external_func_setup_wake3d(p4est_t *p4est,ctx_t *ctx);

/** Deallocates external mpi and wake3d data storage
 *
 * @param [in,out] ctx    ctx_t data structure
 */
void external_func_deallocate_geom(ctx_t *ctx);

#ifdef __cplusplus
}
#endif
#endif /* AMR_EXTERNAL_FUNC_H */