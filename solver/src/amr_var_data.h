/**
 * \file    amr_var_data.h
 * \ingroup amr_group
 * \author  akirby
 *
 * \brief   Simulation context data related to other data types:
 *              mpi, grid, external solver, simulation, wake3d.
 */

#ifndef AMR_VAR_DATA_H
#define AMR_VAR_DATA_H

/* header files */
#include "amr_var_defines.h"
#include "dg4est_solver.hxx"
#include "precision_types.h"

#ifndef DOXYGEN_IGNORE
#  include <limits.h>
#  include <string.h>
#  include <stdlib.h>
#  include <stdio.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

void ctx_default(ctx_t *ctx);
void grid_default(grid_t *grid);
void wake3d_default(wake3d_t *wake3d);
void geometry_default(geometry_t *geom);
void amr_hbox_default(amr_hbox_t *hbox);
void ctx_shutdown(ctx_t *ctx);

#ifdef __cplusplus
}
#endif
#endif /* AMR_VAR_DATA_H */