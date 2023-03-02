/**
 * \file   amr_geometry.h
 * \author mugolott
 */

#ifndef AMR_GEOMETRY_H
#define AMR_GEOMETRY_H

/* header files */
#include "amr_var_defines.h"
#include "precision_types.h"
#include "dg4est_solver.hxx"

#ifndef DOXYGEN_IGNORE
#  include <math.h>
#  include <stdio.h>
#  include <stdlib.h>
#  include <stdint.h>
#  include <string.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Allocate geometry patch data
 *
 * @param [in] geom     geometry data
 */
void geometry_allocate_patch_data(geometry_t *geom);

/** Allocate geometry entity data
 *
 * @param [in] geom     geometry data
 */
void geometry_allocate_entity_data(geometry_t *geom);

/** Allocate geometry facenode data
 *
 * @param [in] geom     geometry data
 * @param num_faces     number of grid faces
 */
void geometry_allocate_facenode_data(geometry_t *geom,size_t num_faces);

/** Deallocate geometry_t
 *
 * @param [in]  ctx     simulation context data
 */
void geometry_deallocate_data(ctx_t *ctx);

/** Deallocate geometry_t face data
 *
 * @param [in]  ctx     simulation context data
 */
void geometry_deallocate_face_data(ctx_t *ctx);

/** Broadcast geometry data to all cores
 *
 * @param [in]  ctx     simulation context data
 */
void geometry_broadcast_data(ctx_t *ctx);
#ifdef __cplusplus
}
#endif
#endif /* AMR_GEOMETRY_H */