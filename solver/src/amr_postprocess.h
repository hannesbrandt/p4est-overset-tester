/**
 * \file   amr_postprocess.h
 * \author akirby
 */

#ifndef AMR_POSTPROCESS_H
#define AMR_POSTPROCESS_H

/* header files */
#include "precision_types.h"
#include "amr_vtk.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Post processing control function.
 *
 * @param [in] ctx      simulation context data
 * @param [in] p4est    p4est tree data structure
 */
void postprocess(ctx_t *ctx,p4est_t **p4est);

#ifdef __cplusplus
}
#endif
#endif /* AMR_POSTPROCESS_H */
