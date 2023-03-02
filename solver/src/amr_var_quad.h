/**
 * \file    amr_var_quad.h
 * \ingroup amr_group
 * \author  akirby
 *
 * \brief   Quadrant solution data storage.
 *          This is the data that is stored, communicated, adapted, visualized,
 *          and checkpoint by p4est.
 */

#ifndef AMR_VAR_QUAD_H
#define AMR_VAR_QUAD_H

/* header files */
#include "precision_types.h"
#include "amr_var_defines.h"

#ifndef DOXYGEN_IGNORE
#  include <stddef.h>
#  include <stdint.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * quad_data_t contains per quadrant p4est related data.
 */
typedef struct {
    char bc[6];     /**< quadrant boundary condition flag (patch index) */
    char type;      /**< quadrant type (patch index) */
    int tag;        /**< grid refinement tag type indicator */
    int level;      /**< quadrant level */
    int iblank;     /**< iblank tag for Tioga overset */
    int qdegree;    /**< quadrant geometry degree */
    int soln_size;  /**< Size [bytes] of geom data */
    int geom_size;  /**< Size [bytes] of geom data */
    Real *soln;     /**< Quadrant solution data */
    Real *geom;     /**< Quadrant geometry data */
}
quad_data_t; /**< Internal data type for per quadratant p4est data */

#ifdef __cplusplus
}
#endif
#endif /* AMR_VAR_QUAD_H */