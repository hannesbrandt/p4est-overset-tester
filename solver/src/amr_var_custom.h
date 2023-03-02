/*
 * \file    amr_var_custom.h
 * \ingroup amr_group
 * \author  akirby
 */

#ifndef AMR_VAR_CUSTOM_H
#define AMR_VAR_CUSTOM_H

/* header files */
#include "precision_types.h"
#include "amr_var_quad.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int soln_size;
    int geom_size;
    Real *soln_address;
    Real *geom_address;
}
data_info_t;

#ifdef __cplusplus
}
#endif
#endif /* AMR_VAR_CUSTOM_H */