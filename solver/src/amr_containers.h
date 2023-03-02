/**
 * \file   amr_containers.h
 * \author mugolott
 */

#ifndef AMR_CONTAINERS_H
#define AMR_CONTAINERS_H

/* header files */
#include "amr_var_defines.h"
#include "precision_types.h"

#ifndef DOXYGEN_IGNORE
#  include <stdio.h>
#  include <stdlib.h>
#  include <stdint.h>
#  include <string.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Dynamically allocatable array of integers
 *
 *  Vector
 */
typedef struct {
    Uint *array;          // Vector of Uint
    Uint size;            // Size of vector based on entries added
    Uint capacity;        // Size of memeory allocated
    int  buffer;          // Increase capacity of buff_add when reached max capacity
}
ivector_t;

/*
 * Initialize dynamic array
 */
void ivector_init(ivector_t *vec, int initialSize);

/*
 * Push back entry
 */
void ivector_pushback(ivector_t *vec, Uint element);

/*
 * Trim capacity
 */
void ivector_trim(ivector_t *vec);

/*
 * Get value by ID
 */
Uint ivector_get_value(ivector_t *vec, int index);

/*
 * Dump
 */
void ivector_dump(ivector_t *vec);

/*
 * Free array
 */
void ivector_free(ivector_t *vec);

#ifdef __cplusplus
}
#endif
#endif /* AMR_GRIDFILE_H */