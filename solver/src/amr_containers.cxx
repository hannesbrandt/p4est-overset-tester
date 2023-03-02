/**
 * \file    amr_containers.cxx
 * \author  mugolott
 *
 * \brief Grid-file IO
 */

/* header files */
#include "amr_containers.h"

/*
 * Extend array
 */
void ivector_extend(ivector_t *vec, int new_capacity) {

    /* Check if the new_capacity is smaller than size */
    if (new_capacity < vec->size){
        printf("\x1B[1;31m[ amr ] Container capacity cannot be less than size! Implementation error. \x1B[0m\n");
        exit(EXIT_FAILURE);
    }
    /* Create a temo array */
    Uint *tmp_array = (Uint *) malloc(vec->size);
    memset(tmp_array,0,vec->size);

    /* Copy array to temp array*/
    memcpy(tmp_array,vec->array,vec->size);

    /* Reallocate vec->array */
    vec->array = (Uint *) realloc(vec->array, new_capacity * sizeof(Uint));

    /* Copy back to vec->array */
    memcpy(vec->array,tmp_array,vec->size);

    /* Free tmp_array */
    free(tmp_array);

    /* Change capacity */
    vec->capacity = new_capacity;
}

/*
 * Initialize dynamic array
 */
void ivector_init(ivector_t *vec, int initialSize) {
    vec->array = (Uint *) malloc(initialSize * sizeof(Uint));
    vec->size = 0;
    vec->capacity = initialSize;
    vec->buffer = 10;
}

/*
 * Push back entry
 */
void ivector_pushback(ivector_t *vec, Uint element) {
    if (vec->size == vec->capacity) {
        /* Extend array with buffer size */
        ivector_extend(vec,vec->capacity+vec->buffer);
    }
    vec->array[vec->size++] = element;
}

/*
 * Trim capacity
 */
void ivector_trim(ivector_t *vec) {
    /* Trim array with array size */
    ivector_extend(vec,vec->size);
}

/*
 * Get value by ID
 */
Uint ivector_get_value(ivector_t *vec, int index) {
    if (index > vec->size) {
        printf("\x1B[1;31m[ amr ] Index exceeds vector size! \x1B[0m\n");
        exit(EXIT_FAILURE);
    } else {
        return vec->array[index];
    }
}

/*
 * dump
 */
void ivector_dump(ivector_t *vec){
    Uint i;
    for ( i=0;i<vec->size;i++){
        printf("index %llu, value %llu \n",i,vec->array[i]);
    }
}

/*
 * Free array
 */
void ivector_free(ivector_t *vec) {
    free(vec->array);
    vec->array = NULL;
    vec->size = vec->capacity = 0;
}