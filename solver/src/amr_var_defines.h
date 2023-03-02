/**
 * \file    amr_var_defines.h
 * \author  akirby
 *
 * \brief   Global variable data related to simulation static constants:
 *          ghost fringe size, number of sub-elements in each amr quadrant,
 *          number of fields (unknowns) per point, and external
 *          communication constants.
 */

#ifndef AMR_VAR_DEFINES_H
#define AMR_VAR_DEFINES_H

/* header files */
#include "main.h"

/*============================================================================*/
#ifndef P4_TO_P8
#  define DGVDIM(x)   ((x) * (x))
#  define DGSDIM(x)    (x)
#  define ARG3D(x)
#  define ARG3Dl(x)
#  define E3D(x)
#  define X3D(x)
#  define OPDIM(x2d,x3d) (x2d)
#  define IFCON3D(x)
#else
#  define DGVDIM(x)   ((x) * (x) * (x))
#  define DGSDIM(x)   ((x) * (x))
#  define ARG3D(x)      x,
#  define ARG3Dl(x)     ,x
#  define E3D(x)        x;
#  define X3D(x)        x
#  define OPDIM(x2d,x3d) (x3d)
#  define IFCON3D(x)    &&x
#endif

#ifdef DEBUG
#  define DEBUG_MESG(x)    x;
#else
#  define DEBUG_MESG(x)
#endif

/*============================================================================*/
#define BUFF_SIZE 1024          /**< Number of string characters allowed      */
#define BASE 1                  /**< Base 1 to shift all indices for Tioga    */

#ifndef P4_TO_P8
#  define DIM 2
#else
#  define DIM 3
#endif

#endif /* AMR_VAR_DEFINES_H */