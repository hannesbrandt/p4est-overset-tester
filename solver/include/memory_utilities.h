/**
 * \file    memory_utilities.h
 * \ingroup cartdg_group
 * \author  akirby
 * Created on April 18, 2019, 3:47 PM
 */

#ifndef MEMORY_UTILITIES_H
#define MEMORY_UTILITIES_H

/* data header files */
#include "precision_types.h"

/* system header files */
#ifndef DOXYGEN_IGNORE
#  include <math.h>
#  include <stdio.h>
#  include <stdlib.h>
# ifndef __APPLE__
#    include <malloc.h>
# endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Utility function to get current memory usage
 *
 * @param [in] mpi_rank     this mpi rank
 * @param [in] display      flag indicator to display usage to console
 * @return                  current memory usage
 */
static inline
Real memory_usage(int mpi_rank,int display){
#ifndef __APPLE__
    /* get malloc info structure */
    struct mallinfo my_mallinfo = mallinfo();

    /*total memory reserved by the system for malloc currently */
    Real reserved_mem = my_mallinfo.arena;

    /* get all the memory currently allocated to user by malloc, etc. */
    Real used_mem = my_mallinfo.hblkhd
                  + my_mallinfo.usmblks
                  + my_mallinfo.uordblks;

    /* get memory not currently allocated to user but malloc controls */
    Real free_mem = my_mallinfo.fsmblks
                  + my_mallinfo.fordblks;

    /* get number of items currently allocated */
    /* double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks; */

    /* Print out concise malloc info line */
    if (display) {
        DGOUT(stdout,"Rank[%d]: %f MB(%.0f) malloc: %f MB reserved (%.0f unused)\n",
              mpi_rank,
              used_mem / (1024.0 * 1024.0),
              used_mem,
              reserved_mem / (1024.0 * 1024.0),
              free_mem);
    }
    return used_mem;
#else
    return 0.0;
#endif
}

static inline
Real memory_state(Real mem1,Real mem2,int mpi_rank,int display){
    Real mem_diff= mem2 - mem1;

    if (display) {
        DGOUT(stdout,"\033[1;32m[Memory] Rank[%d]: %f MB(%.0f) difference \033[0m\n",
              mpi_rank,
              mem_diff / (1024.0 * 1024.0),
              mem_diff);
    }
    return mem_diff;
}

#ifdef __cplusplus
}
#endif
#endif /* MEMORY_UTILITIES_H */