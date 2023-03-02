/**
 * \file   amr_mpi.h
 * \author akirby
 */

#ifndef AMR_MPI_H
#define AMR_MPI_H

/* header files */
#include "precision_types.h"
#include "amr_var_custom.h"
#include "dg4est_solver.hxx"

#ifndef P4_TO_P8
#  include <p4est_bits.h>
#else
#  include <p8est_bits.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** MPI and p4est initialization function wrapper given a mpi communicator
 *
 * @param [in]    comm              mpi communicator
 * @param [inout] ctx               context data
 * @param [in]    log_threshold     logging threshold flag
 */
void mpi_init_from_comm(MPI_Comm *comm,ctx_t *ctx,int log_threshold);

/** MPI finalization function wrapper
 *
 * @return mpi finalize success
 */
int mpi_finalize();

#ifdef __cplusplus
}
#endif
#endif /* AMR_MPI_H */