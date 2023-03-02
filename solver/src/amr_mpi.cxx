/** MPI wrapper functions for AMR code module.
 *
 * \file    amr_mpi.cxx
 * \author  akirby
 *
 * \brief   Wrapper functions around mpi and p4est initialization and finalize.
 */

/* header files */
#include "amr_mpi.h"

void mpi_init_from_comm(MPI_Comm *comm,ctx_t *ctx,int log_threshold){
    (void) log_threshold;
    /* ============================================ */
    /* initialize MPI information from communicator */
    /* ============================================ */
    ctx->comm = *comm;

    /* ======================= */
    /* initialize sc and p4est */
    /* ======================= */
    sc_set_log_defaults(ctx->log_io,NULL,SC_LP_ERROR);
    sc_init(ctx->comm,1,1,NULL,SC_LP_ERROR);
    p4est_init(NULL,SC_LP_ERROR);

    /* ======================== */
    /* assign MPI rank and size */
    /* ======================== */
    MPI_Comm_rank(ctx->comm,&ctx->rank);
    MPI_Comm_size(ctx->comm,&ctx->nranks);
}

int mpi_finalize(){
    return MPI_Finalize();
}