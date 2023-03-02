/**
 * \file    amr_transfer_data.c
 * \ingroup amr_group
 * \author  akirby
 *
 * \brief   Communication of solution and geometry data for repartitioning.
 *          Step 1: Check if quad came from another mpi rank:
 *                  Tell old mpi rank to send that quad to new mpi rank (mpi_put)
 *          Step 2: After finishing mpi_put, check if this rank needs to send old
 *                  quad data. If so, grab the old quad data and send it.
 *          Step 3: After all the mpi sends/recvs, unpack the data.
 */

/* header files */
#include "amr_transfer_data.h"
#include <sc_search.h>

/** Given target, find index p such that gfq[p] <= target < gfq[p + 1].
 * \param [in] nmemb    Number of entries in array MINUS ONE.
 */
static int
p4est_bsearch_partition(p4est_gloidx_t target,
                        const p4est_gloidx_t *gfq, int nmemb){
  size_t              res;

  P4EST_ASSERT (nmemb > 0);
  P4EST_ASSERT (gfq[0] <= target);
  P4EST_ASSERT (target < gfq[nmemb]);

  res = sc_bsearch_range (&target, gfq, (size_t) nmemb,
                          sizeof (p4est_gloidx_t), p4est_gloidx_compare);
  P4EST_ASSERT (res < (size_t) nmemb);

  return (int) res;
}

transfer_info_t transfer_data_save_prepartition(p4est_t *p4est){
    p4est_topidx_t first_local_tree = p4est->first_local_tree;
    p4est_topidx_t last_local_tree = p4est->last_local_tree;
    sc_array_t *trees = p4est->trees;

    p4est_topidx_t t;
    p4est_locidx_t si;
    p4est_locidx_t n_quads;

    sc_array_t *quadrants;
    p4est_quadrant_t *quad;
    p4est_tree_t *tree;
    quad_data_t *data;
    int count;

    /* ============================= */
    /* Save data before partitioning */
    /* ============================= */
    transfer_info_t pre_info;

    pre_info.local_num_quadrants = p4est->local_num_quadrants;
    pre_info.pre_gfq = (p4est_gloidx_t *) malloc((p4est->mpisize+1)*sizeof(p4est_gloidx_t));
    pre_info.pre_data = (data_info_t *) malloc(p4est->local_num_quadrants * sizeof(data_info_t));

    /* save quad partition map */
    memcpy(pre_info.pre_gfq,
           p4est->global_first_quadrant,
           (p4est->mpisize+1)*sizeof(p4est_gloidx_t));

    /* save data sizes and addresses */
    count = 0;
    for (t = first_local_tree; t <= last_local_tree; ++t) {
        tree = p4est_tree_array_index (trees, t);
        quadrants = &(tree->quadrants);
        n_quads = (p4est_locidx_t) quadrants->elem_count;

        for (si = 0; si < n_quads; ++si) {
            quad = p4est_quadrant_array_index(quadrants, si);
            data = (quad_data_t *) quad->p.user_data;

            pre_info.pre_data[count].soln_size = data->soln_size;
            pre_info.pre_data[count].geom_size = data->geom_size;

            pre_info.pre_data[count].soln_address = data->soln;
            pre_info.pre_data[count].geom_address = data->geom;
            count++;
        }
    }
    return pre_info;
}

void transfer_data_destroy_prepartition(transfer_info_t *pre_info){
    if(pre_info->pre_gfq) {free(pre_info->pre_gfq);} pre_info->pre_gfq = NULL;
    if(pre_info->pre_data) {free(pre_info->pre_data);} pre_info->pre_data = NULL;
    pre_info->local_num_quadrants = -1;
}

Real transfer_data(p4est_t *p4est,transfer_info_t *pre_info){
    p4est_topidx_t first_local_tree = p4est->first_local_tree;
    p4est_topidx_t last_local_tree = p4est->last_local_tree;
    p4est_gloidx_t *post_gfq = p4est->global_first_quadrant;
    p4est_gloidx_t *pre_gfq = pre_info->pre_gfq;
    sc_array_t *trees = p4est->trees;

    p4est_topidx_t t;
    p4est_locidx_t si;
    p4est_locidx_t n_quads;

    sc_array_t *quadrants;
    p4est_quadrant_t *quad;
    p4est_tree_t *tree;
    quad_data_t *data;

    p4est_gloidx_t pre_begin, pre_end;
    p4est_gloidx_t post_begin, post_end;
    p4est_gloidx_t gbegin, gend;

    int q, num_senders, num_receivers;
    int first_sender, last_sender;
    int first_receiver, last_receiver;

    quad_data_t **user_data_ptrs = NULL, **ud;
    char **recv_soln_data = NULL;
    char **recv_geom_data = NULL;
    char **send_soln_data = NULL;
    char **send_geom_data = NULL;
    int *recv_process_map = NULL, *rpm;

    MPI_Request *request = NULL;

    int recv_proc_count = 0;
    int send_proc_count = 0;
    int nrequest = 0;
    int nrecv_proc;
    int nrecv,nsend;

    Real t1,t2;
    int i;

    /* =============================================================== */
    /* We follow p4est_transfer_fixed to determine which MPI processes */
    /* communicate after repartitioning the quadrants for regridding.  */
    /* =============================================================== */
    // See: github.com/cburstedde/p4est/blob/master/src/p4est_communication.c#L982

    /* start timer */
    t1 = amr_utilities_timer();

    /* =============================== */
    /* Step 0: Allocate Process Arrays */
    /* =============================== */
    user_data_ptrs = ud = (quad_data_t **) malloc(p4est->local_num_quadrants*sizeof(quad_data_t *));
    request = (MPI_Request *) malloc(sizeof(MPI_Request)*
                                     (2*p4est->local_num_quadrants +
                                      2*pre_info->local_num_quadrants));

    recv_geom_data = (char **) malloc(p4est->local_num_quadrants*sizeof(char *));
    send_geom_data = (char **) malloc(pre_info->local_num_quadrants*sizeof(char *));

    recv_process_map = rpm = (int *) malloc(2*p4est->local_num_quadrants*sizeof(int));

    /* ======================================== */
    /* Step 1: linear index user_data addresses */
    /* ======================================== */
    for (t = first_local_tree; t <= last_local_tree; ++t) {
        tree = p4est_tree_array_index (trees, t);
        quadrants = &(tree->quadrants);
        n_quads = (p4est_locidx_t) quadrants->elem_count;

        for (si = 0; si < n_quads; ++si) {
            quad = p4est_quadrant_array_index(quadrants, si);
            *ud++ = (quad_data_t *) quad->p.user_data;
        }
    }

    /* ======================================== */
    /* Step 2: grab local partition information */
    /* ======================================== */
    pre_begin = pre_gfq[p4est->mpirank];   // first global quad id on this rank before partitioning
    pre_end   = pre_gfq[p4est->mpirank+1]; // first global quad id on rank+1 before partitioning

    post_begin = post_gfq[p4est->mpirank];   // first global quad id on this rank after partitioning
    post_end   = post_gfq[p4est->mpirank+1]; // first global quad id on rank+1 after partitioning
#if 0
    if (p4est->mpirank == 0) {
        for (i = 0; i < p4est->mpisize; i++) {
            DGOUT(ctx->log_io,"Rank %d: Pre[%3ld] Post[%3ld]\n",i,pre_gfq[i],post_gfq[i]);
        }
        DGOUT(ctx->log_io,"        Pre[%3ld] Post[%3ld]\n",
              pre_gfq[p4est->mpisize],post_gfq[p4est->mpisize]);
    }
#endif

    /* =============================================== */
    /* Step 3: determine MPI processes to receive from */
    /* =============================================== */
    recv_proc_count = 0;
    if (post_begin < post_end) {

        /* our process as the receiver is not empty */
        first_sender = p4est_bsearch_partition(post_begin, pre_gfq, p4est->mpisize);
        last_sender  = p4est_bsearch_partition(post_end-1, pre_gfq, p4est->mpisize);
        num_senders  = last_sender - first_sender + 1;

        /* go through sender processes and post receive calls */
        gend = post_begin;
        int local_quad_id_begin = 0;
        for (q = first_sender; q <= last_sender; ++q) {

            /* prepare positions for the sender process q */
            gbegin = gend;
            gend = pre_gfq[q + 1];
            if(gend > post_end) gend = post_end;

            /* choose how to treat the sender process */
            if (gbegin == gend) {
                /* the sender process is empty; we need no message */
                P4EST_ASSERT (first_sender < q && q < last_sender);
            } else {

                /* nonzero message from this sender */
                int recv_nquad = gend - gbegin;
                if (q != p4est->mpirank) {
                    *rpm++ = local_quad_id_begin; // local quad id to start filling
                    *rpm++ = recv_nquad;          // nquads to recv

                    /* count bytes and allocate memory */
                    size_t geom_mem = 0;
                    for (nrecv = 0; nrecv < recv_nquad; ++nrecv) {
                        data = user_data_ptrs[local_quad_id_begin + nrecv];

                        geom_mem += data->geom_size;
                    }
                    recv_geom_data[recv_proc_count] = (char *) malloc(geom_mem);
#if 0
                    DGOUT(stdout,
                          "[Rank: %d] GFQ: Receiving: "
                          "number of quadrants=%d from rank %d (Start Quad id %d)"
                          " geom_size: %ld\n",
                          p4est->mpirank,recv_nquad,q,local_quad_id_begin,
                          geom_mem);
#endif
                    /* recv geometry data from old process */
                    MPI_Irecv(recv_geom_data[recv_proc_count], // recv buffer
                              geom_mem,                        // number of bytes
                              MPI_BYTE,                        // data type
                              q,                               // from rank
                              1,                               // geom tag
                              p4est->mpicomm,
                              &request[nrequest++]);

                    recv_proc_count++;
                }
                local_quad_id_begin += recv_nquad;
            }
        }
    }

    /* ========================================== */
    /* Step 4: determine MPI processes to send to */
    /* ========================================== */
    if (pre_begin < pre_end) {

        /* our process as the sender is not empty */
        first_receiver = p4est_bsearch_partition(pre_begin, post_gfq, p4est->mpisize);
        last_receiver  = p4est_bsearch_partition(pre_end-1, post_gfq, p4est->mpisize);
        num_receivers  = last_receiver - first_receiver + 1;

        /* go through receiver processes and post send calls */
        gend = pre_begin;
        int local_quad_id_begin = 0;
        for (q = first_receiver; q <= last_receiver; ++q) {
            /* prepare positions for the receiver process q */
            gbegin = gend;
            gend = post_gfq[q + 1];
            if(gend > pre_end) gend = pre_end;

            /* choose how to treat the receiver process */
            if (gbegin == gend) {
                /* the receiver process is empty; we need no message */
                P4EST_ASSERT (first_receiver < q && q < last_receiver);
            } else {

                /* nonzero message for this receiver */
                int send_nquad = gend - gbegin;
                if (q != p4est->mpirank) {
                    data_info_t *pre_data = &pre_info->pre_data[local_quad_id_begin];

                    /* count bytes and allocate memory */
                    size_t geom_mem = 0;
                    for (nsend = 0; nsend < send_nquad; ++nsend) {
                        geom_mem += pre_data[nsend].geom_size;
                    }

                    send_geom_data[send_proc_count] = (char *) malloc(geom_mem);

                    /* convenience pointers */
                    char *send_geom = send_geom_data[send_proc_count];

                    /* pack geometry data into contiguous messages */
                    geom_mem = 0;
                    for (nsend = 0; nsend < send_nquad; ++nsend) {
                        memcpy(&send_geom[geom_mem],
                               pre_data[nsend].geom_address,
                               pre_data[nsend].geom_size);
                        geom_mem += pre_data[nsend].geom_size;

                        /* free old memory */
                        if (pre_data[nsend].geom_address) {
                            free(pre_data[nsend].geom_address);
                            pre_data[nsend].geom_address = NULL;
                        }
                    }
#if 0
                    DGOUT(stdout,
                          "[Rank: %d] GFQ: Sending: "
                          "number of quadrants=%d to rank %d (Start Quad id %d)"
                          " geom_size: %ld\n",
                          p4est->mpirank,send_nquad,q,local_quad_id_begin,
                          geom_mem);
#endif
                    /* send geometry data to new process */
                    MPI_Isend(send_geom_data[send_proc_count], // send buffer
                              geom_mem,                        // number of bytes
                              MPI_BYTE,                        // data type
                              q,                               // to rank
                              1,                               // geom tag
                              p4est->mpicomm,
                              &request[nrequest++]);

                    send_proc_count++;
                }
                local_quad_id_begin += send_nquad;
            }
        }
    }

    /* ================ */
    /* Step 5: wait all */
    /* ================ */
    MPI_Waitall(nrequest,request,MPI_STATUSES_IGNORE);
    nrequest = 0;

    /* =================== */
    /* Step 6: unpack data */
    /* =================== */
    for (nrecv_proc = 0; nrecv_proc < recv_proc_count; ++nrecv_proc) {
        int local_quad_id_begin = recv_process_map[2*nrecv_proc+0];
        int recv_nquad          = recv_process_map[2*nrecv_proc+1];

        /* convenience pointers */
        char *recv_geom = recv_geom_data[nrecv_proc];

        size_t geom_mem = 0;
        for (nrecv = 0; nrecv < recv_nquad; ++nrecv) {
            data = user_data_ptrs[local_quad_id_begin + nrecv];
#if 0
            DGOUT(ctx->log_io,
                  "{Rank %d} UnPack Quad[%d]: "
                  "from ind %ld, geom size %d\n",
                  p4est->mpirank,
                  local_quad_id_begin + nrecv,
                  data->geom_size);
#endif
            data->geom = (Real *) malloc(data->geom_size);
            memcpy(data->geom,&recv_geom[geom_mem],data->geom_size);
            geom_mem += data->geom_size;
        }
    }

    /* =============== */
    /* Step 7: Cleanup */
    /* =============== */
    if(recv_process_map) {free(recv_process_map);} recv_process_map = NULL;
    if(user_data_ptrs) {free(user_data_ptrs);} user_data_ptrs = NULL;
    if(request) {free(request);} request = NULL;

    for (i = 0; i < recv_proc_count; ++i) {
        free(recv_geom_data[i]); recv_geom_data[i] = NULL;
    }
    if(recv_geom_data) {free(recv_geom_data);} recv_geom_data = NULL;

    for (i = 0; i < send_proc_count; ++i) {
        free(send_geom_data[i]); send_geom_data[i] = NULL;
    }
    if(send_geom_data) {free(send_geom_data);} send_geom_data = NULL;

    /* ================== */
    /* Step 8: Stop Timer */
    /* ================== */
    t2 = amr_utilities_timer();
    return amr_utilities_mpireducemax_real(p4est->mpirank,p4est->mpicomm,t2-t1);
}