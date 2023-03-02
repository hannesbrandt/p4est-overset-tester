/**
 * \file:   amr_transfer_data.h
 * \author: akirby
 */

#ifndef AMR_TRANSFER_DATA_H
#define AMR_TRANSFER_DATA_H

/* header files */
#include "precision_types.h"
#include "amr_utilities.h"
#include "amr_initialize.h"
#include "amr_var_custom.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    p4est_locidx_t local_num_quadrants; /**< total quadrants on this mpi process */
    p4est_gloidx_t *pre_gfq;            /**< [mpisize+1]: global_first_quadrant p4est map */
    data_info_t *pre_data;              /**< [local_num_quadrants]: */
}
transfer_info_t;

/** Saves pre-partitioning info: global first quadrant map, solution/geometry sizes and addresses.
 *
 * @param  [in] p4est   p4est tree data structure
 * @return pre-partitioned information
 */
transfer_info_t transfer_data_save_prepartition(p4est_t *p4est);

/** Destroy data from transfer_data_save_prepartition.
 *
 * @param [in/out] pre_info     per-partitioned info to be destroyed
 */
void transfer_data_destroy_prepartition(transfer_info_t *pre_info);

/** Packs and communicates custom user data not stored in the p4est internals via All to All
 *
 * @param [in]  p4est       p4est tree data structure
 * @param [in]  pre_info    pre-partitioned data info from transfer_data_save_prepartition
 * @return max wall-clock time to perform data move
 */
Real transfer_data(p4est_t *p4est,transfer_info_t *pre_info);


#ifdef __cplusplus
}
#endif
#endif /* AMR_TRANSFER_DATA_H */