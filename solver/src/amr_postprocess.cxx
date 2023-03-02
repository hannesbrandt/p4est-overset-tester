/**
 * \file    amr_postprocess.cxx
 * \author  akirby
 *
 * \brief   Solution post processing controller.
 */

/* header files */
#include "amr_postprocess.h"

void postprocess(ctx_t *ctx,p4est_t **p4est){
    /* visualization solution */
    if (ctx->visualize) {
        ctx->nvisualize++;
            vtk_write_solution(ctx->d_simulation.visualization_path,*p4est,
                               1,
                               ctx->high_order_viz,
                               ctx->high_order_viz_npts);
    }
}
