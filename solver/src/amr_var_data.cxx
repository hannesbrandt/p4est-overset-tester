/**
 * \file    amr_var_data.c
 * \ingroup amr_group
 * \author  akirby
 */

/* header files */
#include "amr_var_data.h"

void ctx_default(ctx_t *ctx){
    ctx->version = VERSION;
    ctx->comm = MPI_COMM_WORLD;
    ctx->log_io = stdout;
    ctx->rank = 0;
    ctx->nranks = 1;
    ctx->initial = 0;
    ctx->visualize = 0;
    ctx->nvisualize = 0;
    ctx->high_order_viz = 0;
    ctx->high_order_viz_npts = 3;

    grid_default(&ctx->d_grid);
    wake3d_default(&ctx->d_wake3d);
    geometry_default(&ctx->d_geometry);
    amr_hbox_default(&ctx->d_amr_hboxes);
}

void geometry_default(geometry_t *geom){
    geom->npatch = 0;
    geom->nelem_face = 0;
    geom->nentity_vol = 0;
    geom->nentity_face = 0;
    geom->patch_name = NULL;
}

void grid_default(grid_t *grid){
    grid->nlevels = 1;
    grid->min_level = 0;
    grid->max_level = 0;
    grid->max_level_pmax = 0;
    grid->max_level_global = 0;
    grid->nelem[0] = 1;
    grid->nelem[1] = 1;
    grid->nelem[2] = 1;
    grid->periodic[0] = 0;
    grid->periodic[1] = 0;
    grid->periodic[2] = 0;
    grid->construct_grid = 1;
    grid->xlo[0] = 0.0;
    grid->xlo[1] = 0.0;
    grid->xlo[2] = 0.0;
    grid->xhi[0] = 1.0;
    grid->xhi[1] = 1.0;
    grid->xhi[2] = 1.0;
    grid->min_dx[0] = 1.0;
    grid->min_dx[1] = 1.0;
    grid->min_dx[2] = 1.0;
    grid->max_dx[0] = 1.0;
    grid->max_dx[1] = 1.0;
    grid->max_dx[2] = 1.0;
    grid->domain_volume = 0.0;
    grid->total_quads = 1;
    grid->total_dofs = 1;
}


void wake3d_default(wake3d_t *wake3d){
    wake3d->igbp_scale = 1.0;
    wake3d->igbp_scale_cf = 1.0;
    wake3d->wake3d_dt = 1.0;
    wake3d->nnodes = 0;
    wake3d->nobc = 0;
}

void amr_hbox_default(amr_hbox_t *hbox){
    hbox->amr_nhboxes = 0;
}

void ctx_shutdown(ctx_t *ctx){}
