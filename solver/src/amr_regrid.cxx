/**
 * \file    amr_regrid.cxx
 * \author  akirby
 *
 * \brief   Regridding functions for the AMR code module.
 */

/* header files */
#include "amr_regrid.h"

Real regrid_feature(ctx_t *ctx,p4est_t *p4est,int initial,int ref_max_level){
    grid_t *grid = &ctx->d_grid;

    int recursive;
    int l;
    Real t1,t2;

    if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
        DGOUT(ctx->log_io,
              "[regrid] Into feature regridding...: max level global: %d\n",
              grid->max_level_global);
    }

    t1 = amr_utilities_timer(); /* start timer */
    if (grid->max_level_global > 0) {
        if (initial) {
            recursive = 0;

            /* build one level at a time */
            for (l = grid->min_level+1; l <= grid->max_level_global; ++l) {
                tag_reset(p4est,no_tag,0);

                /* refine tagged quadrants */
                p4est_refine_ext(p4est,recursive,l,tag_feature,
                                 NULL,regrid_replace_quads);

                regrid_balance(p4est);
                regrid_partition(p4est);
            }
            tag_reset_max_level(p4est,feature_tag);
        } else {
            recursive = 1;
            tag_reset(p4est,no_tag,1);
            p4est_refine_ext(p4est,recursive,grid->max_level_global,tag_feature,
                             NULL,regrid_replace_quads);
        }
    }
    t2 = amr_utilities_timer(); /* stop timer */

    if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
        DGOUT(ctx->log_io,"[regrid] Done feature regridding...\n");
    }
    return amr_utilities_mpireducemax_real(p4est->mpirank,p4est->mpicomm,t2-t1);
}

Real regrid_level(ctx_t *ctx,p4est_t *p4est){
    grid_t *grid = &ctx->d_grid;

    int initial = 1; // This should only be called during initialization
    int recursive;
    int l;
    Real t1,t2;

    if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
        printf("[regrid] Into level-based regridding...\n");
    }

    t1 = amr_utilities_timer(); /* start timer */
    if (grid->min_level > 0) {
        if (initial) {
            recursive = 0;

            /* build one level at a time to reach min_level*/
            for (l = 1; l <= grid->min_level; ++l) {
                tag_reset(p4est,no_tag,0);

                /* refine tagged quadrants */
                p4est_refine_ext(p4est,recursive,l,tag_level,
                                 NULL,regrid_replace_quads);

                regrid_balance(p4est);
                regrid_partition(p4est);
            }
        } else {
            recursive = 1;
            p4est_refine_ext(p4est,recursive,grid->min_level,tag_level,
                             NULL,regrid_replace_quads);
        }
    }
    t2 = amr_utilities_timer(); /* stop timer */

    if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
        printf("[regrid] Done level-based regridding...\n");
    }
    return amr_utilities_mpireducemax_real(p4est->mpirank,p4est->mpicomm,t2-t1);
}

Real regrid_coarsen(p4est_t *p4est,int initial){
    int recursive = 0;
    int callbackorphans = 0;
    Real t1,t2;

    /* coarsen if not tagged */
    t1 = amr_utilities_timer(); /* start timer */
    if (!initial) {
        p4est_coarsen_ext(p4est,recursive,callbackorphans,tag_coarsen,NULL,
                          regrid_replace_quads);
    }
    t2 = amr_utilities_timer(); /* stop timer */
    return amr_utilities_mpireducemax_real(p4est->mpirank,p4est->mpicomm,t2-t1);
}

Real regrid_balance(p4est_t *p4est){
    dg4est_t *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t       *ctx = &dg4est->d_ctx;
    grid_t      *grd = &ctx->d_grid;
    Real t1,t2;

    /* Do not use initialize quadrant data --
     *   The coarsen operator destroys the 2:1 balance quadrants because they
     *   are not tagged. When re-balance occurs, regrid_replace_quads is used.
     *   Therefore the grid changes even if we don't change the solution.
     */

    /* balance mesh 2:1 element ratio */
    t1 = amr_utilities_timer(); /* start timer */
    if (grd->max_level_global != grd->min_level) {
        p4est_balance_ext(p4est,P4EST_CONNECT_FULL,NULL,regrid_replace_quads);
    }
    t2 = amr_utilities_timer(); /* stop timer */
    return amr_utilities_mpireducemax_real(p4est->mpirank,p4est->mpicomm,t2-t1);
}

Real regrid_partition(p4est_t *p4est){
    int partition_for_coarsening = 1;
    transfer_info_t pre_info;
    Real t1,t2;

    t1 = amr_utilities_timer(); /* start timer */

    /* save pre-partition information */
    pre_info = transfer_data_save_prepartition(p4est);


    /* partition mesh based on weighting */
    p4est_partition(p4est,partition_for_coarsening,regrid_load);

    /* move custom data */
    transfer_data(p4est,&pre_info);
    transfer_data_destroy_prepartition(&pre_info);

    t2 = amr_utilities_timer(); /* stop timer */
    return amr_utilities_mpireducemax_real(p4est->mpirank,p4est->mpicomm,t2-t1);
}

int regrid_load(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant){
    (void) p4est;
    (void) which_tree;
    return 1;
}

void quad_to_user_data(quad_data_t *quad_data,quad_user_data_t *user_data){
    int i;
    user_data->qdegree   = quad_data->qdegree;
    user_data->geom_size = quad_data->geom_size;
    user_data->quad_geom = quad_data->geom;
    user_data->type      = quad_data->type;
    for(i = 0; i < 6; i++) user_data->bc[i] = quad_data->bc[i];
}

void user_to_quad_data(quad_user_data_t *user_data,quad_data_t *quad_data){
    int i;
    quad_data->qdegree   = user_data->qdegree;
    quad_data->geom_size = user_data->geom_size;
    quad_data->geom      = user_data->quad_geom;
    quad_data->type      = user_data->type;
    for(i = 0; i < 6; i++) quad_data->bc[i] = user_data->bc[i];
}

void regrid_replace_quads(p4est_t *p4est,p4est_topidx_t which_tree,
                          int num_outgoing,p4est_quadrant_t *outgoing[],
                          int num_incoming,p4est_quadrant_t *incoming[]){
    char hole_mask_flag;
    int tag_flag;
    int i;

    dg4est_t  *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t        *ctx = &dg4est->d_ctx;
    quad_data_t *parent_data;
    quad_data_t *child_data;

    quad_properties_t coarse_quad_prop;
    quad_user_data_t parent_user_data;
    quad_user_data_t children_user_data[P4EST_CHILDREN];

    /* assign coarse element properties */
    coarse_quad_prop.dim          = ctx->d_simulation.dim;
    coarse_quad_prop.unstructured = ctx->d_simulation.unstructured_flag;

    tag_flag = 0;
    hole_mask_flag = 1;
    if (num_outgoing > 1) {

        /* coarsen */
        parent_data = (quad_data_t *) incoming[0]->p.user_data;

        /* assign data pointers */
        for (i = 0; i < P4EST_CHILDREN; ++i) {
            child_data = (quad_data_t *) outgoing[i]->p.user_data;
            quad_to_user_data(child_data,&children_user_data[i]);

            if(child_data->tag) tag_flag = 1;
            if(child_data->iblank == hole_mask) hole_mask_flag = 0;
        }

        /* convert children quadrants to pmax_children p-degree */
        for (i = 0; i < P4EST_CHILDREN; ++i) {
            child_data = (quad_data_t *) outgoing[i]->p.user_data;
        }

        /* set defaults */
        parent_data->level    = incoming[0]->level;
        parent_data->tag      = (ctx->initial) ? no_tag:tag_flag;
        parent_data->type     = children_user_data[0].type;
        parent_data->iblank   = field_mask & hole_mask_flag;
        parent_data->qdegree  = children_user_data[0].qdegree;
        parent_data->geom_size = children_user_data[0].geom_size;
        parent_data->geom     = (Real *) malloc(parent_data->geom_size);

        /* set defaults bcs */
        for(i = 0; i < 6; i++) parent_data->bc[i] = NO_PATCH;

        /* transfer bcs */
        /*                parent               children
        *                 y_hi(3)
        *              ___________           ___________
        *             |           |         |  2  |  3  |
        *     x_lo(0) |           | x_hi(1) |_____|_____|
        *   ^ y       |           |         |  0  |  1  |
        *   |         |___________|         |_____|_____|
        *   |             y_lo(2)
        *    ---> x
        */
        parent_data->bc[0] = children_user_data[0].bc[0]; // x_lo
        parent_data->bc[1] = children_user_data[1].bc[1]; // x_hi
        parent_data->bc[2] = children_user_data[0].bc[2]; // y_lo
        parent_data->bc[3] = children_user_data[2].bc[3]; // y_hi
    X3D(parent_data->bc[4] = children_user_data[0].bc[4];)// z_lo
    X3D(parent_data->bc[5] = children_user_data[4].bc[5];)// z_hi

        /* assign data pointers */
        quad_to_user_data(parent_data,&parent_user_data);

        /* get quad geometry: xyz */
        p4est_utilities_quad_coordinates(p4est,which_tree,outgoing[i],parent_data->geom,0);

        /* free children data */
        for (i = 0; i < (1 << DIM); ++i) {
            child_data = (quad_data_t *) outgoing[i]->p.user_data;
            if(child_data->geom) {free(child_data->geom);} child_data->geom = NULL;
        }

    } else
    if (num_incoming > 1) {

        /* refine */
        parent_data = (quad_data_t *) outgoing[0]->p.user_data;

        /* assign data pointers */
        quad_to_user_data(parent_data,&parent_user_data);

        /* set defaults */
        for (i = 0; i < P4EST_CHILDREN; ++i) {
            child_data = (quad_data_t *) incoming[i]->p.user_data;

            child_data->level   = incoming[i]->level;
            child_data->tag     = (ctx->initial) ? no_tag:parent_data->tag;
            child_data->type    = parent_data->type;
            child_data->iblank  = parent_data->iblank;
            child_data->qdegree = parent_data->qdegree;
            child_data->geom_size =  parent_data->geom_size;
            child_data->geom = (Real *) malloc(child_data->geom_size);

            /* fill geometry */
            p4est_utilities_quad_coordinates(p4est,which_tree,incoming[i],child_data->geom,0);

            /* copy bcs */
            memcpy(child_data->bc,parent_data->bc,6*sizeof(char));
        }

        /* free parent data */
        if(parent_data->geom) {free(parent_data->geom);} parent_data->geom = NULL;
    }
}