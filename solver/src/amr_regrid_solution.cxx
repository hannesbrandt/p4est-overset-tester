/**
 * \file    amr_regrid_solution.c
 * \ingroup amr_group
 * \author  akirby
 *
 * \brief   Regridding functions for the AMR code module.
 */

/* header files */
#include "amr_regrid_solution.h"

Real regrid_solution(ctx_t *ctx,p4est_t *p4est,int initial){
    grid_t *grid = &ctx->d_grid;
    Real regrid_time;

    if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
        DGOUT(ctx->log_io,"\n[regrid] Into regridding...\n");
    }

    regrid_time = regrid_callback(ctx,p4est,initial,grid->max_level);

    if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
        DGOUT(ctx->log_io,"[regrid] Done regridding...\n");
    }
    return regrid_time;
}

Real regrid_2minlevel(ctx_t *ctx,p4est_t *p4est,int current_level){
    grid_t *grid = &ctx->d_grid;

    Real feature_regrid_time = 0.0;
    Real feature_spread_time = 0.0;
    Real buffer_regrid_time  = 0.0;
    Real point_regrid_time   = 0.0;
    Real coarsen_time        = 0.0;
    Real balance_time        = 0.0;
    Real pdegree_time        = 0.0;
    Real part_time           = 0.0;
    Real ship_time           = 0.0;
    Real max_time            = 0.0;
    Real level_regrid_time   = 0.0;
    Real t1,t2;

    if (current_level < grid->min_level) {
        if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
            DGOUT(ctx->log_io,"\n[regrid] Into regridding to min AMR level...\n");
        }

        t1 = amr_utilities_timer(); /* start timer */
            tag_reset(p4est,no_tag,0);
            level_regrid_time = regrid_level(ctx,p4est);
            balance_time = regrid_balance(p4est);
            part_time = regrid_partition(p4est);
        t2 = amr_utilities_timer(); /* stop timer*/
        max_time = amr_utilities_mpireducemax_real(ctx->rank,ctx->comm,t2-t1);

        p4est_utilities_mesh_stats(p4est);
        amr_utilities_write_amr_regrid(ctx,ctx->rank,max_time,
                                       balance_time,part_time,ship_time,
                                       point_regrid_time,feature_regrid_time,
                                       feature_spread_time,coarsen_time,buffer_regrid_time,
                                       pdegree_time,level_regrid_time,ctx->log_info);

        if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
            DGOUT(ctx->log_io,"[regrid] Done regridding...\n");
        }
        return max_time;
    }
    return 0.0;
}

Real regrid_callback(ctx_t *ctx,p4est_t *p4est,int initial,int ref_max_level){
    simulation_t *sim = &ctx->d_simulation;
    Real feature_regrid_time = 0.0;
    Real coarsen_time        = 0.0;
    Real balance_time        = 0.0;
    Real part_time           = 0.0;
    Real max_time            = 0.0;

    Real t1,t2;

    t1 = amr_utilities_timer(); /* start timer */
        tag_reset(p4est,no_tag,0);

        /* Step 2: feature tagging */
        if (sim->regrid_features) {
            feature_regrid_time = regrid_feature(ctx,p4est,initial,ref_max_level);
        }

        /* Step 4: coarsen quadrants not tagged as IGBP or features */
        if (sim->regrid_coarsen) {
            coarsen_time = regrid_coarsen(p4est,initial);
        }

        /* Step 5: 2:1 balance grid */
        balance_time = regrid_balance(p4est);

        /* Step 9: partition the mesh */
        part_time = regrid_partition(p4est);

    t2 = amr_utilities_timer(); /* stop timer*/
    max_time = amr_utilities_mpireducemax_real(ctx->rank,ctx->comm,t2-t1);

    p4est_utilities_mesh_stats(p4est);
    return max_time;
}
