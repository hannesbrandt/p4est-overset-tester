/**
 * \file    amr_tag.cxx
 * \ingroup amr_group
 * \author  akirby
 *
 * \brief   Tagging functions for regridding for the AMR code module.
 */

/* header files */
#include "amr_tag.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

int tag_reset(p4est_t *p4est,int tag_default,char keep_pt_tag){
    p4est_topidx_t t;
    p4est_topidx_t first_local_tree = p4est->first_local_tree;
    p4est_topidx_t last_local_tree = p4est->last_local_tree;
    sc_array_t *trees = p4est->trees;
    sc_array_t *quadrants;
    p4est_tree_t *tree;

    p4est_quadrant_t *quad;
    quad_data_t *data;

    size_t n_quads;
    size_t si;

    /* loop all quadrants and reset tag field to default */
    if (keep_pt_tag == 1) {
        for (t = first_local_tree; t <= last_local_tree; ++t) {
            tree = p4est_tree_array_index(trees,t);
            quadrants = &(tree->quadrants);
            n_quads = quadrants->elem_count;

            for (si = 0; si < n_quads; ++si) {
                quad = p4est_quadrant_array_index(quadrants,si);
                data = (quad_data_t *) quad->p.user_data;
                if(data->tag != point_tag) data->tag = tag_default;
            }
        }
    } else if (keep_pt_tag == -1) {
        for (t = first_local_tree; t <= last_local_tree; ++t) {
            tree = p4est_tree_array_index(trees,t);
            quadrants = &(tree->quadrants);
            n_quads = quadrants->elem_count;

            for (si = 0; si < n_quads; ++si) {
                quad = p4est_quadrant_array_index(quadrants,si);
                data = (quad_data_t *) quad->p.user_data;
                if(data->tag == point_tag) data->tag = tag_default;
            }
        }
    } else {
        /* all tags get reset */
        for (t = first_local_tree; t <= last_local_tree; ++t) {
            tree = p4est_tree_array_index(trees,t);
            quadrants = &(tree->quadrants);
            n_quads = quadrants->elem_count;

            for (si = 0; si < n_quads; ++si) {
                quad = p4est_quadrant_array_index(quadrants,si);
                data = (quad_data_t *) quad->p.user_data;
                data->tag = tag_default;
            }
        }
    }
    return 0;
}

int tag_reset_max_level(p4est_t *p4est,int tag_default){
    dg4est_t *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t       *ctx = &dg4est->d_ctx;
    grid_t     *grid = &ctx->d_grid;

    p4est_topidx_t t;
    p4est_topidx_t first_local_tree = p4est->first_local_tree;
    p4est_topidx_t last_local_tree = p4est->last_local_tree;
    sc_array_t *trees = p4est->trees;
    sc_array_t *quadrants;
    p4est_tree_t *tree;

    p4est_quadrant_t *quad;
    quad_data_t *data;

    size_t n_quads;
    size_t si;

    /* loop all quadrants and reset tag field to default */
    for (t = first_local_tree; t <= last_local_tree; ++t) {
        tree = p4est_tree_array_index(trees,t);
        quadrants = &(tree->quadrants);
        n_quads = quadrants->elem_count;

        for (si = 0; si < n_quads; ++si) {
            quad = p4est_quadrant_array_index(quadrants,si);
            if (quad->level > grid->max_level_pmax) {
                data = (quad_data_t *) quad->p.user_data;
                data->tag = tag_default;
            }
        }
    }
    return 0;
}

int tag_point(p4est_t *p4est,
              p4est_topidx_t which_tree,
              p4est_quadrant_t *q){

    /* suppress unused variable warning */
    (void) which_tree;

    quad_data_t *data = (quad_data_t *) q->p.user_data;
    dg4est_t  *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t        *ctx = &dg4est->d_ctx;

    if(data->tag == no_tag) return 0;
    if(q->level == ctx->d_grid.max_level) return 0;
    if(data->tag == point_tag) return 1;
    return 0;
}

static
char overlapping1D(bound_t box1,bound_t box2){
    return ((box1.hi >= box2.lo) && (box2.hi >= box1.lo));
}

int tag_point_search(p4est_t *p4est,
                     p4est_topidx_t which_tree,
                     p4est_quadrant_t *quadrant,
                     p4est_locidx_t local_num,
                     void *point){

    dg4est_t  *dg4est = (dg4est_t *) p4est->user_pointer;
    wake3d_t    *wake = &dg4est->d_ctx.d_wake3d;

    quad_data_t *data = (quad_data_t *) quadrant->p.user_data;
    Real          *pt = (Real *) point;

    Real xyz[3*P4EST_CHILDREN];
    int display_quad = 0;
    int quad_nnodes = P4EST_CHILDREN;
    int i;

    box_t box1,box2;
    Real vmin,vmax;
    Real x,y,z,h;
    char hit;

    /* ================================================================= */
    /* BOX 1: fill in axis-aligned bounding box of unstructured quadrant */
    /* ================================================================= */
    p4est_utilities_quad_coordinates(p4est,which_tree,quadrant,xyz,display_quad);

    /* calculate bounding volume: x */
    vmin = xyz[0];
    vmax = xyz[0];
    for (i = 1; i < quad_nnodes; ++i) {
        x = xyz[3*i+0];
        vmin = (vmin < x) ? vmin:x;
        vmax = (vmax > x) ? vmax:x;
    }
    box1.x.lo = vmin;
    box1.x.hi = vmax;

    /* calculate bounding volume: y */
    vmin = xyz[1];
    vmax = xyz[1];
    for (i = 1; i < quad_nnodes; ++i) {
        y = xyz[3*i+1];
        vmin = (vmin < y) ? vmin:y;
        vmax = (vmax > y) ? vmax:y;
    }
    box1.y.lo = vmin;
    box1.y.hi = vmax;

    /* calculate bounding volume: z */
#ifdef P4_TO_P8
    vmin = xyz[2];
    vmax = xyz[2];
    for (i = 1; i < quad_nnodes; ++i) {
        z = xyz[3*i+2];
        vmin = (vmin < z) ? vmin:z;
        vmax = (vmax > z) ? vmax:z;
    }
    box1.z.lo = vmin;
    box1.z.hi = vmax;
#endif

    /* ====================================================== */
    /* BOX 2: fill in axis-aligned bounding box of IGBP point */
    /* ====================================================== */
    /* point mesh size: cube root of pt which is a volume */
    h = wake->igbp_scale*cbrt(pt[P4EST_DIM]);

    /* fill in box1: box around pt */
    box2.x.lo = pt[0] - h;
    box2.x.hi = pt[0] + h;

    box2.y.lo = pt[1] - h;
    box2.y.hi = pt[1] + h;
#ifdef P4_TO_P8
    box2.z.lo = pt[2] - h;
    box2.z.hi = pt[2] + h;
#endif

    /* ================================ */
    /* Check Bounding Box Intersections */
    /* ================================ */
    hit = overlapping1D(box1.x,box2.x) && overlapping1D(box1.y,box2.y);
#ifdef P4_TO_P8
    hit = hit && overlapping1D(box1.z,box2.z);
#endif

    if (hit) {
        /* hit but not a leaf */
        if(local_num == -1) return 1;

        /* hit and is a leaf */
        data->tag = point_tag;
        return 1;
    }
    /* not hit */
    return 0;
}

int tag_level(p4est_t *p4est,
              p4est_topidx_t which_tree,
              p4est_quadrant_t *q){

    /* suppress unused variable warning */
    (void) which_tree;

    /* access our simulation context data */
    dg4est_t *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t       *ctx = &dg4est->d_ctx;
    grid_t     *grid = &ctx->d_grid;

    /* access our per quadrant data pointer */
    quad_data_t *data = (quad_data_t *) q->p.user_data;

    /* check if level is less than min_level */
    if (q->level < grid->min_level) {
        data->tag = level_tag;
        return 1;
    } else{
        return 0;
    }
}

int tag_feature(p4est_t *p4est,
                p4est_topidx_t which_tree,
                p4est_quadrant_t *q){
    Real dx[3] = {0.0};
    int max_level;
    int tag;
    int i;

    quad_properties_t quad_prop;
    quad_user_data_t user_data;

    /* access our simulation context data */
    dg4est_t   *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t         *ctx = &dg4est->d_ctx;
    grid_t       *grid = &ctx->d_grid;
    simulation_t  *sim = &ctx->d_simulation;
    amr_hbox_t *hboxes = &ctx->d_amr_hboxes;

    /* access our per quadrant data pointer */
    quad_data_t *data = (quad_data_t *) q->p.user_data;

    /* assign quad properties */
    quad_prop.dim          = sim->dim;
    quad_prop.unstructured = sim->unstructured_flag;
    quad_prop.level        = q->level;

    /* quadrant length in each direction */
    quad_prop.quad_dx = dx; // zero
    quad_prop.qdegree  = data->qdegree;
    quad_prop.quad_xyz = data->geom;

    /* ---------------------------------------------------------------------- */
    max_level = grid->max_level;

    /* ================================================== */
    /* No Feature but may be in AMR h-adaption box region */
    /* ================================================== */
    if (sim->regrid_hboxes) {
        int inside = 0;
        int force_tag = 0;
        for (i = 0; i < hboxes->amr_nhboxes; ++i) {
            force_tag |= hboxes->amr_hboxes_autotag[i];
            if(hboxes->amr_hboxes_autotag[i] == 0) continue; // skip box if not required

            /* forced refinement detected for this box */
            int box_level = hboxes->amr_hboxes_level[i];
            Real *lo = &hboxes->amr_hboxes_lo[3*i];
            Real *hi = &hboxes->amr_hboxes_hi[3*i];

            bbox2_t box_in;
            box_in.bbox_lo[0] = lo[0];     box_in.bbox_hi[0] = hi[0];
            box_in.bbox_lo[1] = lo[1];     box_in.bbox_hi[1] = hi[1];
        E3D(box_in.bbox_lo[2] = lo[2]) E3D(box_in.bbox_hi[2] = hi[2])

            int inside_loc = p4est_utilities_unst_bbox2_search_func(p4est,which_tree,q,1,&box_in);

            /* skip box if refined above box_level */
            if (inside_loc && q->level >= box_level) {
                data->tag = feature_tag; // override tag to prevent coarsening
                continue;
            }

            if(inside_loc) max_level = MAX(max_level,box_level);
            inside |= inside_loc;
        }

        if (inside && force_tag) {
            /* forced refinement detected: refine if not at max level */
            if (q->level < max_level) {
                data->tag = feature_tag; // override this tag for force tag
                return 1; // refine this quadrant
            } else {
                return 0; // inside but already at max level
            }
        } else {
            return 0; // no feature and no forced refinement
        }
    } else {
        return 0; // no boxes and no feature detected
    }
}

int tag_spread(p4est_t *p4est,
               p4est_topidx_t which_tree,
               p4est_quadrant_t *q){
    quad_data_t *data = (quad_data_t *) q->p.user_data;
    dg4est_t  *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t        *ctx = &dg4est->d_ctx;

    /* suppress compiler warning */
    (void) p4est;
    (void) which_tree;

    if (data->tag == spread_tag) {
        /* on initialization: delete the quadrant data */
        if (ctx->initial) {
            if(data->soln) {free(data->soln);} data->soln = NULL;
            data->soln_size = 0;
        }
        return 1;
    }
    return 0;
}

int tag_buffer(p4est_t *p4est,
               p4est_topidx_t which_tree,
               p4est_quadrant_t *q){
    quad_data_t *data = (quad_data_t *) q->p.user_data;
    dg4est_t  *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t        *ctx = &dg4est->d_ctx;

    /* suppress compiler warning */
    (void) p4est;
    (void) which_tree;

    if (data->tag == buffer_tag) {
        /* on initialization: delete the quadrant data */
        if (ctx->initial) {
            if(data->soln) {free(data->soln);} data->soln = NULL;
            data->soln_size = 0;
        }
        data->tag = no_tag;
        return 1;
    }
    return 0;
}

int tag_coarsen(p4est_t *p4est,
                p4est_topidx_t which_tree,
                p4est_quadrant_t *children[]){
    dg4est_t   *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t         *ctx = &dg4est->d_ctx;
    simulation_t  *sim = &ctx->d_simulation;
    amr_hbox_t *hboxes = &ctx->d_amr_hboxes;

    quad_data_t *data;
    int i;

    /* loop over the children and if any children igbp tagged do not coarsen */
    for (i = 0; i < P4EST_CHILDREN; ++i) {
        if (children[i]) {
            if(children[i]->level == ctx->d_grid.min_level) return 0; /* lowest: do not coarsen */

            if (children[i]->level > ctx->d_grid.max_level) {
                /* above max amr level, verify inside an h-box */
                if (sim->regrid_hboxes) {
                    int inside = 0;
                    for (i = 0; i < hboxes->amr_nhboxes; ++i) {
                        int box_level = hboxes->amr_hboxes_level[i];
                        if(box_level != children[i]->level) continue; // next box

                        /* matching hbox level: check if inside hbox */
                        Real *lo = &hboxes->amr_hboxes_lo[3*i];
                        Real *hi = &hboxes->amr_hboxes_hi[3*i];

                        bbox2_t box_in;
                        box_in.bbox_lo[0] = lo[0];     box_in.bbox_hi[0] = hi[0];
                        box_in.bbox_lo[1] = lo[1];     box_in.bbox_hi[1] = hi[1];
                    E3D(box_in.bbox_lo[2] = lo[2]) E3D(box_in.bbox_hi[2] = hi[2])

                        int inside_loc = p4est_utilities_unst_bbox2_search_func(p4est,which_tree,children[i],-1,&box_in);
                        inside |= inside_loc;
                    }
                    if(!inside) return 1; /* coarsen since not inside any hboxes */
                    return 0; /* all inside hbox: do not coarsen */
                }
                return 1; /* coarsen */
            } else {
                /* min/max amr level range: check if tagged */
                data = (quad_data_t *) children[i]->p.user_data;
                if(data->tag > no_tag) return 0;
            }
        }
    }
    /* all children are tagged 0 */
    return 1;
}
