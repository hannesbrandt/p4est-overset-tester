/**
 * \file    amr_p4est_utilities.cxx
 * \author  akirby
 *
 * \brief   Utility functions for accessing info about p4est data.
 */

/* header files */
#include "amr_p4est_utilities.h"

void p4est_utilities_mesh_stats(p4est_t *p4est){
    dg4est_t *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t        *ctx = &dg4est->d_ctx;

    int quad_total = 0;
    MPI_Reduce(&p4est->local_num_quadrants,
               &quad_total,
               1,MPI_INT,MPI_SUM,0,
               p4est->mpicomm);
    ctx->d_grid.total_quads = quad_total;

    if (ctx->rank == 0 && ctx->log_info < DGLOG_MESH_OFF) {
        DGOUT(ctx->log_io,"\n");
        DGOUT(ctx->log_io,"+====================================================+\n");
        DGOUT(ctx->log_io," [dg4est] Mesh Statistics:\n");
        DGOUT(ctx->log_io,"          ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n");
        DGOUT(ctx->log_io," Mesh Totals: \n");
        DGOUT(ctx->log_io," ‾‾‾‾‾‾‾‾‾‾‾               \n");
        DGOUT(ctx->log_io,"  Cell Total: %d\n",quad_total);
        DGOUT(ctx->log_io,"\n");
    }
}

void p4est_utilities_quad_coordinates(p4est_t *p4est,p4est_topidx_t which_tree,
                                      p4est_quadrant_t *q,Real *xyz,
                                      int display_quad){
    p4est_topidx_t vt[P4EST_CHILDREN];
    p4est_connectivity_t *connectivity;
    const Real *v;
    const p4est_topidx_t *tree_to_vertex;
    const Real intsize = 1.0 / (Real) P4EST_ROOT_LEN;

    int j,k;
    int xi,yi;
    Real h2,eta_x,eta_y,eta_z;

    connectivity = p4est->connectivity;
    v = connectivity->vertices;
    tree_to_vertex = connectivity->tree_to_vertex;

    /* fill in element size */
    h2 = 0.5 * intsize * P4EST_QUADRANT_LEN(q->level);

    for (k = 0; k < P4EST_CHILDREN; ++k) {
        vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
    }

    eta_z = 0.0;
    /* fill in quadrant node coordinates */
    k = 0;
#ifdef P4_TO_P8
    int zi;
    for (zi = 0; zi < 2; ++zi) {
        eta_z = intsize * q->z + h2 * (1.0 + (2*zi - 1));
#endif
        for (yi = 0; yi < 2; ++yi) {
            eta_y = intsize * q->y + h2 * (1.0 + (2*yi - 1));
            for (xi = 0; xi < 2; ++xi) {
                eta_x = intsize * q->x + h2 * (1.0 + (2*xi - 1));

                for (j = 0; j < 3; ++j) {
                    xyz[3*k+j] =
                        ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                               eta_x  * v[3 * vt[1] + j]) +
                                               eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                               eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
                         +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                               eta_x  * v[3 * vt[5] + j]) +
                                               eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                               eta_x  * v[3 * vt[7] + j]))
#endif
                        );
                }
                if (display_quad) {
                    DGOUT(stdout,"Quad geometry: %f %f %f\n",xyz[3*k+0],xyz[3*k+1],xyz[3*k+2]);
                }
                k++;
            }
        }
#ifdef P4_TO_P8
    }
#endif
    if(display_quad) DGOUT(stdout,"----------------\n");
}

void p4est_utilities_init_quad_coordinates(p4est_t *p4est,Real *xyz,
                                           ivector_t *quad_nodes,
                                           int display_quad){
    p4est_connectivity_t *connectivity;
    const Real *v;

    Uint inode,icoord;
    Uint inode_id;

    connectivity = p4est->connectivity;
    v = connectivity->vertices;

    /* Loop through the nodes of the quadrant */
    /* NOTE: the nodes are already z-ordered */
    for (inode = 0; inode < quad_nodes->size; inode++) {
        /* node ID */
        inode_id = ivector_get_value(quad_nodes,inode);

        /* node coordinates from conn->vertices */
        for(icoord = 0; icoord < 3; icoord++) xyz[3*inode+icoord] = v[3*inode_id+icoord];

        /* display node coordinates*/
        if (display_quad) {
            DGOUT(stdout,"Quad node local ID: %llu of %llu \n", inode,quad_nodes->size-1);
            DGOUT(stdout,"Quad node global ID: %llu \n", inode_id);
            DGOUT(stdout,"Quad geometry: %.15e %.15e %.15e\n",xyz[3*inode+0],xyz[3*inode+1],xyz[3*inode+2]);
        }
    }
    if(display_quad) DGOUT(stdout,"--------------------------\n");
}

void p4est_utilities_indep_node_coordinates(Real *xgeom,p4est_t *p4est,
                                            sc_array_t *indeps){
    p4est_connectivity_t *connectivity = p4est->connectivity;
    const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
    const Real *v = connectivity->vertices;

    p4est_indep_t *in;
    p4est_topidx_t vt[P4EST_CHILDREN];
    p4est_topidx_t jt;
    const Real intsize = 1.0 / P4EST_ROOT_LEN;
    Real eta_x, eta_y, eta_z=0.0;
    Real xyz[3];
    int i,j,k;

    /* unique node coordinates */
    for (i = 0; i < indeps->elem_count; ++i) {
        in = (p4est_indep_t *) sc_array_index(indeps,i);
        jt = in->p.which_tree;
        for (k = 0; k < P4EST_CHILDREN; ++k) {
            vt[k] = tree_to_vertex[P4EST_CHILDREN*jt + k];
        }

        /* calculate vertex coordinates */
        eta_x = intsize * in->x;
        eta_y = intsize * in->y;
    E3D(eta_z = intsize * in->z)

        for (j = 0; j < 3; ++j) {
            xyz[j] =
                ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                       eta_x  * v[3 * vt[1] + j]) +
                                       eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                       eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
               +       eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                       eta_x  * v[3 * vt[5] + j]) +
                                       eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                       eta_x  * v[3 * vt[7] + j]))
#endif
            );
        }

        xgeom[i*P4EST_DIM+0] = xyz[0];
        xgeom[i*P4EST_DIM+1] = xyz[1];
    E3D(xgeom[i*P4EST_DIM+2] = xyz[2])
    }
}

void p4est_utilities_get_xlo(p4est_t *p4est,p4est_topidx_t which_tree,
                             p4est_quadrant_t *q,Real xyz[3]){
    p4est_qcoord_to_vertex(p4est->connectivity,which_tree,
                           q->x,q->y,ARG3D(q->z)
                           xyz);
}

p4est_locidx_t p4est_utilities_get_local_id_volume(p4est_iter_volume_info_t *info){
    p4est_t            *p4est;
    p4est_topidx_t     which_tree;
    p4est_locidx_t     local_id;
    p4est_tree_t       *tree;

    p4est = info->p4est;
    which_tree = info->treeid;
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    local_id = info->quadid + tree->quadrants_offset;

    return local_id;
}

p4est_locidx_t p4est_utilities_get_local_id_face_hang(p4est_iter_face_info_t *info,
                                                      int which_side,
                                                      int which_hang){
    p4est_t                 *p4est;
    p4est_iter_face_side_t  *side;
    sc_array_t              *sides;
    p4est_topidx_t          which_tree;
    p4est_locidx_t          local_id;
    p4est_tree_t            *tree;

    p4est = info->p4est;
    sides = &(info->sides);
    side = p4est_iter_fside_array_index_int (sides, which_side);
    which_tree = side->treeid;
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    local_id = side->is.hanging.quadid[which_hang] + tree->quadrants_offset;

    return local_id;
}

p4est_locidx_t p4est_utilities_get_local_id_face_full(p4est_iter_face_info_t *info,
                                                      int which_side){
    p4est_t                 *p4est;
    p4est_iter_face_side_t  *side;
    sc_array_t              *sides;
    p4est_topidx_t          which_tree;
    p4est_locidx_t          local_id;
    p4est_tree_t            *tree;

    p4est = info->p4est;
    sides = &(info->sides);
    side = p4est_iter_fside_array_index_int (sides, which_side);
    which_tree = side->treeid;
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    local_id = side->is.full.quadid + tree->quadrants_offset;

    return local_id;
}

void p4est_utilities_get_quadrant_from_local_id(p4est_t *p4est,
                                                p4est_locidx_t local_id,
                                                p4est_quadrant_t **quad,
                                                p4est_topidx_t *which_tree){
    p4est_topidx_t first_local_tree = p4est->first_local_tree;
    p4est_topidx_t last_local_tree = p4est->last_local_tree;
    p4est_topidx_t t;
    p4est_tree_t *tree;
    sc_array_t *quadrants;
    sc_array_t *trees = p4est->trees;
    size_t n_quads;

    for (t = first_local_tree; t <= last_local_tree; ++t) {
        tree = p4est_tree_array_index(trees, t);
        quadrants = &(tree->quadrants);
        n_quads = quadrants->elem_count;
        if (local_id < tree->quadrants_offset + n_quads) {
            *quad = p4est_quadrant_array_index(quadrants,local_id-tree->quadrants_offset);
            *which_tree = t;
            return;
        }
    }
}

static
void pack_quad_iblank_ghosts(p4est_t *p4est,p4est_ghost_t *ghost,
                             char *mirror_data,void **mirror_data_ptr){
    p4est_topidx_t which_tree;
    p4est_locidx_t which_quad;
    p4est_quadrant_t *mirror;
    p4est_quadrant_t *q;
    p4est_tree_t *tree;
    quad_data_t *data;
    size_t zz;

    /* fill in the mirror data with the AMR quad type specified by the user */
    /* loop all elements that touch processor boundary */
    for (zz = 0; zz < ghost->mirrors.elem_count; ++zz) {
        mirror = p4est_quadrant_array_index(&ghost->mirrors,zz);
        which_tree = mirror->p.piggy3.which_tree;

        tree = p4est_tree_array_index(p4est->trees,which_tree);
        which_quad = mirror->p.piggy3.local_num - tree->quadrants_offset;

        q = p4est_quadrant_array_index (&tree->quadrants, which_quad);
        data = (quad_data_t *) q->p.user_data;

        mirror_data[zz] = (char) data->iblank;
        mirror_data_ptr[zz] = &mirror_data[zz];
    }
}


void iblank_ext_to_p4est_callback(p4est_iter_volume_info_t *info,void *user_data){
    ctx_t             *ctx = (ctx_t *) info->p4est->user_pointer;
    wake3d_t         *wake = &ctx->d_wake3d;
    p4est_quadrant_t *quad = info->quad;
    quad_data_t      *data = (quad_data_t *) quad->p.user_data;

    p4est_locidx_t local_id = p4est_utilities_get_local_id_volume(info);
    data->iblank = wake->iblank_cell[local_id];
}

static
void iblank_p4est_to_ext_callback(p4est_iter_volume_info_t *info,void *user_data){
    ctx_t             *ctx = (ctx_t *) info->p4est->user_pointer;
    wake3d_t         *wake = &ctx->d_wake3d;
    p4est_quadrant_t *quad = info->quad;
    quad_data_t      *data = (quad_data_t *) quad->p.user_data;

    /* suppress compiler warning */
    (void) user_data;

    p4est_locidx_t local_id = p4est_utilities_get_local_id_volume(info);
    wake->iblank_cell[local_id] = data->iblank;
}

static
void iblank_callback(p4est_iter_face_info_t *info,void *user_data){
    p4est_t   *p4est = (p4est_t *) info->p4est;
    dg4est_t *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t       *ctx = &dg4est->d_ctx;

    iblank_data_t *iblank_data = (iblank_data_t *) user_data;
    char *ghost_data_iblanks = iblank_data->ghost_data;
    char correct_flag        = iblank_data->correct_flag;
    char display_flag        = iblank_data->display_flag;

    char ldata_iblank,rdata_iblank,hdata_iblank,fdata_iblank;
    quad_data_t *ldata,*rdata,*hdata,*fdata;
    sc_array_t *sides = &(info->sides);

    p4est_iter_face_side_t *side[2];
    p4est_quadrant_t *quad;
    p4est_topidx_t treeid;
    p4est_locidx_t local_id;

    int hangside,fullside,hanging;
    int j;

    /* non-periodic boundary face */
    if(sides->elem_count == 1) return;
    P4EST_ASSERT(sides->elem_count == 2);

    side[0] = p4est_iter_fside_array_index_int(sides, 0);
    side[1] = p4est_iter_fside_array_index_int(sides, 1);

    if (side[0]->is_hanging) {
        hangside = 0;
        fullside = 1;
        hanging  = 1;
    } else
    if (side[1]->is_hanging) {
        hangside = 1;
        fullside = 0;
        hanging  = 1;
    } else {
        hanging  = 0;
        hangside = 0; /* suppress compiler warning */
        fullside = 0; /* suppress compiler warning */
    }

    if (hanging) {
        if (side[fullside]->is.full.is_ghost) {
            fdata_iblank = ghost_data_iblanks[side[fullside]->is.full.quadid];
        } else {
            fdata = (quad_data_t *) side[fullside]->is.full.quad->p.user_data;
            fdata_iblank = fdata->iblank;
        }

        for (j = 0; j < P4EST_HALF; j++) {
            if (side[hangside]->is.hanging.is_ghost[j]) {
                hdata_iblank = ghost_data_iblanks[side[hangside]->is.hanging.quadid[j]];
            } else {
                hdata = (quad_data_t *) side[hangside]->is.hanging.quad[j]->p.user_data;
                hdata_iblank = hdata->iblank;
            }

            if (hdata && fdata) {
                if ( (hdata_iblank == hole_mask && fdata_iblank == field_mask) ||
                     (fdata_iblank == hole_mask && hdata_iblank == field_mask)) {
#ifdef WARNINGS
                    DGOUT(stderr,"[ amr ] WARNING: Hole next to field on HANG Face; Rank: %d\n",
                           ctx->rank);
#endif
                    if (display_flag) {
                        if (!side[fullside]->is.full.is_ghost) {
                            quad = side[fullside]->is.full.quad;
                            treeid = side[fullside]->treeid;
                            local_id = p4est_utilities_get_local_id_face_full(info,fullside);
#ifdef WARNINGS
                            DGOUT(stderr,"            \t Rank[%d] Full-Side Quad ID: %d\n",ctx->rank,local_id);
#endif
                        }
                        if (!side[hangside]->is.hanging.is_ghost[j]) {
                            quad = side[hangside]->is.hanging.quad[j];
                            treeid = side[hangside]->treeid;
                            local_id = p4est_utilities_get_local_id_face_hang(info,hangside,j);
#ifdef WARNINGS
                            DGOUT(stderr,"            \t Rank[%d] Hang-Side[%d] Quad ID: %d\n",ctx->rank,j,local_id);
#endif
                        }
                    }

                    if (correct_flag) {
                        if (!side[fullside]->is.full.is_ghost &&
                            fdata->iblank == field_mask){
                            fdata->iblank = receptor_mask;
                        }
                        if (!side[hangside]->is.hanging.is_ghost[j] &&
                            hdata->iblank == field_mask){
                            hdata->iblank = receptor_mask;
                        }
                    }
                }
            }
        }

    } else {

        if (side[0]->is.full.is_ghost) {
            ldata_iblank = ghost_data_iblanks[side[0]->is.full.quadid];
        } else {
            ldata = (quad_data_t *) side[0]->is.full.quad->p.user_data;
            ldata_iblank = ldata->iblank;
        }

        if (side[1]->is.full.is_ghost) {
            rdata_iblank = ghost_data_iblanks[side[1]->is.full.quadid];
        } else {
            rdata = (quad_data_t *) side[1]->is.full.quad->p.user_data;
            rdata_iblank = rdata->iblank;
        }

        if (ldata && rdata) {
            if ( (ldata_iblank == hole_mask && rdata_iblank == field_mask) ||
                 (rdata_iblank == hole_mask && ldata_iblank == field_mask)) {
#ifdef WARNINGS
                DGOUT(stderr,"[ amr ] WARNING: Hole next to field on FULL Face; Rank: %d\n",
                       ctx->rank);
#endif
                if (display_flag) {
                    if (!side[0]->is.full.is_ghost) {
                        quad = side[0]->is.full.quad;
                        treeid = side[0]->treeid;
                        local_id = p4est_utilities_get_local_id_face_full(info,0);
#ifdef WARNINGS
                        DGOUT(stderr,"            \t Rank[%d] Full-Side Quad ID: %d\n",ctx->rank,local_id);
#endif
                    }
                    if (!side[1]->is.full.is_ghost) {
                        quad = side[1]->is.full.quad;
                        treeid = side[1]->treeid;
                        local_id = p4est_utilities_get_local_id_face_full(info,1);
#ifdef WARNINGS
                        DGOUT(stderr,"            \t Rank[%d] Full-Side Quad ID: %d\n",ctx->rank,local_id);
#endif
                    }
                }

                if (correct_flag) {
                    if (!side[0]->is.full.is_ghost &&
                        ldata->iblank == field_mask) {
                        ldata->iblank = receptor_mask;
                    }
                    if (!side[1]->is.full.is_ghost &&
                        rdata->iblank == field_mask) {
                        rdata->iblank = receptor_mask;
                    }
                }
            }
        }
    }
}

char intersect_box(vec3 origin,vec3 dir,const box b){
    Real txMin,txMax;
    Real tyMin,tyMax;
    Real tzMin,tzMax;
    Real tMin,tMax;

    Real eps = 1.0e-14;
    Real tmp;
    int i;

    char inside = 1;
    char pt = 1;

    /* check if ray origin is inside box */
    for (i = 0; i < P4EST_DIM; ++i) {
        if(origin.x[i] < b.min.x[i] || origin.x[i] > b.max.x[i]) inside = 0;
    }
    if(inside) return 1;

    /* compute the end point of ray and check if the ray is just a point */
    for (i = 0; i < P4EST_DIM; ++i) {
        dir.x[i] = dir.x[i] - origin.x[i]; // make relative to origin
        if(dir.x[i] > eps) pt = 0;
    }
    /* if ray is just a point:
     *    note: the origin is not in the box -> ray is not in the box
     */
    if(pt) return 0;

    /* compute parameterized distances to the box in each direction */
    /* x-direction */
    if (dir.x[0] != 0.0) {
        const Real oneOdir = 1.0 / dir.x[0];
        txMin = (b.min.x[0] - origin.x[0]) * oneOdir;
        txMax = (b.max.x[0] - origin.x[0]) * oneOdir;
    } else {
        txMin = -1.0;
        txMax = -1.0;
    }
    if(txMax < txMin) {tmp = txMax; txMax = txMin; txMin = tmp;}

    /* y-direction */
    if (dir.x[1] != 0.0) {
        const Real oneOdir = 1.0 / dir.x[1];
        tyMin = (b.min.x[1] - origin.x[1]) * oneOdir;
        tyMax = (b.max.x[1] - origin.x[1]) * oneOdir;
    } else {
        tyMin = -1.0;
        tyMax = -1.0;
    }
    if(tyMax < tyMin) {tmp = tyMax; tyMax = tyMin; tyMin = tmp;}

    /* z-direction */
    tzMin = -DBL_MAX;
    tzMax =  DBL_MAX;
    if (P4EST_DIM == 3) {
        if (dir.x[2] != 0.0) {
            const Real oneOdir = 1.0 / dir.x[2];
            tzMin = (b.min.x[2] - origin.x[2]) * oneOdir;
            tzMax = (b.max.x[2] - origin.x[2]) * oneOdir;
        } else {
            tzMin = -1.0;
            tzMax = -1.0;
        }
        if(tzMax < tzMin) {tmp = tzMax; tzMax = tzMin; tzMin = tmp;}
    }

    /* find the largest minimum and the smallest maximum t values */
    tMin = (txMin > tyMin) ? txMin : tyMin;
    tMax = (txMax < tyMax) ? txMax : tyMax;

    if(txMin > tyMax || tyMin > txMax) return 0;
    if(tMin  > tzMax || tzMin > tMax)  return 0;
    if(tzMin > tMin) tMin = tzMin;
    if(tzMax < tMax) tMax = tzMax;

    /* if point has min t value parameterized between 0 and 1, then hit */
    return (tMin > 0.0 && tMin <= 1.0 && (tMin < tMax));
}

static
char overlapping1D(bound_t box1,bound_t box2){
    return ((box1.hi >= box2.lo) && (box2.hi >= box1.lo));
}

static inline
void linear_phi_1d(Real r,Real phi[2]){
    /* 1D linear basis function */
    phi[0] = 0.5*(1.0 - r);
    phi[1] = 0.5*(1.0 + r);
}

static inline
void linear_dphi_1d(Real phi[2]){
    /* 1D linear basis derivative function */
    phi[0] = -0.5;
    phi[1] =  0.5;
}

static inline
void linear_phi_2d(Real r,Real s,Real phi[2][2]){
    Real phi1[2];
    Real phi2[2];
    int i,j;

    /* evaluate basis functions */
    linear_phi_1d(r,phi1);
    linear_phi_1d(s,phi2);

    /* form tensor-product basis */
    for (j = 0; j < 2; ++j) {
        for (i = 0; i < 2; ++i) {
            phi[j][i] = phi1[i]*phi2[j];
        }
    }
}

static inline
void linear_dphi_2d(Real r,Real s,char deriv_index,Real phi[2][2]){
    Real phi1[2];
    Real phi2[2];
    int i,j;

    /* evaluate basis functions */
    switch (deriv_index) {
        case(1):
            linear_dphi_1d(phi1);
            linear_phi_1d(s,phi2);
            break;
        case(2):
            linear_phi_1d(r,phi1);
            linear_dphi_1d(phi2);
            break;
    }

    /* form tensor-product basis */
    for (j = 0; j < 2; ++j) {
        for (i = 0; i < 2; ++i) {
            phi[j][i] = phi1[i]*phi2[j];
        }
    }
}

static inline
void linear_phi_3d(Real r,Real s,Real t,Real phi[2][2][2]){
    Real phi1[2];
    Real phi2[2];
    Real phi3[2];
    int i,j,k;

    /* evaluate basis functions */
    linear_phi_1d(r,phi1);
    linear_phi_1d(s,phi2);
    linear_phi_1d(t,phi3);

    /* form tensor-product basis */
    for (k = 0; k < 2; ++k) {
        for (j = 0; j < 2; ++j) {
            for (i = 0; i < 2; ++i) {
                phi[k][j][i] = phi1[i]*phi2[j]*phi3[k];
            }
        }
    }
}

static inline
void linear_dphi_3d(Real r,Real s,Real t,char deriv_index,Real phi[2][2][2]){
    Real phi1[2];
    Real phi2[2];
    Real phi3[2];
    int i,j,k;

    /* evaluate basis functions */
    switch (deriv_index) {
        case(1):
            linear_dphi_1d(phi1);
            linear_phi_1d(s,phi2);
            linear_phi_1d(t,phi3);
            break;
        case(2):
            linear_phi_1d(r,phi1);
            linear_dphi_1d(phi2);
            linear_phi_1d(t,phi3);
            break;
        case(3):
            linear_phi_1d(r,phi1);
            linear_phi_1d(s,phi2);
            linear_dphi_1d(phi3);
    }

    /* form tensor-product basis */
    for (k = 0; k < 2; ++k) {
        for (j = 0; j < 2; ++j) {
            for (i = 0; i < 2; ++i) {
                phi[k][j][i] = phi1[i]*phi2[j]*phi3[k];
            }
        }
    }
}

Real p4est_utilities_face_area_calc(int face,Real *quad_geom){
#ifdef P4_TO_P8
    Real xyz_quad[4][3];
    Real dphidr[2][2];
    Real dphids[2][2];

    Real detJ;
    Real dxdr,dxds;
    Real dydr,dyds;
    Real r,s;
    int xco,yco;
    int i,j;

    switch (face) {
        case(0): /* xlo face */
            xco = 2;
            yco = 1;
            for(i = 0; i < 3; ++i) xyz_quad[0][i] = quad_geom[3*0+i]; // node 0
            for(i = 0; i < 3; ++i) xyz_quad[1][i] = quad_geom[3*4+i]; // node 4
            for(i = 0; i < 3; ++i) xyz_quad[2][i] = quad_geom[3*2+i]; // node 2
            for(i = 0; i < 3; ++i) xyz_quad[3][i] = quad_geom[3*6+i]; // node 6
            break;
        case(1): /* xhi face */
            xco = 2;
            yco = 1;
            for(i = 0; i < 3; ++i) xyz_quad[0][i] = quad_geom[3*1+i]; // node 1
            for(i = 0; i < 3; ++i) xyz_quad[1][i] = quad_geom[3*5+i]; // node 5
            for(i = 0; i < 3; ++i) xyz_quad[2][i] = quad_geom[3*3+i]; // node 3
            for(i = 0; i < 3; ++i) xyz_quad[3][i] = quad_geom[3*7+i]; // node 7
            break;
        case(2): /* ylo face */
            xco = 0;
            yco = 2;
            for(i = 0; i < 3; ++i) xyz_quad[0][i] = quad_geom[3*0+i]; // node 0
            for(i = 0; i < 3; ++i) xyz_quad[1][i] = quad_geom[3*1+i]; // node 1
            for(i = 0; i < 3; ++i) xyz_quad[2][i] = quad_geom[3*4+i]; // node 4
            for(i = 0; i < 3; ++i) xyz_quad[3][i] = quad_geom[3*5+i]; // node 5
            break;
        case(3): /* yhi face */
            xco = 0;
            yco = 2;
            for(i = 0; i < 3; ++i) xyz_quad[0][i] = quad_geom[3*2+i]; // node 2
            for(i = 0; i < 3; ++i) xyz_quad[1][i] = quad_geom[3*3+i]; // node 3
            for(i = 0; i < 3; ++i) xyz_quad[2][i] = quad_geom[3*6+i]; // node 6
            for(i = 0; i < 3; ++i) xyz_quad[3][i] = quad_geom[3*7+i]; // node 7
            break;
        case(4): /* zlo face */
            xco = 0;
            yco = 1;
            for(i = 0; i < 3; ++i) xyz_quad[0][i] = quad_geom[3*0+i]; // node 0
            for(i = 0; i < 3; ++i) xyz_quad[1][i] = quad_geom[3*1+i]; // node 1
            for(i = 0; i < 3; ++i) xyz_quad[2][i] = quad_geom[3*2+i]; // node 2
            for(i = 0; i < 3; ++i) xyz_quad[3][i] = quad_geom[3*3+i]; // node 3
            break;
        case(5): /* zhi face */
            xco = 0;
            yco = 1;
            for(i = 0; i < 3; ++i) xyz_quad[0][i] = quad_geom[3*4+i]; // node 4
            for(i = 0; i < 3; ++i) xyz_quad[1][i] = quad_geom[3*5+i]; // node 5
            for(i = 0; i < 3; ++i) xyz_quad[2][i] = quad_geom[3*6+i]; // node 6
            for(i = 0; i < 3; ++i) xyz_quad[3][i] = quad_geom[3*7+i]; // node 7
            break;
    }
    /* center point for transformation */
    r = 0.0;
    s = 0.0;

    /* form Jacobian to get volume */
    linear_dphi_2d(r,s,1,dphidr);
    linear_dphi_2d(r,s,2,dphids);

    dxdr = 0.0; dxds = 0.0;
    dydr = 0.0; dyds = 0.0;

    int ind = 0;
    for (j = 0; j < 2; ++j) {
        for (i = 0; i < 2; ++i) {
            dxdr += dphidr[j][i]*xyz_quad[ind][xco];
            dxds += dphids[j][i]*xyz_quad[ind][xco];
            dydr += dphidr[j][i]*xyz_quad[ind][yco];
            dyds += dphids[j][i]*xyz_quad[ind][yco];
            ind++;
        }
    }

    detJ = (dxdr*dyds - dxds*dydr);
    return 4.0*detJ;
#else
    Real xyz_quad[2][2];
    Real dx,dy;
    int i;

    switch (face) {
        case(0): /* xlo face */
            for(i = 0; i < 2; ++i) xyz_quad[0][i] = quad_geom[3*0+i]; // node 0
            for(i = 0; i < 2; ++i) xyz_quad[1][i] = quad_geom[3*2+i]; // node 2
            break;
        case(1): /* xhi face */
            for(i = 0; i < 2; ++i) xyz_quad[0][i] = quad_geom[3*1+i]; // node 1
            for(i = 0; i < 2; ++i) xyz_quad[1][i] = quad_geom[3*3+i]; // node 3
            break;
        case(2): /* ylo face */
            for(i = 0; i < 2; ++i) xyz_quad[0][i] = quad_geom[3*0+i]; // node 0
            for(i = 0; i < 2; ++i) xyz_quad[1][i] = quad_geom[3*1+i]; // node 1
            break;
        case(3): /* yhi face */
            for(i = 0; i < 2; ++i) xyz_quad[0][i] = quad_geom[3*2+i]; // node 2
            for(i = 0; i < 2; ++i) xyz_quad[1][i] = quad_geom[3*3+i]; // node 3
            break;
    }

    dx = xyz_quad[1][0] - xyz_quad[0][0];
    dy = xyz_quad[1][1] - xyz_quad[0][1];
    return sqrt((dx*dx + dy*dy));
#endif
}

Real p4est_utilities_surface_area_calc(Real *xyz_quad){
    Real sarea = 0.0;
    int nface  = 2*DIM;
    int f;

    for (f = 0; f < nface; ++f) {
        sarea += p4est_utilities_face_area_calc(f,xyz_quad);
    }
    return sarea;
}

Real p4est_utilities_volume_calc(Real *xyz_quad){
#ifdef P4_TO_P8
    Real dphidr[2][2][2];
    Real dphids[2][2][2];
    Real dphidt[2][2][2];

    Real detJ;
    Real quad_volume;
    Real dxdr,dxds,dxdt;
    Real dydr,dyds,dydt;
    Real dzdr,dzds,dzdt;
    Real r,s,t;
    int i,j,k;

    /* center point for transformation */
    r = 0.0;
    s = 0.0;
    t = 0.0;

    /* form Jacobian to get volume */
    linear_dphi_3d(r,s,t,1,dphidr);
    linear_dphi_3d(r,s,t,2,dphids);
    linear_dphi_3d(r,s,t,3,dphidt);

    dxdr = 0.0; dxds = 0.0; dxdt = 0.0;
    dydr = 0.0; dyds = 0.0; dydt = 0.0;
    dzdr = 0.0; dzds = 0.0; dzdt = 0.0;

    int ind = 0;
    for (k = 0; k < 2; ++k) {
        for (j = 0; j < 2; ++j) {
            for (i = 0; i < 2; ++i) {
                const int node = 3*ind;
                dxdr += dphidr[k][j][i]*xyz_quad[node+0];
                dxds += dphids[k][j][i]*xyz_quad[node+0];
                dxdt += dphidt[k][j][i]*xyz_quad[node+0];
                dydr += dphidr[k][j][i]*xyz_quad[node+1];
                dyds += dphids[k][j][i]*xyz_quad[node+1];
                dydt += dphidt[k][j][i]*xyz_quad[node+1];
                dzdr += dphidr[k][j][i]*xyz_quad[node+2];
                dzds += dphids[k][j][i]*xyz_quad[node+2];
                dzdt += dphidt[k][j][i]*xyz_quad[node+2];
                ind++;
            }
        }
    }

    detJ = dxdr*(dyds*dzdt-dzds*dydt) -
           dxds*(dydr*dzdt-dzdr*dydt) +
           dxdt*(dydr*dzds-dzdr*dyds);

    quad_volume = 8.0*detJ;
    return quad_volume;
#else
    Real dphidr[2][2];
    Real dphids[2][2];

    Real detJ;
    Real quad_volume;
    Real dxdr,dxds;
    Real dydr,dyds;
    Real r,s;
    int i,j;

    /* center point for transformation */
    r = 0.0;
    s = 0.0;

    /* form Jacobian to get volume */
    linear_dphi_2d(r,s,1,dphidr);
    linear_dphi_2d(r,s,2,dphids);

    dxdr = 0.0; dxds = 0.0;
    dydr = 0.0; dyds = 0.0;
    int ind = 0;
    for (j = 0; j < 2; ++j) {
        for (i = 0; i < 2; ++i) {
            const int node = 3*ind;
            dxdr += dphidr[j][i]*xyz_quad[node+0];
            dxds += dphids[j][i]*xyz_quad[node+0];
            dydr += dphidr[j][i]*xyz_quad[node+1];
            dyds += dphids[j][i]*xyz_quad[node+1];
            ind++;
        }
    }

    detJ = (dxdr*dyds - dxds*dydr);

    quad_volume = 4.0*detJ;
    return quad_volume;
#endif
}

Real p4est_utilities_quad_volume(p4est_t *p4est,p4est_topidx_t which_tree,
                                 p4est_quadrant_t *quad){
    Real xyz_quad[3*P4EST_CHILDREN];
    Real volume;

    /* get quadrant geometry coordinates */
    p4est_utilities_quad_coordinates(p4est,which_tree,quad,xyz_quad,0);

    /* get quadrant volume */
    volume = p4est_utilities_volume_calc(xyz_quad);
    return volume;
}

int p4est_pt_rst_coordinates_2d(int qdegree,Real *geom,Real *xyz,Real *rst){
    /* =================================================== */
    /* FIXME: We are assuming bilinear geometry quadrants. */
    /* When high order geometry is incorporated,           */
    /* we will replace the built-in linear finite elements */
    /* with the high order elements from the flow solver.  */
    /* =================================================== */
    (void) qdegree;

    Real phi[2][2];
    Real dphidr[2][2];
    Real dphids[2][2];

    Real J_inv[2][2];
    Real R[2];

    Real tol = 1.0E-15;
    Real newton_tol = 1.0E-13;
    int max_iter = 100;
    int iter;

    Real res_norm;
    Real oneOdetj;
    Real dxdr,dxds;
    Real dydr,dyds;
    Real r,s;
    int i,j,d;
    int ind;

    /* initial guess for natural coordinates */
    r = 0.0;
    s = 0.0;

    /* Newton's method to find natural coordinates */
    for (iter = 0; iter < max_iter; ++iter) {

        /* evaluate basis and derivatives at rst */
        linear_phi_2d(r,s,phi);
        linear_dphi_2d(r,s,1,dphidr);
        linear_dphi_2d(r,s,2,dphids);

        /* form residual -R = -(phi_i(r,s)*xg_i - x) */
        for(d = 0; d < 2; ++d) R[d] = xyz[d];

        for (j = 0; j < 2; ++j) {
            for (i = 0; i < 2; ++i) {
                const int node = 2*j + i;

                for (d = 0; d < 2; ++d) {
                    R[d] -= phi[j][i]*geom[3*node + d]; // geom contains x,y,z
                }
            }
        }

        /* norm of residual */
        res_norm = 0.0;
        for (d = 0; d < 2; ++d) res_norm += fabs(R[d]);
        res_norm /= (Real) 2;

        /* break from Newton loop if converged */
        if (res_norm < newton_tol) break;

        /* form Jacobian = \partial(R) / \partial(r,s) */
        /*               = \dphi_i(r,s)*xg_i           */
        dxdr = 0.0; dxds = 0.0;
        dydr = 0.0; dyds = 0.0;

        ind = 0;
        for (j = 0; j < 2; ++j) {
            for (i = 0; i < 2; ++i) {
                const int node = 3*(2*j + i);
                dxdr += dphidr[j][i]*geom[node+0];
                dxds += dphids[j][i]*geom[node+0];
                dydr += dphidr[j][i]*geom[node+1];
                dyds += dphids[j][i]*geom[node+1];
            }
        }

        /* form inverse Jacobian */
        /* 1st row */
        J_inv[0][0] =  dyds; /* drdx */
        J_inv[0][1] = -dxds; /* drdy */

        /* 2nd row */
        J_inv[1][0] = -dydr; /* dsdx */
        J_inv[1][1] =  dxdr; /* dsdy */

        oneOdetj = 1.0/(dxdr*dyds - dxds*dydr);

        /* rst += matmul(J_inv,R) */
        r += oneOdetj*(J_inv[0][0]*R[0] + J_inv[0][1]*R[1]);
        s += oneOdetj*(J_inv[1][0]*R[0] + J_inv[1][1]*R[1]);
    }

    /* check if values are NAN (Newton diverged) */
    if(isnan(r)) return 0;
    if(isnan(s)) return 0;

    /* save natural coordinates */
    rst[0] = r;
    rst[1] = s;

    /* check bounds */
    if(fabs(r) - 1.0 > tol) return 0;
    if(fabs(s) - 1.0 > tol) return 0;

    /* found point in quadrant */
    return 1;
}

int p4est_pt_rst_coordinates_3d(int qdegree,Real *geom,Real *xyz,Real *rst){
    /* =================================================== */
    /* FIXME: We are assuming bilinear geometry quadrants. */
    /* When high order geometry is incorporated,           */
    /* we will replace the built-in linear finite elements */
    /* with the high order elements from the flow solver.  */
    /* =================================================== */
    (void) qdegree;

    Real phi[2][2][2];
    Real dphidr[2][2][2];
    Real dphids[2][2][2];
    Real dphidt[2][2][2];

    Real J_inv[3][3];
    Real R[3];

    Real tol = 1.0E-3;
    Real newton_tol = 1.0E-13;
    int max_iter = 10;
    int iter;

    Real res_norm;
    Real oneOdetj;
    Real dxdr,dxds,dxdt;
    Real dydr,dyds,dydt;
    Real dzdr,dzds,dzdt;
    Real r,s,t;
    int i,j,k,d;
    int ind;

    /* initial guess for natural coordinates */
    r = 0.0;
    s = 0.0;
    t = 0.0;

    /* Newton's method to find natural coordinates */
    for (iter = 0; iter < max_iter; ++iter) {

        /* evaluate basis and derivatives at rst */
        linear_phi_3d(r,s,t,phi);
        linear_dphi_3d(r,s,t,1,dphidr);
        linear_dphi_3d(r,s,t,2,dphids);
        linear_dphi_3d(r,s,t,3,dphidt);

        /* form residual -R = -(phi_i(r,s,t)*xg_i - x) */
        for (d = 0; d < 3; ++d) R[d] = 0.0;
        for (d = 0; d < P4EST_DIM; ++d) R[d] = xyz[d];

        for (k = 0; k < 2; ++k) {
            for (j = 0; j < 2; ++j) {
                for (i = 0; i < 2; ++i) {
                    const int node = 3*(4*k + 2*j + i);
                    for (d = 0; d < P4EST_DIM; ++d) {
                        R[d] -= phi[k][j][i]*geom[node + d];
                    }
                }
            }
        }

        /* norm of residual */
        res_norm = 0.0;
        for (d = 0; d < P4EST_DIM; ++d) res_norm += fabs(R[d]);
        res_norm /= (Real) P4EST_DIM;

        /* break from Newton loop if converged */
        if (res_norm < newton_tol) break;

        /* form Jacobian = \partial(R) / \partial(r,s,t) */
        /*               = \dphi_i(r,s,t)*xg_i           */
        dxdr = 0.0; dxds = 0.0; dxdt = 0.0;
        dydr = 0.0; dyds = 0.0; dydt = 0.0;
        dzdr = 0.0; dzds = 0.0; dzdt = 0.0;
        ind = 0;
        for (k = 0; k < 2; ++k) {
            for (j = 0; j < 2; ++j) {
                for (i = 0; i < 2; ++i) {
                    const int node = 3*(4*k + 2*j + i);
                    dxdr += dphidr[k][j][i]*geom[node+0];
                    dxds += dphids[k][j][i]*geom[node+0];
                    dxdt += dphidt[k][j][i]*geom[node+0];
                    dydr += dphidr[k][j][i]*geom[node+1];
                    dyds += dphids[k][j][i]*geom[node+1];
                    dydt += dphidt[k][j][i]*geom[node+1];
                    dzdr += dphidr[k][j][i]*geom[node+2];
                    dzds += dphids[k][j][i]*geom[node+2];
                    dzdt += dphidt[k][j][i]*geom[node+2];
                }
            }
        }

        /* form inverse Jacobian */
        /* 1st row */
        J_inv[0][0] =  dyds*dzdt - dydt*dzds; /* drdx */
        J_inv[0][1] = -dxds*dzdt + dxdt*dzds; /* drdy */
        J_inv[0][2] =  dxds*dydt - dxdt*dyds; /* drdz */

        /* 2nd row */
        J_inv[1][0] = -dydr*dzdt + dydt*dzdr; /* dsdx */
        J_inv[1][1] =  dxdr*dzdt - dxdt*dzdr; /* dsdy */
        J_inv[1][2] = -dxdr*dydt + dxdt*dydr; /* dsdz */

        /* 3rd row */
        J_inv[2][0] =  dydr*dzds - dyds*dzdr; /* dtdx */
        J_inv[2][1] = -dxdr*dzds + dxds*dzdr; /* dtdy */
        J_inv[2][2] =  dxdr*dyds - dxds*dydr; /* dtdz */

        oneOdetj = 1.0 /
                  (dxdr*(dyds*dzdt-dzds*dydt)
                  -dxds*(dydr*dzdt-dzdr*dydt)
                  +dxdt*(dydr*dzds-dzdr*dyds));

        /* rst += matmul(J_inv,R) */
        r += oneOdetj*(J_inv[0][0]*R[0] + J_inv[0][1]*R[1] + J_inv[0][2]*R[2]);
        s += oneOdetj*(J_inv[1][0]*R[0] + J_inv[1][1]*R[1] + J_inv[1][2]*R[2]);
        t += oneOdetj*(J_inv[2][0]*R[0] + J_inv[2][1]*R[1] + J_inv[2][2]*R[2]);
    }

    /* check if values are NAN (Newton diverged) */
    if(isnan(r)) return 0;
    if(isnan(s)) return 0;
E3D(if(isnan(t)) return 0;)

    /* save natural coordinates */
    rst[0] = r;
    rst[1] = s;
E3D(rst[2] = t)

    /* check bounds */
    if(fabs(r) - 1.0 > tol) return 0;
    if(fabs(s) - 1.0 > tol) return 0;
E3D(if(fabs(t) - 1.0 > tol) return 0)

    /* found point in quadrant */
    return 1;
}

int p4est_utilities_unst_pt_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                        p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                        void *pt_in){
    pt2_t *pt = (pt2_t *) pt_in;
    quad_data_t *data;

    Real rst_out[3] = {0.0,0.0,0.0};
    Real xyz[3] = {0.0,0.0,0.0};
    int found;
    int i;

    /* suppress compiler warning */
    (void) p4est;
    (void) which_tree;

    /* data held by quadrant */
    data = (quad_data_t *) quad->p.user_data;

    /* copy into 3D coordinate vector */
    for(i = 0; i < DIM; i++) xyz[i] = pt->x[i];

    /* search point */
    //TODO
    exit(0);
    //dg4est_solver_donor_inclusion_test(xyz,data->pdegree,data->geom,&found,rst_out);

    /* return if not found */
    if(!found) return 0;

    /* if not a leaf, return true */
    if(local_num == -1) return 1;

    /* hit and is a leaf */
    pt->quad_id = local_num;
    return 1;
}

int p4est_utilities_structured_pt_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                              p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                              void *pt_in){
    pt2_t *pt = (pt2_t *) pt_in;
    Real xyz[3*P4EST_CHILDREN];
    int display_quad = 0;
    int i;

    Real x,y,z;
    Real vmin;
    Real vmax;
    box_t box1;

    /* fill in box using axis-aligned bounding volume of this unstructured quadrant */
    p4est_utilities_quad_coordinates(p4est,which_tree,quad,xyz,display_quad);

    /* calculate bounding volume: x */
    vmin = xyz[0];
    vmax = xyz[0];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        x = xyz[3*i+0];
        vmin = (vmin < x) ? vmin:x;
        vmax = (vmax > x) ? vmax:x;
    }
    box1.x.lo = vmin;
    box1.x.hi = vmax;

    /* calculate bounding volume: y */
    vmin = xyz[1];
    vmax = xyz[1];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
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
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        z = xyz[3*i+2];
        vmin = (vmin < z) ? vmin:z;
        vmax = (vmax > z) ? vmax:z;
    }
    box1.z.lo = vmin;
    box1.z.hi = vmax;
#endif

    /* check intersection with axis-aligned bounding box point */
    if(pt->x[0] < box1.x.lo) return 0;
    if(pt->x[0] > box1.x.hi) return 0;
    if(pt->x[1] < box1.y.lo) return 0;
    if(pt->x[1] > box1.y.hi) return 0;
E3D(if(pt->x[2] < box1.z.lo) return 0)
E3D(if(pt->x[2] > box1.z.hi) return 0)

    /* hit but not a leaf */
    if(local_num == -1) return 1;

    /* hit and is a leaf */
    pt->quad_id = local_num;
    return 1;
}

int p4est_utilities_unst_bbox_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                          p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                          void *box_in){
    bbox_t *box = (bbox_t *) box_in;
    box_t box1;
    box_t box2;

    Real xyz[3*P4EST_CHILDREN];
    int display_quad = 0;
    int i;

    Real x,y,z;
    Real vmin;
    Real vmax;
    char hit;

    /* fill in box1 from bbox: need to sort to make box is oriented */
    vmin = (box->bbox_lo[0] < box->bbox_hi[0]) ? box->bbox_lo[0]:box->bbox_hi[0];
    vmax = (box->bbox_lo[0] < box->bbox_hi[0]) ? box->bbox_hi[0]:box->bbox_lo[0];
    box1.x.lo = vmin;
    box1.x.hi = vmax;

    vmin = (box->bbox_lo[1] < box->bbox_hi[1]) ? box->bbox_lo[1]:box->bbox_hi[1];
    vmax = (box->bbox_lo[1] < box->bbox_hi[1]) ? box->bbox_hi[1]:box->bbox_lo[1];
    box1.y.lo = vmin;
    box1.y.hi = vmax;

E3D(vmin = (box->bbox_lo[2] < box->bbox_hi[2]) ? box->bbox_lo[2]:box->bbox_hi[2])
E3D(vmax = (box->bbox_lo[2] < box->bbox_hi[2]) ? box->bbox_hi[2]:box->bbox_lo[2])
E3D(box1.z.lo = vmin)
E3D(box1.z.hi = vmax)

    /* fill in box2 using axis-aligned bounding volume of this unstructured quadrant */
    p4est_utilities_quad_coordinates(p4est,which_tree,quad,xyz,display_quad);

    /* calculate bounding volume: x */
    vmin = xyz[0];
    vmax = xyz[0];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        x = xyz[3*i+0];
        vmin = (vmin < x) ? vmin:x;
        vmax = (vmax > x) ? vmax:x;
    }
    box2.x.lo = vmin;
    box2.x.hi = vmax;

    /* calculate bounding volume: y */
    vmin = xyz[1];
    vmax = xyz[1];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        y = xyz[3*i+1];
        vmin = (vmin < y) ? vmin:y;
        vmax = (vmax > y) ? vmax:y;
    }
    box2.y.lo = vmin;
    box2.y.hi = vmax;

    /* calculate bounding volume: z */
#ifdef P4_TO_P8
    vmin = xyz[2];
    vmax = xyz[2];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        z = xyz[3*i+2];
        vmin = (vmin < z) ? vmin:z;
        vmax = (vmax > z) ? vmax:z;
    }
    box2.z.lo = vmin;
    box2.z.hi = vmax;
#endif

    /* check intersections of the boxes */
    hit = overlapping1D(box1.x,box2.x) && overlapping1D(box1.y,box2.y);
E3D(hit = hit && overlapping1D(box1.z,box2.z))

    if (hit) {
        /* hit but not a leaf */
        if(local_num == -1) return 1;

        /* hit and is a leaf */
        box->quad_hit[local_num] = 1;
        return 1;
    }
    /* not hit */
    return 0;
}

int p4est_utilities_unst_bbox2_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                           p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                           void *box_in){
    bbox2_t *box = (bbox2_t *) box_in;
    box_t box1;
    box_t box2;

    Real xyz[3*P4EST_CHILDREN];
    int display_quad = 0;
    int i;

    Real x,y,z;
    Real vmin;
    Real vmax;
    char hit;

    /* fill in box1 from bbox: need to sort to make box is oriented */
    vmin = (box->bbox_lo[0] < box->bbox_hi[0]) ? box->bbox_lo[0]:box->bbox_hi[0];
    vmax = (box->bbox_lo[0] < box->bbox_hi[0]) ? box->bbox_hi[0]:box->bbox_lo[0];
    box1.x.lo = vmin;
    box1.x.hi = vmax;

    vmin = (box->bbox_lo[1] < box->bbox_hi[1]) ? box->bbox_lo[1]:box->bbox_hi[1];
    vmax = (box->bbox_lo[1] < box->bbox_hi[1]) ? box->bbox_hi[1]:box->bbox_lo[1];
    box1.y.lo = vmin;
    box1.y.hi = vmax;

E3D(vmin = (box->bbox_lo[2] < box->bbox_hi[2]) ? box->bbox_lo[2]:box->bbox_hi[2])
E3D(vmax = (box->bbox_lo[2] < box->bbox_hi[2]) ? box->bbox_hi[2]:box->bbox_lo[2])
E3D(box1.z.lo = vmin)
E3D(box1.z.hi = vmax)

    /* fill in box2 using axis-aligned bounding volume of this unstructured quadrant */
    p4est_utilities_quad_coordinates(p4est,which_tree,quad,xyz,display_quad);

    /* calculate bounding volume: x */
    vmin = xyz[0];
    vmax = xyz[0];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        x = xyz[3*i+0];
        vmin = (vmin < x) ? vmin:x;
        vmax = (vmax > x) ? vmax:x;
    }
    box2.x.lo = vmin;
    box2.x.hi = vmax;

    /* calculate bounding volume: y */
    vmin = xyz[1];
    vmax = xyz[1];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        y = xyz[3*i+1];
        vmin = (vmin < y) ? vmin:y;
        vmax = (vmax > y) ? vmax:y;
    }
    box2.y.lo = vmin;
    box2.y.hi = vmax;

    /* calculate bounding volume: z */
#ifdef P4_TO_P8
    vmin = xyz[2];
    vmax = xyz[2];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        z = xyz[3*i+2];
        vmin = (vmin < z) ? vmin:z;
        vmax = (vmax > z) ? vmax:z;
    }
    box2.z.lo = vmin;
    box2.z.hi = vmax;
#endif

    /* check intersections of the boxes */
    hit = overlapping1D(box1.x,box2.x) && overlapping1D(box1.y,box2.y);
E3D(hit = hit && overlapping1D(box1.z,box2.z))

    if (hit) {
        /* hit but not a leaf */
        if(local_num == -1) return 1;

        /* hit and is a leaf */
        box->box_hit = 1;
        return 1;
    }
    /* not hit */
    return 0;
}

int p4est_utilities_structured_bbox_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                                p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                                void *box_in){
    dg4est_t *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t       *ctx = &dg4est->d_ctx;
    grid_t     *grid = &ctx->d_grid;
    bbox2_t     *box = (bbox2_t *) box_in;

    Real xlo[3],xhi[3];
    Real dx,dy,dz;

    /* quadrant mesh size */
    dx = grid->max_dx[0] * (Real) P4EST_QUADRANT_LEN (quad->level) / (Real) P4EST_ROOT_LEN;
    dy = grid->max_dx[1] * (Real) P4EST_QUADRANT_LEN (quad->level) / (Real) P4EST_ROOT_LEN;
    dz = grid->max_dx[2] * (Real) P4EST_QUADRANT_LEN (quad->level) / (Real) P4EST_ROOT_LEN;

    /* lower corner of quadrant */
    p4est_utilities_get_xlo(p4est,which_tree,quad,xlo);

    /* upper corner of quadrant */
    xhi[0] = xlo[0] + dx;
    xhi[1] = xlo[1] + dy;
    xhi[2] = xlo[2] + dz;

    if(xhi[0] < box->bbox_lo[0]) return 0;
    if(xlo[0] > box->bbox_hi[0]) return 0;
    if(xhi[1] < box->bbox_lo[1]) return 0;
    if(xlo[1] > box->bbox_hi[1]) return 0;
E3D(if(xhi[2] < box->bbox_lo[2]) return 0)
E3D(if(xlo[2] > box->bbox_hi[2]) return 0)

    /* hit but not a leaf */
    if(local_num == -1) return 1;

    /* hit and is a leaf */
    box->box_hit = 1;
    return 1;
}

char inside(bound_t b,bound_t q){
    return (b.lo <= q.lo && q.hi <= b.hi);
}

int count_bits(uint8_t b){
    b = b - ((b >> 1) & 0x55);
    b = (b & 0x33) + ((b >> 2) & 0x33);
    return (((b + (b >> 4)) & 0x0F) * 0x01);
}

void p4est_box_intersection_mask_set(box_t bbox,box_t qbox,char *mask){
    /* =================================================================== */
    /* qbox intersects or exists within bbox. We want to determine if qbox */
    /* intersects the surfaces of bbox. If so, what side of qbox does it   */
    /* intersect bbox (outward normal side from bbox).                     */
    /* Note: this only works for p4est brick.                              */
    /* =================================================================== */
    char inx = inside(bbox.x,qbox.x); // for checking corners
    char iny = inside(bbox.y,qbox.y); // for checking corners
    char inz = (DIM==3) ? inside(bbox.z,qbox.z) : 1; // for checking corners

    if(qbox.x.lo < bbox.x.lo && qbox.x.hi >= bbox.x.lo && ((qbox.x.lo == bbox.x.lo) ? (iny && inz):1)) mask[0] |= xlo_mask; // xlo intersect
    if(qbox.x.hi > bbox.x.hi && qbox.x.lo <= bbox.x.hi && ((qbox.x.hi == bbox.x.hi) ? (iny && inz):1)) mask[0] |= xhi_mask; // xhi intersect
    if(qbox.y.lo < bbox.y.lo && qbox.y.hi >= bbox.y.lo && ((qbox.y.lo == bbox.y.lo) ? (inx && inz):1)) mask[0] |= ylo_mask; // ylo intersect
    if(qbox.y.hi > bbox.y.hi && qbox.y.lo <= bbox.y.hi && ((qbox.y.hi == bbox.y.hi) ? (inx && inz):1)) mask[0] |= yhi_mask; // yhi intersect
E3D(if(qbox.z.lo < bbox.z.lo && qbox.z.hi >= bbox.z.lo && ((qbox.z.lo == bbox.z.lo) ? (inx && iny):1)) mask[0] |= zlo_mask) // zlo intersect
E3D(if(qbox.z.hi > bbox.z.hi && qbox.z.lo <= bbox.z.hi && ((qbox.z.hi == bbox.z.hi) ? (inx && iny):1)) mask[0] |= zhi_mask) // zhi intersect
}

int p4est_utilities_bbox_side_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                          p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                          void *box_in){
    bbox_t *box = (bbox_t *) box_in;
    box_t box1;
    box_t box2;

    Real xyz[3*P4EST_CHILDREN];
    int display_quad = 0;
    int i;

    Real x,y,z;
    Real vmin;
    Real vmax;
    char hit;

    /* fill in box1 from bbox: need to sort to make box is oriented */
    vmin = (box->bbox_lo[0] < box->bbox_hi[0]) ? box->bbox_lo[0]:box->bbox_hi[0];
    vmax = (box->bbox_lo[0] < box->bbox_hi[0]) ? box->bbox_hi[0]:box->bbox_lo[0];
    box1.x.lo = vmin;
    box1.x.hi = vmax;

    vmin = (box->bbox_lo[1] < box->bbox_hi[1]) ? box->bbox_lo[1]:box->bbox_hi[1];
    vmax = (box->bbox_lo[1] < box->bbox_hi[1]) ? box->bbox_hi[1]:box->bbox_lo[1];
    box1.y.lo = vmin;
    box1.y.hi = vmax;

E3D(vmin = (box->bbox_lo[2] < box->bbox_hi[2]) ? box->bbox_lo[2]:box->bbox_hi[2])
E3D(vmax = (box->bbox_lo[2] < box->bbox_hi[2]) ? box->bbox_hi[2]:box->bbox_lo[2])
E3D(box1.z.lo = vmin)
E3D(box1.z.hi = vmax)

    /* fill in box2 using axis-aligned bounding volume of this unstructured quadrant */
    p4est_utilities_quad_coordinates(p4est,which_tree,quad,xyz,display_quad);

    /* calculate bounding volume: x */
    vmin = xyz[0];
    vmax = xyz[0];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        x = xyz[3*i+0];
        vmin = (vmin < x) ? vmin:x;
        vmax = (vmax > x) ? vmax:x;
    }
    box2.x.lo = vmin;
    box2.x.hi = vmax;

    /* calculate bounding volume: y */
    vmin = xyz[1];
    vmax = xyz[1];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        y = xyz[3*i+1];
        vmin = (vmin < y) ? vmin:y;
        vmax = (vmax > y) ? vmax:y;
    }
    box2.y.lo = vmin;
    box2.y.hi = vmax;

    /* calculate bounding volume: z */
#ifdef P4_TO_P8
    vmin = xyz[2];
    vmax = xyz[2];
    for (i = 1; i < P4EST_CHILDREN; ++i) {
        z = xyz[3*i+2];
        vmin = (vmin < z) ? vmin:z;
        vmax = (vmax > z) ? vmax:z;
    }
    box2.z.lo = vmin;
    box2.z.hi = vmax;
#endif

    /* check intersections of the boxes */
    hit = overlapping1D(box1.x,box2.x) && overlapping1D(box1.y,box2.y);
E3D(hit = hit && overlapping1D(box1.z,box2.z))

    if (hit) {
        /* hit but not a leaf */
        if(local_num == -1) return 1;

        /* hit and is a leaf */
        p4est_box_intersection_mask_set(box1,box2,&box->quad_hit[local_num]);
        return 1;
    }
    /* not hit */
    return 0;
}

static
int p4est_utilities_pt_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                              p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                              void *pt_in){
    dg4est_t *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t       *ctx = &dg4est->d_ctx;
    grid_t     *grid = &ctx->d_grid;
    pt_t         *pt = (pt_t *) pt_in;

    Real dx[P4EST_DIM];
    Real xlo[3];

    /* calculate the length ratio between this level and level 0 */
    const Real base_length =
        (Real) P4EST_QUADRANT_LEN (quad->level) /
        (Real) P4EST_ROOT_LEN;

    /* quadrant mesh size */
    dx[0] = grid->max_dx[0] * base_length;
    dx[1] = grid->max_dx[1] * base_length;
E3D(dx[2] = grid->max_dx[2] * base_length)

    /* get lower coordinate of this quadrant */
    p4est_utilities_get_xlo(p4est,which_tree,quad,xlo);

    /* check intersection with box around each point */
    if(pt->x[0] < xlo[0])       return 0;
    if(pt->x[0] > xlo[0]+dx[0]) return 0;
    if(pt->x[1] < xlo[1])       return 0;
    if(pt->x[1] > xlo[1]+dx[1]) return 0;
E3D(if(pt->x[2] < xlo[2])       return 0)
E3D(if(pt->x[2] > xlo[2]+dx[2]) return 0)

    /* hit but not a leaf */
    if(local_num == -1) return 1;

    /* hit and is a leaf */
    pt->quad_hit[local_num] = 1;
    return 1;
}

static
int p4est_utilities_ray_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                    p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                    void *ray_in){
    dg4est_t *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t       *ctx = &dg4est->d_ctx;
    grid_t     *grid = &ctx->d_grid;
    ray_t       *ray = (ray_t *) ray_in;

    Real dx[P4EST_DIM];
    Real xlo[3];

    vec3 d_origin;
    vec3 d_dir;
    box bb;
    char hit;
    int i;

    /* calculate the length ratio between this level and level 0 */
    const Real base_length =
        (Real) P4EST_QUADRANT_LEN (quad->level) /
        (Real) P4EST_ROOT_LEN;

    /* quadrant mesh size */
    // FIXME: this doesn't work for unstructured grids
    dx[0] = grid->max_dx[0] * base_length;
    dx[1] = grid->max_dx[1] * base_length;
E3D(dx[2] = grid->max_dx[2] * base_length)

    /* get lower coordinate of this quadrant */
    p4est_utilities_get_xlo(p4est,which_tree,quad,xlo);

    /* fill in intersect_box variables */
    for (i = 0; i < P4EST_DIM; ++i) {
        bb.min.x[i]   = xlo[i];
        bb.max.x[i]   = xlo[i] + dx[i];
        d_origin.x[i] = ray->ray_origin[i];
        d_dir.x[i]    = ray->ray_end[i];
    }

    hit = intersect_box(d_origin,d_dir,bb);

    if (hit) {
        /* hit but not a leaf */
        if(local_num == -1) return 1;

        /* hit and is a leaf */
        ray->quad_hit[local_num] = 1;
        return 1;
    }
    /* not hit */
    return 0;
}

quad_list_t * p4est_utilities_points_search(p4est_t *p4est,int npts,Real *xyz){
    quad_list_t *quad_id_list;
    sc_array_t *pt_data;
    pt_t *pt;

    int nquads = p4est->local_num_quadrants;
    int nquads_hit;
    int counter;
    int i,j;

    /* allocate and initialize ray quadrant hit values */
    char *quad_hit = (char *) malloc(nquads*sizeof(char));
    for(i = 0; i < nquads; ++i){
        quad_hit[i] = 0;
    }

    /* allocate the ray structure */
    pt_data = sc_array_new_size(sizeof(pt_t),npts);

    /* fill in ray data */
    for (j = 0; j < npts; ++j) {
        pt = (pt_t *) sc_array_index(pt_data,j);
        for (i = 0; i < P4EST_DIM; ++i) {
            pt->x[i] = xyz[P4EST_DIM*j + i];
        }
        /* assign the same quad_hit pointer for single list (no plurality) */
        pt->quad_hit = quad_hit;
    }

    /* search the forest of trees on this mpi rank for the points */
    p4est_search(p4est,NULL,p4est_utilities_pt_search_func,pt_data);

    /* count the number of quadrants hit by ray */
    nquads_hit = 0;
    for (i = 0; i < nquads; ++i) {
        if(quad_hit[i]) nquads_hit++;
    }

    /* allocate the return array of quad _ids */
    quad_id_list = (quad_list_t *) malloc(sizeof(quad_list_t));
    quad_id_list->nquads = nquads_hit;
    quad_id_list->quad_list = (int *) malloc(nquads_hit*sizeof(int));

    /* copy quad_id_list to the return values */
    counter = 0;
    for (i = 0; i < nquads; ++i) {
        if(quad_hit[i]) quad_id_list->quad_list[counter++] = i;
    }

    /* deallocate pt_data */
    free(quad_hit);
    sc_array_destroy(pt_data);

    /* return the quadrants list */
    return quad_id_list;
}

quad_list_t * p4est_utilities_ray_search(p4est_t *p4est,Real *pt_a,Real *pt_b){
    quad_list_t *quad_id_list;
    sc_array_t *pt_data;
    ray_t *pt;

    int nquads = p4est->local_num_quadrants;
    int nquads_hit;
    int counter;
    int i;

    /* allocate the ray structure */
    pt_data = sc_array_new_size(sizeof(ray_t),1);

    /* initialize the ray data from user inputs */
    pt = (ray_t *) sc_array_index(pt_data,0);
    for (i = 0; i < P4EST_DIM; ++i) {
        pt->ray_origin[i] = pt_a[i];
        pt->ray_end[i]    = pt_b[i];
    }

    /* allocate  ray quadrant hist list */
    pt->quad_hit = (char *) malloc(nquads*sizeof(char));

    /* initialize ray quadrant hit values */
    for (i = 0; i < nquads; ++i) {
        pt->quad_hit[i] = 0;
    }

    /* search the forest of trees on this mpi rank for the points */
    p4est_search(p4est,NULL,p4est_utilities_ray_search_func,pt_data);

    /* count the number of quadrants hit by ray */
    nquads_hit = 0;
    for (i = 0; i < nquads; ++i) {
        if(pt->quad_hit[i]) nquads_hit++;
    }

    /* allocate the return array of quad_ids */
    quad_id_list = (quad_list_t *) malloc(sizeof(quad_list_t));
    quad_id_list->nquads = nquads_hit;
    quad_id_list->quad_list = (int *) malloc(nquads_hit*sizeof(int));

    /* copy quad_id_list to the return values */
    counter = 0;
    for (i = 0; i < nquads; ++i) {
        if(pt->quad_hit[i]) quad_id_list->quad_list[counter++] = i;
    }

    /* deallocate pt_data */
    free(pt->quad_hit);
    sc_array_destroy(pt_data);

    /* return the quadrants list */
    return quad_id_list;
}

quad_list_t * p4est_utilities_bbox_search(p4est_t *p4est,Real *bbox_lo,Real *bbox_hi){
    quad_list_t *quad_id_list;
    bbox_t *bb;
    sc_array_t *bb_data;

    int nquads = p4est->local_num_quadrants;
    int nquads_hit;
    int counter;
    int i;

    /* allocate the ray structure */
    bb_data = sc_array_new_size(sizeof(bbox_t),1);

    /* initialize the ray data from user inputs */
    bb = (bbox_t *) sc_array_index(bb_data,0);
    for (i = 0; i < P4EST_DIM; ++i) {
        bb->bbox_lo[i] = bbox_lo[i];
        bb->bbox_hi[i] = bbox_hi[i];
    }

    /* allocate quadrant hit list */
    bb->quad_hit = (char *) malloc(nquads*sizeof(char));

    /* initialize bbox quadrant hit values */
    for(i = 0; i < nquads; ++i) bb->quad_hit[i] = 0;

    /* search the forest of trees on this mpi rank for bounding box */
    p4est_search(p4est,NULL,p4est_utilities_unst_bbox_search_func,bb_data);

    /* count the number of quadrants intersecting the bounding box */
    nquads_hit = 0;
    for (i = 0; i < nquads; ++i) {
        if(bb->quad_hit[i]) nquads_hit++;
    }

    /* allocate the return array of quad_ids */
    quad_id_list = (quad_list_t *) malloc(sizeof(quad_list_t));
    quad_id_list->nquads = nquads_hit;
    quad_id_list->quad_list = (int *) malloc(nquads_hit*sizeof(int));

    /* copy quad_id_list to the return values */
    counter = 0;
    for (i = 0; i < nquads; ++i) {
        if(bb->quad_hit[i]) quad_id_list->quad_list[counter++] = i;
    }

    /* deallocate pt_data */
    free(bb->quad_hit);
    sc_array_destroy(bb_data);

    /* return the quadrants list */
    return quad_id_list;
}

quad_side_list_t * p4est_utilities_bbox_surface_search(p4est_t *p4est,Real *bbox_lo,Real *bbox_hi){
    quad_side_list_t *quad_id_list;
    bbox_t *bb;
    sc_array_t *bb_data;

    int nquads = p4est->local_num_quadrants;
    int nquads_hit;
    int counter;
    int i;

    /**< boundary bit mask:
     *     0-interior of domain
     *     1-intersects boundary
     * 8 bits:
     *     bit 0: xlo face  bit 1: xhi face
     *     bit 2: ylo face  bit 3: yhi face
     *     bit 4: zlo face  bit 5: zhi face
     *     bit 6: free      bit 7: free
     */

    /* allocate the ray structure */
    bb_data = sc_array_new_size(sizeof(bbox_t),1);

    /* initialize the ray data from user inputs */
    bb = (bbox_t *) sc_array_index(bb_data,0);
    for (i = 0; i < P4EST_DIM; ++i) {
        bb->bbox_lo[i] = bbox_lo[i];
        bb->bbox_hi[i] = bbox_hi[i];
    }

    /* allocate quadrant hit list */
    bb->quad_hit = (char *) malloc(nquads*sizeof(char));

    /* initialize bbox quadrant hit bit-mask values */
    for(i = 0; i < nquads; ++i) bb->quad_hit[i] = 0;

    /* search the forest of trees on this mpi rank for bounding box */
    p4est_search(p4est,NULL,p4est_utilities_bbox_side_search_func,bb_data);

    /* count the number of quadrants intersecting the bounding box */
    nquads_hit = 0;
    for (i = 0; i < nquads; ++i) {
        if(bb->quad_hit[i]) nquads_hit++;
    }

    /* allocate the return array of quad_ids */
    quad_id_list = (quad_side_list_t *) malloc(sizeof(quad_side_list_t));
    quad_id_list->nquads = nquads_hit;
    quad_id_list->quad_list = (int *) malloc(nquads_hit*sizeof(int));
    quad_id_list->side_list = (char *) malloc(nquads_hit*sizeof(char));
    quad_id_list->nsides = 0;

    /* copy quad_id_list to the return values */
    counter = 0;
    for (i = 0; i < nquads; ++i) {
        if (bb->quad_hit[i]) {
            quad_id_list->quad_list[counter] = i;
            quad_id_list->side_list[counter] = bb->quad_hit[i];
            quad_id_list->nsides += count_bits(bb->quad_hit[i]);
            counter++;
        }
    }

    /* deallocate pt_data */
    free(bb->quad_hit);
    sc_array_destroy(bb_data);

    /* return the quadrants list */
    return quad_id_list;
}

quad_side_list_t * p4est_utilities_surface_name_list(p4est_t *p4est,
                                                     geometry_t *geom,
                                                     char *surface_patch_name){
    quad_side_list_t *quad_id_list;

    p4est_topidx_t   first_local_tree = p4est->first_local_tree;
    p4est_topidx_t   last_local_tree  = p4est->last_local_tree;
    p4est_topidx_t   t;
    sc_array_t       *trees = p4est->trees;
    p4est_tree_t     *tree;
    sc_array_t       *quadrants;
    p4est_quadrant_t *quad;
    quad_data_t      *data;
    p4est_locidx_t    si,n_quads;

    int nquads = p4est->local_num_quadrants;
    int nquads_hit;
    int counter;
    int i,j;

    /**< boundary bit mask:
     *     0-interior of domain
     *     1-intersects boundary
     * 8 bits:
     *     bit 0: xlo face  bit 1: xhi face
     *     bit 2: ylo face  bit 3: yhi face
     *     bit 4: zlo face  bit 5: zhi face
     *     bit 6: free      bit 7: free
     */
    char bitmasks[6] = {xlo_mask,xhi_mask,
                        ylo_mask,yhi_mask,
                        zlo_mask,zhi_mask};

    /* allocate quadrant hit list */
    char *quad_hit = (char *) malloc(nquads*sizeof(char));

    /* initialize quadrant hit bit-mask values */
    for(i = 0; i < nquads; ++i) quad_hit[i] = 0;

    /* search the forest of trees on this mpi rank for quadrants touching surface */
    counter = 0;
    for (t = first_local_tree; t <= last_local_tree; ++t) {
        tree = p4est_tree_array_index(trees,t);
        quadrants = &(tree->quadrants);
        n_quads = (p4est_locidx_t) quadrants->elem_count;

        for (si = 0; si < n_quads; ++si) {
            quad = p4est_quadrant_array_index(quadrants,si);
            data = (quad_data_t *) quad->p.user_data;

            /* check if this quadrant touches the surface patch quadrant list */
            for(j = 0; j < 2*DIM; j++){
                if(data->bc[j] >= 0 && strstr(surface_patch_name,geom->patch_name[data->bc[j]])){
                    quad_hit[counter] |= bitmasks[j];
                }
            }
            counter++;
        }
    }

    /* count the number of quadrants intersecting the bounding box */
    nquads_hit = 0;
    for (i = 0; i < nquads; ++i) {
        if(quad_hit[i]) nquads_hit++;
    }

    /* allocate the return array of quad_ids */
    quad_id_list = (quad_side_list_t *) malloc(sizeof(quad_side_list_t));
    quad_id_list->nquads = nquads_hit;
    quad_id_list->quad_list = (int *) malloc(nquads_hit*sizeof(int));
    quad_id_list->side_list = (char *) malloc(nquads_hit*sizeof(char));
    quad_id_list->nsides = 0;

    /* copy quad_id_list to the return values */
    counter = 0;
    for (i = 0; i < nquads; ++i) {
        if (quad_hit[i]) {
            quad_id_list->quad_list[counter] = i;
            quad_id_list->side_list[counter] = quad_hit[i];
            quad_id_list->nsides += count_bits(quad_hit[i]);
            counter++;
        }
    }

    /* deallocate local memory */
    free(quad_hit);

    /* return the quadrants list */
    return quad_id_list;
}

quad_list_t * p4est_utilities_volume_name_list(p4est_t *p4est,
                                               geometry_t *geom,
                                               char *volume_patch_names){
    quad_list_t *quad_id_list;

    p4est_topidx_t   first_local_tree = p4est->first_local_tree;
    p4est_topidx_t   last_local_tree  = p4est->last_local_tree;
    p4est_topidx_t   t;
    sc_array_t       *trees = p4est->trees;
    p4est_tree_t     *tree;
    sc_array_t       *quadrants;
    p4est_quadrant_t *quad;
    quad_data_t      *data;
    p4est_locidx_t    si,n_quads;

    int nquads = p4est->local_num_quadrants;
    int nquads_hit;
    int counter;
    int i;

    /* allocate quadrant hit list */
    char *quad_hit = (char *) malloc(nquads*sizeof(char));

    /* initialize quadrant hit values */
    for(i = 0; i < nquads; ++i) quad_hit[i] = 0;

    nquads_hit = 0;
    for (t = first_local_tree; t <= last_local_tree; ++t) {
        tree = p4est_tree_array_index(trees,t);
        quadrants = &(tree->quadrants);
        n_quads = (p4est_locidx_t) quadrants->elem_count;

        for (si = 0; si < n_quads; ++si) {
            quad = p4est_quadrant_array_index(quadrants,si);
            data = (quad_data_t *) quad->p.user_data;

            /* check if this quadrant is in the volume patch quadrant list */
            if(data->type >= 0 && strstr(volume_patch_names,geom->patch_name[data->type])){
                quad_hit[nquads_hit++] = 1;
            }
        }
    }

    /* allocate the return array of quad_ids */
    quad_id_list = (quad_list_t *) malloc(sizeof(quad_list_t));
    quad_id_list->nquads = nquads_hit;
    quad_id_list->quad_list = (int *) malloc(nquads_hit*sizeof(int));

    /* copy quad_id_list to the return values */
    counter = 0;
    for (i = 0; i < nquads; ++i) {
        if(quad_hit[i]) quad_id_list->quad_list[counter++] = i;
    }

    /* deallocate local memory */
    free(quad_hit);

    /* return the quadrants list */
    return quad_id_list;
}
