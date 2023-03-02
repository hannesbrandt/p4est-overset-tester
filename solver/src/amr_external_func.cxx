/**
 * \file    amr_external_func.cxx
 * \author  akirby
 *
 * \brief   External wrapper functions to call solver_interface code module.
 */

/* header files */
#include "amr_external_func.h"

#define GEOM_IND(jmp,x,y,z) \
    ((z)*(jmp)*(jmp) +      \
     (y)*(jmp)       +      \
     (x)             )

static
void copy_coordinates(Real *dst,Real *src){
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
}

static
void fill_geom_corner(quad_data_t *data,Real *node_geom,int corner){
    int qtm = data->qdegree + 1;
    int node;

    switch (corner) {
        case 0: node = GEOM_IND(qtm,            0,            0,            0); break;
        case 1: node = GEOM_IND(qtm,data->qdegree,            0,            0); break;
        case 2: node = GEOM_IND(qtm,            0,data->qdegree,            0); break;
        case 3: node = GEOM_IND(qtm,data->qdegree,data->qdegree,            0); break;
        case 4: node = GEOM_IND(qtm,            0,            0,data->qdegree); break;
        case 5: node = GEOM_IND(qtm,data->qdegree,            0,data->qdegree); break;
        case 6: node = GEOM_IND(qtm,            0,data->qdegree,data->qdegree); break;
        case 7: node = GEOM_IND(qtm,data->qdegree,data->qdegree,data->qdegree); break;
    }
    copy_coordinates(node_geom,&data->geom[3*node]);
}

void setup_wake3d_allnodes(p4est_t *p4est,ctx_t *ctx){
    wake3d_t   *wke  = &ctx->d_wake3d;
    geometry_t *geom = &ctx->d_geometry;

    p4est_topidx_t first_local_tree = p4est->first_local_tree;
    p4est_topidx_t last_local_tree = p4est->last_local_tree;
    p4est_tree_t *tree;
    sc_array_t *quadrants;
    p4est_quadrant_t *q;
    p4est_nodes_t *nodes;
    sc_array_t *indeps;

    int lqid,jt,zz,corner;
    int numz_quadrants;
    int indep_node_id;
    int i,j;
    const int display = 0;

    /* get quadrant nodes for wake3d data */
    /*
     * NOTE:
     *  * nodes->local_nodes is a list of node indices (orders as local quadrants and corners)
     *    that locate than local node into indeps->array
     *  * nodes->local_nodes contains only grid support nodes (not high-order nodes)
     *  * indeps->array contains a list of UNIQUE local nodes
     */

    nodes = p4est_nodes_new(p4est,NULL);
    indeps = &nodes->indep_nodes;

    dg::memory<int> &ndc = (DIM==2) ? wke->ndcQuad:wke->ndcHex;
    wke->nnodes = nodes->num_owned_indeps;
    wke->iblank.malloc(nodes->num_owned_indeps);
    wke->iblank_cell.malloc(p4est->local_num_quadrants);
    wke->xgeom.malloc(3*wke->nnodes);
    ndc.malloc(P4EST_CHILDREN*p4est->local_num_quadrants);

    /* local independent nodes tag: tag=0 (i.e. "not processed")*/
    std::vector<char> ln(wke->nnodes,0);

    /* initialize wake3d values */
    for(i = 0; i < p4est->local_num_quadrants; ++i) wke->iblank_cell[i] = field_mask;
    for(i = 0; i < nodes->num_owned_indeps; ++i) wke->iblank[i] = field_mask;

    /* note: order of nodes modified for Tioga */
    for (i = 0; i < nodes->num_local_quadrants; ++i) {
        const int ind = i*P4EST_CHILDREN;

        ndc[ind+0] = nodes->local_nodes[ind+0] + BASE;
        ndc[ind+1] = nodes->local_nodes[ind+1] + BASE;
        ndc[ind+2] = nodes->local_nodes[ind+2] + BASE;
        ndc[ind+3] = nodes->local_nodes[ind+3] + BASE;
#ifdef P4_TO_P8
        ndc[ind+4] = nodes->local_nodes[ind+4] + BASE;
        ndc[ind+5] = nodes->local_nodes[ind+5] + BASE;
        ndc[ind+6] = nodes->local_nodes[ind+6] + BASE;
        ndc[ind+7] = nodes->local_nodes[ind+7] + BASE;
#endif
    }

    /* unique grid node coordinates */
    for (jt = first_local_tree, lqid = 0; jt <= last_local_tree; ++jt) {
        tree = p4est_tree_array_index (p4est->trees, jt);
        quadrants = &tree->quadrants;
        numz_quadrants = quadrants->elem_count;

        /* loop through quadrants */
        for (zz = 0; zz < numz_quadrants; ++zz, ++lqid) {
            q = p4est_quadrant_array_index (quadrants, zz);
            quad_data_t *data = (quad_data_t *) q->p.user_data;

            /* loop through corners (real element vertices) */
            for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
                /* get independent node */
                indep_node_id = nodes->local_nodes[lqid*P4EST_CHILDREN+corner];

                /* if already processed continue to next node */
                if(ln[indep_node_id]) continue;

                /* process corner geometry */
                fill_geom_corner(data,&wke->xgeom[3*indep_node_id],corner);

                /* mark the node as processed */
                ln[indep_node_id] = 1;
            }
        }
    }

    /* Print out node numbering */
    if (display) {
        for (i = 0; i < nodes->num_local_quadrants; ++i) {
            printf("------------------------------------\n");

            const int ind = i*P4EST_CHILDREN;
            for (j = 0; j < P4EST_CHILDREN; ++j) {
                Real x = wke->xgeom[3*nodes->local_nodes[ind+j]+0];
                Real y = wke->xgeom[3*nodes->local_nodes[ind+j]+1];
                Real z = (DIM==3) ? wke->xgeom[3*nodes->local_nodes[ind+j]+2]:0.0;;
                printf("NODE [rank%d]: "
                       "Local quadrant %d; "
                       "Corner %d; "
                       "Node unique ID %d, "
                       "Node coords(x,y,z) [",
                       ctx->rank,i,j,
                       nodes->local_nodes[ind+j]);

                (x<0.0) ? printf("%f ",x):printf(" %f ",x);
                (y<0.0) ? printf("%f ",y):printf(" %f ",y);
                (z<0.0) ? printf("%f ",z):printf(" %f ",z);
                printf("]\n");
            }
        }
    }

    /*
     * FIND WALL AND OUTER BOUNDARY NODES
     */
    if (ctx->d_simulation.unstructured_flag) {
        ivector_t wke_wbc_nodes;
        ivector_t wke_obc_nodes;
        p4est_indep_t *in;

        int il,iface,report;

        /* Initialize vector of nodes wall and outer tags */
        ivector_init(&wke_wbc_nodes,wke->nnodes);
        ivector_init(&wke_obc_nodes,wke->nnodes);

        /* initialize the local nodes with tag 0 (i.e. "not processed") */
        for(il = 0; il < wke->nnodes; ++il) ln[il] = 0;

        /* loop over all local trees */
        for (jt = first_local_tree, lqid = 0; jt <= last_local_tree; ++jt) {
            tree = p4est_tree_array_index (p4est->trees, jt);
            quadrants = &tree->quadrants;
            numz_quadrants = quadrants->elem_count;

            /* loop through quadrants */
            for (zz = 0; zz < numz_quadrants; ++zz, ++lqid) {
                q = p4est_quadrant_array_index (quadrants, zz);
                quad_data_t *data = (quad_data_t *) q->p.user_data;

                /* loop through corners (real element vertices) */
                for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
                    int c_wbc = 0;
                    int c_obc = 0;

                    /* get independent node */
                    indep_node_id = nodes->local_nodes[lqid*P4EST_CHILDREN+corner];
                    in = (p4est_indep_t *) sc_array_index(indeps,indep_node_id);

                    /* assert indep_node_id == in->p.piggy3.local_num */
                    if (indep_node_id != in->p.piggy3.local_num) {
                        DGOUT(stderr,
                              "\033[1;31m>>>>> [DG4EST] "
                              "indep_node_id does not match in->p.piggy3.local_num!"
                              "\033[0m (%s:%d)\n",__FILE__,__LINE__);
                        exit(EXIT_FAILURE);
                    }

                    /* if already processed continue to next node */
                    if(ln[indep_node_id]) continue;

                    /* for each quadrant face */
                    if(in->x == 0              && data->bc[0] != NO_PATCH && geom->wake3d_bc[data->bc[0]] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;}
                    if(in->x == P4EST_ROOT_LEN && data->bc[1] != NO_PATCH && geom->wake3d_bc[data->bc[1]] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;}
                    if(in->y == 0              && data->bc[2] != NO_PATCH && geom->wake3d_bc[data->bc[2]] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;}
                    if(in->y == P4EST_ROOT_LEN && data->bc[3] != NO_PATCH && geom->wake3d_bc[data->bc[3]] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;}
                X3D(if(in->z == 0              && data->bc[4] != NO_PATCH && geom->wake3d_bc[data->bc[4]] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;})
                X3D(if(in->z == P4EST_ROOT_LEN && data->bc[5] != NO_PATCH && geom->wake3d_bc[data->bc[5]] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;})

                    if(in->x == 0              && data->bc[0] != NO_PATCH && geom->wake3d_bc[data->bc[0]] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;}
                    if(in->x == P4EST_ROOT_LEN && data->bc[1] != NO_PATCH && geom->wake3d_bc[data->bc[1]] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;}
                    if(in->y == 0              && data->bc[2] != NO_PATCH && geom->wake3d_bc[data->bc[2]] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;}
                    if(in->y == P4EST_ROOT_LEN && data->bc[3] != NO_PATCH && geom->wake3d_bc[data->bc[3]] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;}
                X3D(if(in->z == 0              && data->bc[4] != NO_PATCH && geom->wake3d_bc[data->bc[4]] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;})
                X3D(if(in->z == P4EST_ROOT_LEN && data->bc[5] != NO_PATCH && geom->wake3d_bc[data->bc[5]] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;})

                    /* mark the node as processed */
                    ln[indep_node_id] = 1;

                    if (display && report) {
                        DGOUT(ctx->log_io,
                              "BOUNDARY NODE [rank%d]: "
                              "Tree: %d, "
                              "Quadrant: %d, "
                              "Corner: %d, "
                              "Face: %d, "
                              "Node unique ID %d, "
                              "Wake3D_wbc: %d, "
                              "Wake3D_obc: %d, "
                              "Node_coords: [%.15f,%.15f,%.15f] \n",
                              ctx->rank,jt,zz,corner,iface,
                              indep_node_id,c_wbc,c_obc,
                              wke->xgeom[3*indep_node_id+0],
                              wke->xgeom[3*indep_node_id+1],
                              wke->xgeom[3*indep_node_id+2]);
                        report = 0;
                    }
                } // corner
            } // quadrant
        } // tree

        /* trim arrays */
        ivector_trim(&wke_wbc_nodes);
        ivector_trim(&wke_obc_nodes);

        /* number of wall and outer nodes */
        wke->nwbc = wke_wbc_nodes.size;
        wke->nobc = wke_obc_nodes.size;

        /* allocate wall and outer bc */
        wke->iwbcnode.malloc(wke->nwbc);
        wke->iobcnode.malloc(wke->nobc);

        /* copy array */
        for(il = 0; il < wke->nwbc; il++) wke->iwbcnode[il] = wke_wbc_nodes.array[il];
        for(il = 0; il < wke->nobc; il++) wke->iobcnode[il] = wke_obc_nodes.array[il];

        /* deallocate nodes list */
        ivector_free(&wke_wbc_nodes);
        ivector_free(&wke_obc_nodes);
    } else {
        ivector_t wke_wbc_nodes;
        ivector_t wke_obc_nodes;
        p4est_indep_t *in;

        char filename[] = "input.dg4est.v3";
        int bc[6],wbc[6],obc[6];
        int il;

        /* Initialize vector of nodes wall and outer tags */
        ivector_init(&wke_wbc_nodes,wke->nnodes);
        ivector_init(&wke_obc_nodes,wke->nnodes);

        /* initialize the local nodes with tag 0 (i.e. "not processed") */
        for(il = 0; il < wke->nnodes; ++il) ln[il] = 0;

        amr_utilities_find_keyword_integer(filename,"bc_xlo:",&bc[0]);
        amr_utilities_find_keyword_integer(filename,"bc_xhi:",&bc[1]);
        amr_utilities_find_keyword_integer(filename,"bc_ylo:",&bc[2]);
        amr_utilities_find_keyword_integer(filename,"bc_yhi:",&bc[3]);
        amr_utilities_find_keyword_integer(filename,"bc_zlo:",&bc[4]);
        amr_utilities_find_keyword_integer(filename,"bc_zhi:",&bc[5]);

        /**< Boundary conditions */
#define BC_FREESTREAM     1
#define BC_SLIPWALL       2
#define BC_NONSPLIPWALL   3
#define BC_CHARACTERISTIC 4
#define BC_TOTALINLET     5
#define BC_STATICOUTLET   6

        wbc[0] = (bc[0] == BC_SLIPWALL || bc[0] == BC_NONSPLIPWALL) ? W3D_WBC:W3D_NO;
        wbc[1] = (bc[1] == BC_SLIPWALL || bc[1] == BC_NONSPLIPWALL) ? W3D_WBC:W3D_NO;
        wbc[2] = (bc[2] == BC_SLIPWALL || bc[2] == BC_NONSPLIPWALL) ? W3D_WBC:W3D_NO;
        wbc[3] = (bc[3] == BC_SLIPWALL || bc[3] == BC_NONSPLIPWALL) ? W3D_WBC:W3D_NO;
        wbc[4] = (bc[4] == BC_SLIPWALL || bc[4] == BC_NONSPLIPWALL) ? W3D_WBC:W3D_NO;
        wbc[5] = (bc[5] == BC_SLIPWALL || bc[5] == BC_NONSPLIPWALL) ? W3D_WBC:W3D_NO;

        obc[0] = (bc[0] == BC_CHARACTERISTIC || bc[0] == BC_FREESTREAM) ? W3D_OBC:W3D_NO;
        obc[1] = (bc[1] == BC_CHARACTERISTIC || bc[1] == BC_FREESTREAM) ? W3D_OBC:W3D_NO;
        obc[2] = (bc[2] == BC_CHARACTERISTIC || bc[2] == BC_FREESTREAM) ? W3D_OBC:W3D_NO;
        obc[3] = (bc[3] == BC_CHARACTERISTIC || bc[3] == BC_FREESTREAM) ? W3D_OBC:W3D_NO;
        obc[4] = (bc[4] == BC_CHARACTERISTIC || bc[4] == BC_FREESTREAM) ? W3D_OBC:W3D_NO;
        obc[5] = (bc[5] == BC_CHARACTERISTIC || bc[5] == BC_FREESTREAM) ? W3D_OBC:W3D_NO;

        /* loop over all local trees */
        for (jt = first_local_tree, lqid = 0; jt <= last_local_tree; ++jt) {
            tree = p4est_tree_array_index (p4est->trees, jt);
            quadrants = &tree->quadrants;
            numz_quadrants = quadrants->elem_count;

            /* loop through quadrants */
            for (zz = 0; zz < numz_quadrants; ++zz, ++lqid) {
                q = p4est_quadrant_array_index (quadrants, zz);
                quad_data_t *data = (quad_data_t *) q->p.user_data;

                /* loop through corners (real element vertices) */
                for (corner = 0; corner < P4EST_CHILDREN; ++corner) {
                    int report = 0;
                    int c_wbc = 0;
                    int c_obc = 0;

                    /* get independent node */
                    indep_node_id = nodes->local_nodes[lqid*P4EST_CHILDREN+corner];
                    in = (p4est_indep_t *) sc_array_index(indeps,indep_node_id);

                    /* if already processed continue to next node */
                    if(ln[indep_node_id]) continue;

                    /* for each quadrant face */
                    if(in->x == 0              && data->bc[0] != NO_PATCH && wbc[0] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;}
                    if(in->x == P4EST_ROOT_LEN && data->bc[1] != NO_PATCH && wbc[1] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;}
                    if(in->y == 0              && data->bc[2] != NO_PATCH && wbc[2] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;}
                    if(in->y == P4EST_ROOT_LEN && data->bc[3] != NO_PATCH && wbc[3] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;}
                X3D(if(in->z == 0              && data->bc[4] != NO_PATCH && wbc[4] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;})
                X3D(if(in->z == P4EST_ROOT_LEN && data->bc[5] != NO_PATCH && wbc[5] == W3D_WBC && c_wbc == 0) {ivector_pushback(&wke_wbc_nodes,indep_node_id + BASE); c_wbc = report = 1;})

                    if(in->x == 0              && data->bc[0] != NO_PATCH && obc[0] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;}
                    if(in->x == P4EST_ROOT_LEN && data->bc[1] != NO_PATCH && obc[1] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;}
                    if(in->y == 0              && data->bc[2] != NO_PATCH && obc[2] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;}
                    if(in->y == P4EST_ROOT_LEN && data->bc[3] != NO_PATCH && obc[3] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;}
                X3D(if(in->z == 0              && data->bc[4] != NO_PATCH && obc[4] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;})
                X3D(if(in->z == P4EST_ROOT_LEN && data->bc[5] != NO_PATCH && obc[5] == W3D_OBC && c_obc == 0) {ivector_pushback(&wke_obc_nodes,indep_node_id + BASE); c_obc = report = 1;})

                    /* mark the node as processed */
                    ln[indep_node_id] = 1;

                    if (display && report) {
                        DGOUT(ctx->log_io,
                              "BOUNDARY NODE [rank%d]: "
                              "Tree: %d, "
                              "Quadrant: %d, "
                              "Corner: %d, "
                              "Node unique ID %d, "
                              "Wake3D_wbc: %d, "
                              "Wake3D_obc: %d, "
                              "Node_coords: [%.15f,%.15f,%.15f] \n",
                              ctx->rank,jt,zz,corner,
                              indep_node_id,c_wbc,c_obc,
                              wke->xgeom[3*indep_node_id+0],
                              wke->xgeom[3*indep_node_id+1],
                              wke->xgeom[3*indep_node_id+2]);
                    }
                    report = 0;
                } // corner
            } // quadrant
        } // tree

        /* trim arrays */
        ivector_trim(&wke_wbc_nodes);
        ivector_trim(&wke_obc_nodes);

        /* number of wall and outer nodes */
        wke->nwbc = wke_wbc_nodes.size;
        wke->nobc = wke_obc_nodes.size;

        /* allocate wall and outer bc */
        wke->iwbcnode.malloc(wke->nwbc);
        wke->iobcnode.malloc(wke->nobc);

        /* copy array */
        for(il = 0; il < wke->nwbc; il++) wke->iwbcnode[il] = wke_wbc_nodes.array[il];
        for(il = 0; il < wke->nobc; il++) wke->iobcnode[il] = wke_obc_nodes.array[il];

        /* deallocate nodes list */
        ivector_free(&wke_wbc_nodes);
        ivector_free(&wke_obc_nodes);
    }

    /* free nodes */
    p4est_nodes_destroy(nodes);
}

void external_func_setup_wake3d(p4est_t *p4est,ctx_t *ctx){
    wake3d_t *wke = &ctx->d_wake3d;

    int total_nnodes;
    int total_nobc;
    int total_nwbc;

    /* fill ndc8 and xgeom */
    setup_wake3d_allnodes(p4est,ctx);

    /* count total nobc */
    MPI_Reduce(&wke->nobc,&total_nobc,1,MPI_INT,MPI_SUM,0,ctx->comm);
    MPI_Reduce(&wke->nwbc,&total_nwbc,1,MPI_INT,MPI_SUM,0,ctx->comm);
    MPI_Reduce(&wke->nnodes,&total_nnodes,1,MPI_INT,MPI_SUM,0,ctx->comm);

    if (ctx->rank == 0 && ctx->log_info < DGLOG_STATUS_OFF) {
        DGOUT(ctx->log_io,
              "[ amr ] WAKE3D data structure construction: \n"
              "        DG4EST nnodes: %d,  nobc: %d,  nwbc: %d\n",
              total_nnodes,total_nobc,total_nwbc);
    }
}

void external_func_deallocate_geom(ctx_t *ctx){
    geometry_deallocate_data(ctx);
}