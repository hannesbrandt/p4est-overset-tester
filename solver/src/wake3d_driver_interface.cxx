/**
 * \file    wake3d_driver_interface.c
 * \ingroup amr_group
 * \author  akirby
 *
 * \brief   WAKE3D driver interface functions.
 */

/* header files */
#include "wake3d_driver_interface.h"

dg4est_t dg4est;
p4est_t *p4est;
p4est_connectivity_t *conn;

void driver_interface_initialize_group_mpi(MPI_Fint *group_comm){
    ctx_t *ctx = &dg4est.d_ctx;
    MPI_Comm comm;

    /* default context data */
    ctx_default(ctx);

    /* dg4est file stream for logging */
    ctx->log_io = stdout;//fopen("log.dg4est", "w");

    /* convert Fortran MPI comm to C MPI comm */
    comm = MPI_Comm_f2c(*group_comm);

    /* initialize MPI data */
    mpi_init_from_comm(&comm,ctx,ctx->log_info);
}

void driver_interface_initialize(void){
    ctx_t *ctx = &dg4est.d_ctx;
    simulation_t *sim = &ctx->d_simulation;

    /* input file name */
    char filename[] = "input.dg4est.v3";

    /* initialize dg4est and grid */
    initialize_amr_file(filename,&dg4est,&p4est,&conn);

    /* setup wake3d buffers */
    external_func_setup_wake3d(p4est,ctx);

    if(sim->visualization_interval < INT_MAX) ctx->visualize = 1;
}

int driver_interface_shutdown(void){
    ctx_t *ctx = &dg4est.d_ctx;

    /* deallocate simulation data */
    external_func_deallocate_geom(ctx);

    /* clean up p4est data */
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);
    sc_finalize();

    /* close logging file */
    fclose(ctx->log_io);
    return 0;
}

int driver_interface_mpi_finalize(){
    return mpi_finalize();
}

void driver_interface_output_solution(int *index){
    ctx_t        *ctx = &dg4est.d_ctx;
    simulation_t *sim = &ctx->d_simulation;

    /* move the wake3d iblank_cell to the p4est quadrant storage */
    p4est_iterate(p4est,NULL,NULL,
                  iblank_ext_to_p4est_callback,
                  NULL,
            ARG3D(NULL)
                  NULL);

    /* visualization solution */
    if (ctx->visualize) {
        ctx->nvisualize++;
        vtk_write_solution(sim->visualization_path,
                           p4est,*index,
                           ctx->high_order_viz,
                           ctx->high_order_viz_npts);
    }
}

void driver_interface_point_inclusion(int *npts,Real *pts,int *cell_id){
    sc_array_t *pt_data;
    pt2_t *pt;
    int i,j;

    /* allocate the pt 2 structure */
    pt_data = sc_array_new_size(sizeof(pt2_t),*npts);

    /* initialize the point data */
    for (i = 0; i < *npts; ++i) {
        pt = (pt2_t *) sc_array_index(pt_data,i);

        pt->quad_id = -1;
        for (j = 0; j < DIM; ++j) {
            pt->x[j] = pts[DIM*i+j];
        }
    }

    /* point search using Newton's method */
    /* search if any of the points intersect the leafs */
    p4est_search(p4est,NULL,p4est_utilities_unst_pt_search_func,pt_data);

    /* move sc array flag to pointer */
    for (i = 0; i < *npts; ++i) {
        pt = (pt2_t *) sc_array_index(pt_data,i);
        cell_id[i] = pt->quad_id + BASE;
    }

    /* deallocate pt_data */
    sc_array_destroy(pt_data);
}


void driver_interface_bounding_box_intersection(int *nbb,
                                                Real *bb_xlo,Real *bb_xhi,
                                                int *bb_flag){
    sc_array_t *bb_data;
    bbox2_t *bb;
    int i,j;

    /* allocate the bounding box 2 structure */
    bb_data = sc_array_new_size(sizeof(bbox2_t),*nbb);

    /* initialize the bounding box data */
    for (i = 0; i < *nbb; ++i) {
        bb = (bbox2_t *) sc_array_index(bb_data,i);

        bb->box_hit = 0;
        for(j = 0; j < P4EST_DIM; ++j){
            bb->bbox_lo[j] = bb_xlo[P4EST_DIM*i+j];
            bb->bbox_hi[j] = bb_xhi[P4EST_DIM*i+j];
        }
    }

    /* search if bounding boxes intersect the leafs */
    /* NOTE: this search function works with linear element (NOT high-order) */
    p4est_search(p4est,NULL,p4est_utilities_unst_bbox2_search_func,bb_data);

    /* move sc array flag to pointer */
    for (i = 0; i < *nbb; ++i) {
        bb = (bbox2_t *) sc_array_index(bb_data,i);
        bb_flag[i] = bb->box_hit;
    }

    /* deallocate bb_data */
    sc_array_destroy(bb_data);
}

void driver_interface_get_data(int *body_tag,int *nnode,int *nwbc,int *nobc,
                               int *ntri,int *nquad,
                               int *ntet,int *npyramid,
                               int *nprism,int *nhex){
    ctx_t *ctx = &dg4est.d_ctx;
   *body_tag = 1;
   *nnode    = ctx->d_wake3d.nnodes;
   *nobc     = ctx->d_wake3d.nobc;
   *nwbc     = ctx->d_wake3d.nwbc;
   *ntri     = 0;
   *nquad    = (DIM==2) ? p4est->local_num_quadrants:0;
   *ntet     = 0;
   *npyramid = 0;
   *nprism   = 0;
   *nhex     = (DIM==3) ? p4est->local_num_quadrants:0;
}

void driver_interface_set_pointers(Real **ug_ptr,
                                   Real **xgeom_ptr,
                                   int **iblank_ptr,
                                   int **iwbcnode_ptr,
                                   int **iobcnode_ptr,
                                   int **ndcTri_ptr,
                                   int **ndcQuad_ptr,
                                   int **ndcTet_ptr,
                                   int **ndcPyr_ptr,
                                   int **ndcPrism_ptr,
                                   int **ndcHex_ptr,
                                   int **iblank_cell_ptr,
                                   void (**corn_ptr)(int*,int*),
                                   void (**crrn_ptr)(int*,int*,Real*),
                                   void (**dit_ptr) (int*,Real*,int*,Real*),
                                   void (**cdf_ptr) (int*,Real*,int*,int*,Real*,Real*,int*),
                                   void (**crc_ptr) (int*,int*,Real*,int*,int*,Real*)){
    ctx_t      *ctx = &dg4est.d_ctx;
    wake3d_t   *wke = &ctx->d_wake3d;

   *ug_ptr          = NULL;
   *xgeom_ptr       = wke->xgeom.ptr();
   *iblank_ptr      = wke->iblank.ptr();
   *iobcnode_ptr    = wke->iobcnode.ptr();
   *iwbcnode_ptr    = wke->iwbcnode.ptr();
   *ndcTri_ptr      = wke->ndcTri.ptr();
   *ndcQuad_ptr     = wke->ndcQuad.ptr();
   *ndcTet_ptr      = wke->ndcTet.ptr();
   *ndcPyr_ptr      = wke->ndcPyr.ptr();
   *ndcPrism_ptr    = wke->ndcPrism.ptr();
   *ndcHex_ptr      = wke->ndcHex.ptr();
   *iblank_cell_ptr = wke->iblank_cell.ptr();

   *corn_ptr = wake3d_count_receptor_nodes;
   *crrn_ptr = wake3d_create_receptor_nodes;
   *dit_ptr  = wake3d_donor_inclusion_test;
   *cdf_ptr  = wake3d_create_donor_frac;
   *crc_ptr  = wake3d_convert_to_receptor_coefficients;
}

/****************************************************/
/*               callback functions                 */
/****************************************************/

void wake3d_count_receptor_nodes(int *cell_index,int *num_receptor_nodes){
    p4est_quadrant_t *quad;
    p4est_locidx_t local_id;
    p4est_topidx_t which_tree;
    quad_data_t *data;

    local_id = *cell_index - BASE;
    p4est_utilities_get_quadrant_from_local_id(p4est,local_id,&quad,&which_tree);

    data = (quad_data_t *) quad->p.user_data;
    *num_receptor_nodes = 1;// FIXME
}

void wake3d_create_receptor_nodes(int *cell_index,int *num_receptor_nodes,
                                  Real *receptor_nodes){
    p4est_quadrant_t *quad;
    p4est_locidx_t local_id;
    p4est_topidx_t which_tree;
    quad_data_t *data;

    /* suppress compiler warning */
    (void) num_receptor_nodes;

    /* convert local id to c indexing */
    local_id = *cell_index - BASE;

    /* find the quadrant by searching over trees */
    p4est_utilities_get_quadrant_from_local_id(p4est,local_id,&quad,&which_tree);

    /* data held by quadrant */
    data = (quad_data_t *) quad->p.user_data;

    /* create receptor nodes */
    //FIXME
}

void wake3d_donor_inclusion_test(int *cell_index,Real *pt_xyz,
                                 int *pass_flag,Real *rst_out){
    p4est_quadrant_t *quad;
    p4est_locidx_t local_id;
    p4est_topidx_t which_tree;
    quad_data_t *data;

    /* convert local id to c indexing */
    local_id = *cell_index - BASE;

    /* find the quadrant by searching over trees */
    p4est_utilities_get_quadrant_from_local_id(p4est,local_id,&quad,&which_tree);

    /* data held by quadrant */
    data = (quad_data_t *) quad->p.user_data;

    /* donor inclusion test */
    //FIXME
}

void wake3d_create_donor_frac(int *cell_index,Real *xyz,int *nfrac,int *index,
                              Real *frac,Real *rst_in,int *dimfrac){
    p4est_quadrant_t *quad;
    p4est_locidx_t local_id;
    p4est_topidx_t which_tree;
    quad_data_t *data;

    /* suppress compiler warning */
    (void) xyz;
    (void) dimfrac;

    /* convert local id to c indexing */
    local_id = *cell_index - BASE;

    /* find the quadrant by searching over trees */
    p4est_utilities_get_quadrant_from_local_id(p4est,local_id,&quad,&which_tree);

    /* data held by quadrant */
    data = (quad_data_t *) quad->p.user_data;

    /* total number of modes in 2D/3D */
    *nfrac = 1;// FIXME

    /* solver create donor fraction weights */
    //FIXME

    /* index of solution array */
    *index = 0; //FIXME
}

void wake3d_convert_to_receptor_coefficients(int *cell_index,int *npts,Real *f,
                                             int *tm,int *index_out,Real *a_out){
    p4est_quadrant_t *quad;
    p4est_locidx_t local_id;
    p4est_topidx_t which_tree;
    quad_data_t  *data;

    /* convert local id to c indexing */
    local_id = *cell_index - BASE;

    /* find the quadrant by searching over trees */
    p4est_utilities_get_quadrant_from_local_id(p4est,local_id,&quad,&which_tree);

    /* data held by quadrant */
    data = (quad_data_t *) quad->p.user_data;

    /* convert receptor coefficients */
    //FIXME

    /* index of solution array to insert a_out */
    *index_out = 0;// FIXME
}

/* Internal functions for testing */
ctx_t * driver_interface_get_ctx(){
    return &dg4est.d_ctx;
}
p4est_t * driver_interface_get_p4est(){
    return p4est;
}
p4est_connectivity_t * driver_interface_get_p4est_conn(){
    return conn;
}
void driver_interface_set_logging(int log_threshold){
    ctx_t *ctx = &dg4est.d_ctx;
    ctx->log_info = log_threshold;
}