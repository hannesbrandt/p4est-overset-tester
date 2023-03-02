/**
 * \file    amr_vtk.cxx
 * \author  akirby
 *
 * \brief   Visualization functions for the AMR code module.
 */

/* header files */
#include "amr_vtk.h"
#include "amr_var_data.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

void vtk_visualize_grid(p4est_t *p4est){
    p4est_vtk_write_file(p4est,NULL,P4EST_STRING "_grid");
}

Real vtk_write_solution(char *filepath,p4est_t *p4est,int timestep,
                        int HO_VIZ,int HO_NPTS){
    dg4est_t  *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t        *ctx = &dg4est->d_ctx;
    simulation_t *sim = &ctx->d_simulation;
    Real t1,t2;

    t1 = amr_utilities_timer(); /* start timer */

    char solution_name[BUFF_SIZE] = {'\0'};
    (HO_VIZ) ? snprintf(solution_name,19,"HOVTU.soln_%07d",timestep):
               snprintf(solution_name,13,"soln_%07d",timestep);

    char *filename = (char *) malloc(strlen(filepath) + strlen(solution_name) + 1);
    strcpy(filename,filepath);
    strcat(filename,solution_name);

    p4est_locidx_t numquads = p4est->local_num_quadrants;
    vtk_plot_fluid_t fluid_data;

    const int NPLOTPTS = (HO_VIZ) ? pow(HO_NPTS,DIM):P4EST_CHILDREN;
    fluid_data.high_order_viz = HO_VIZ;
    fluid_data.nplotpts1D = (HO_VIZ) ? HO_NPTS:2;
    fluid_data.nplotpts = NPLOTPTS;

    if(HO_VIZ) fluid_data.positions = sc_array_new_size(sizeof(Real),numquads*NPLOTPTS*DIM);

    fluid_data.tag_id  = sc_array_new_size(sizeof(Real),numquads);
    fluid_data.type_id = sc_array_new_size(sizeof(Real),numquads);
    fluid_data.quad_id = sc_array_new_size(sizeof(Real),numquads);
    fluid_data.iblank  = sc_array_new_size(sizeof(Real),numquads);

    p4est_iterate(p4est,
                  NULL,
        (void *) &fluid_data,
                  interpolate_solution,
                  NULL,
            ARG3D(NULL)
                  NULL);

    /* create VTK output context and set its parameters */
    p4est_vtk_context_t *context = p4est_vtk_context_new(p4est,filename);
    p4est_vtk_context_set_scale(context,1.0);
    p4est_vtk_context_set_continuous(context,1);

    /* begin writing the output files */
    context = (HO_VIZ) ? p4est_vtk_write_header_ho(context,fluid_data.positions,HO_NPTS):
                         p4est_vtk_write_header(context);

    SC_CHECK_ABORT(context != NULL,P4EST_STRING "_vtk: Error writing vtk header");

    context = p4est_vtk_write_cell_dataf(context,
                                         1,  /* write tree */
                                         1,  /* write level */
                                         1,  /* write mpi rank */
                                         0,  /* do not wrap the mpi rank */
                                         4,  /* num of quad scalar fields */
                                         0,  /* num of quad vector fields. */
                                         "AMR Tag",fluid_data.tag_id,
             ((sim->unstructured_flag) ? "Physical ID":"BC"),fluid_data.type_id,
                                         "Quad ID",fluid_data.quad_id,
                                         "Iblank",fluid_data.iblank,
                                         context); /* mark the end of the variable quad data. */

    SC_CHECK_ABORT(context != NULL,P4EST_STRING "_vtk: Error writing quad data");

    const int retval = p4est_vtk_write_footer (context);
    SC_CHECK_ABORT(!retval, P4EST_STRING "_vtk: Error writing footer");


    if(HO_VIZ) sc_array_destroy(fluid_data.positions);
    if(filename) free(filename); filename = NULL;

    t2 = amr_utilities_timer(); /* stop timer */
    return amr_utilities_mpireducemax_real(p4est->mpirank,p4est->mpicomm,t2-t1);
}

void interpolate_solution(p4est_iter_volume_info_t *info,void *user_data){
    vtk_plot_fluid_t *soln_viz = (vtk_plot_fluid_t *) user_data;
    p4est_t *p4est = info->p4est;
    p4est_quadrant_t *q = info->quad;
    p4est_topidx_t which_tree = info->treeid;
    p4est_locidx_t local_id = info->quadid;

    p4est_tree_t  *tree;
    p4est_locidx_t arrayoffset;

    quad_data_t *data = (quad_data_t *) q->p.user_data;
    dg4est_t  *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t        *ctx = &dg4est->d_ctx;
    simulation_t *sim = &ctx->d_simulation;

    Real *this_quad_int_ptr;

    int nqp;
    int i,j,k;

    const int NPLOTPTS1D = soln_viz->nplotpts1D;
    const int NPLOTPTS = soln_viz->nplotpts;

    Real GLL[NPLOTPTS1D];
    Real plot_pts[P4EST_DIM*NPLOTPTS];
    Real xyz[data->geom_size];

    /* calculate quadrant geometry coordinates */
    // p4est_utilities_quad_coordinates(p4est,which_tree,q,xyz,0);
    memcpy(xyz,data->geom,data->geom_size);

    tree = p4est_tree_array_index(p4est->trees, which_tree);

    /* now the id is relative to the MPI process */
    local_id += tree->quadrants_offset;

    /* each local quadrant has NPLOTPTS values in soln_viz */
    arrayoffset = NPLOTPTS * local_id;

    /* fill in plot points to project solution */
    nqp = soln_viz->nplotpts1D;

    /* fill GLL points */
    GLL_points(NPLOTPTS1D,GLL);

    /* fill plot points */
#ifdef P4_TO_P8
    int npt = 0;
    for (k = 0; k < nqp; ++k) {
        for (j = 0; j < nqp; ++j) {
            for (i = 0; i < nqp; ++i) {
                Real *pt = &(plot_pts[P4EST_DIM*npt]);
                pt[0] = GLL[i]; // x
                pt[1] = GLL[j]; // y
                pt[2] = GLL[k]; // z
                npt++;
            }
        }
    }
#else
    int npt = 0;
    for (j = 0; j < nqp; ++j) {
        for (i = 0; i < nqp; ++i) {
            Real *pt = &(plot_pts[P4EST_DIM*npt]);
            pt[0] = GLL[i]; // x
            pt[1] = GLL[j]; // y
            npt++;
        }
    }
#endif

    /* fill in cell values */
    this_quad_int_ptr = (Real *) sc_array_index(soln_viz->quad_id,local_id);
    this_quad_int_ptr[0] = local_id;

    this_quad_int_ptr = (Real *) sc_array_index(soln_viz->tag_id,local_id);
    this_quad_int_ptr[0] = data->tag;

    int bc = 0;
    if(sim->unstructured_flag){
        // data->bc[i] contains physical boundary tag
        for(i = 0; i < 2*DIM; ++i){
            if(data->bc[i] != -1) {bc = data->bc[i]; break;}
        }
    } else {
        /* structured mesh */
#define BC_CHARACTERISTIC 4
        for(i = 0; i < 2*DIM; ++i){
            if(data->bc[i] != BC_CHARACTERISTIC) bc = MAX(bc,data->bc[i]);
        }
        if(bc == 0) bc = BC_CHARACTERISTIC;
    }
    this_quad_int_ptr = (Real *) sc_array_index(soln_viz->type_id,local_id);
    this_quad_int_ptr[0] = bc;

    this_quad_int_ptr = (Real *) sc_array_index(soln_viz->iblank,local_id);
    this_quad_int_ptr[0] = data->iblank;
}