/**
 * \file    amr_initialize.cxx
 * \author  akirby
 *
 * \brief Initialization functions for the AMR code module.
 */

/* header files */
#include "amr_initialize.h"

void initialize_amr_file(char *filename,
                         dg4est_t *dg4est,
                         p4est_t **p4est,
                         p4est_connectivity_t **conn){
    ctx_t        *ctx = &dg4est->d_ctx;
    simulation_t *sim = &ctx->d_simulation;

    char regrid_now;
    char partition_now;
    char visualize_now;
    char noinput = 0;
    char reread_inputs = 0;

    sim->dim = P4EST_DIM;

    /* copy input file name */
    strcpy(sim->input_file,filename);

    /* build wrk directories */
    amr_utilities_create_directories(ctx->rank,ctx->comm,ctx->log_info);

    /* read the input file */
    initialize_inputs_from_file(sim->input_file,dg4est,noinput,reread_inputs);

    /* construct the grid and fill in data into quadrants */
    visualize_now = 0;
    partition_now = 1;
    regrid_now = 1;
    initialize_grid_and_solver(dg4est,p4est,conn,regrid_now,partition_now,visualize_now);
}

void initialize_inputs_from_file(char *filename,dg4est_t *dg4est,
                                 char noinput,char reread_inputs){
    ctx_t        *ctx = &dg4est->d_ctx;
    grid_t       *grd = &ctx->d_grid;
    gridfile_t *gfile = &ctx->d_gridfile;
    simulation_t *sim = &ctx->d_simulation;
    amr_hbox_t  *hbox = &ctx->d_amr_hboxes;
    multigrid_t   *mg = &dg4est->d_multigrid;

    char box_name_level[40];
    char box_name_auto[36];
    char box_name_lo[80];
    char box_name_hi[80];
    Real box_xlo[3];
    Real box_xhi[3];

    struct stat file_stat;
    int n;

    /* file state for checking runtime modifications */
    stat(filename,&file_stat);
    sim->input_file_mod_time = file_stat.st_mtime;

    /* default all inputs */
    initialize_default_inputs(ctx);

    /* read AMR inputs from file */
    if (!noinput) {
        amr_utilities_find_keyword_integer(filename,"log_display:",
                &ctx->log_info);

        if (!reread_inputs) {
            amr_utilities_find_keyword_integer(filename,"unstructured_flag:",
                    &sim->unstructured_flag);

            amr_utilities_find_keyword_string(filename,"unstructured_file:",
                    gfile->gridfile_name);

            amr_utilities_find_keyword_three_integers(filename,"periodic_flag:",
                    &grd->periodic[0],&grd->periodic[1],&grd->periodic[2]);

            amr_utilities_find_keyword_three_reals(filename,"domain_lo:",
                    &grd->xlo[0],&grd->xlo[1],&grd->xlo[2]);

            amr_utilities_find_keyword_three_reals(filename,"domain_hi:",
                    &grd->xhi[0],&grd->xhi[1],&grd->xhi[2]);

            amr_utilities_find_keyword_three_integers(filename,"nelem:",
                    &grd->nelem[0],&grd->nelem[1],&grd->nelem[2]);

            amr_utilities_find_keyword_integer(filename,"max_amr_level:",
                    &grd->max_level);

            amr_utilities_find_keyword_integer(filename,"min_amr_level:",
                    &grd->min_level);

            /* check if min/max are accidentally switched by the user. If so, flip them */
            if (grd->min_level > grd->max_level) {
                int max_level  = grd->min_level;
                grd->min_level = grd->max_level;
                grd->max_level = max_level;
            }

            grd->nlevels = grd->max_level + 1;
        }

        /* Read wake3d wall and outer patch lists */
        amr_utilities_find_keyword_string2upper(filename,"wake3d_wbc:",ctx->d_geometry.wake3d_wbc);
        amr_utilities_find_keyword_string2upper(filename,"wake3d_obc:",ctx->d_geometry.wake3d_obc);

        mg->p_levels = 1;

        /* I/O inputs */
        amr_utilities_find_keyword_integer(filename,"visualization_interval:",&sim->visualization_interval);

        /* check bounds on high order plotting */
        ctx->high_order_viz = 0;
        if(ctx->high_order_viz_npts < 1)  ctx->high_order_viz_npts = 2;
        if(ctx->high_order_viz_npts > 10) ctx->high_order_viz_npts = 10;

        /* regrid inputs */
        sim->regrid_strategy = 1;
        amr_utilities_find_keyword_integer(filename,"regrid_interval:",&sim->regrid_interval);

        /* amr_hboxes inputs */
        hbox->amr_nhboxes = 0;
        if (sim->regrid_interval > 0) {
            amr_utilities_find_keyword_integer(filename,"amr_nhboxes:",&hbox->amr_nhboxes);
            if (hbox->amr_nhboxes > 0) {
                sim->regrid_hboxes = 1;

                /* allocate memory */
                hbox->amr_hboxes_autotag.malloc(hbox->amr_nhboxes);
                hbox->amr_hboxes_level.malloc(hbox->amr_nhboxes);
                hbox->amr_hboxes_lo.malloc(3*hbox->amr_nhboxes);
                hbox->amr_hboxes_hi.malloc(3*hbox->amr_nhboxes);

                /* read boxes */
                for (n = 0; n < hbox->amr_nhboxes; ++n) {
                    sprintf(box_name_level,"amr_hbox%d_level:",n);
                    sprintf(box_name_lo,"amr_hbox%d_lo:",n);
                    sprintf(box_name_hi,"amr_hbox%d_hi:",n);
                    sprintf(box_name_auto,"amr_hbox%d_autotag:",n);

                    amr_utilities_find_keyword_integer(filename,box_name_auto,&hbox->amr_hboxes_autotag[n]);
                    amr_utilities_find_keyword_integer(filename,box_name_level,&hbox->amr_hboxes_level[n]);
                    amr_utilities_find_keyword_three_reals(filename,box_name_lo,
                            &box_xlo[0],&box_xlo[1],&box_xlo[2]);
                    amr_utilities_find_keyword_three_reals(filename,box_name_hi,
                            &box_xhi[0],&box_xhi[1],&box_xhi[2]);

                    hbox->amr_hboxes_lo[3*n+0] = box_xlo[0];
                    hbox->amr_hboxes_lo[3*n+1] = box_xlo[1];
                    hbox->amr_hboxes_lo[3*n+2] = box_xlo[2];
                    hbox->amr_hboxes_hi[3*n+0] = box_xhi[0];
                    hbox->amr_hboxes_hi[3*n+1] = box_xhi[1];
                    hbox->amr_hboxes_hi[3*n+2] = box_xhi[2];
                }
            }
        }

        /* reset sim->regrid_interval to max if no AMR */
        if(hbox->amr_nhboxes == 0 && (grd->max_level - grd->min_level == 0)) sim->regrid_interval = INT_MAX;

        if (sim->regrid_interval > 0) {
            amr_utilities_find_keyword_integer(filename,"regrid_features:",&sim->regrid_features);
            amr_utilities_find_keyword_integer(filename,"regrid_coarsen:",&sim->regrid_coarsen);
        } else {
            sim->regrid_features = 0;
            sim->regrid_coarsen = 0;
            sim->regrid_buffer = 0;
            sim->regrid_hboxes = 0;
        }
        if(sim->visualization_interval <= 0) sim->visualization_interval = INT_MAX;
        if(sim->regrid_interval <= 0) sim->regrid_interval = INT_MAX;

        /* set global max level */
        grd->max_level_global = grd->max_level;
        if (sim->regrid_hboxes) {
            for (int i = 0; i < hbox->amr_nhboxes; ++i) {
                grd->max_level_global = std::max(grd->max_level_global,hbox->amr_hboxes_level[i]);
            }
        }

        amr_utilities_inputs_amr_message(ctx,ctx->rank);
    }
}

void initialize_default_inputs(ctx_t *ctx){
    gridfile_t *gfile = &ctx->d_gridfile;
    simulation_t *sim = &ctx->d_simulation;

    strcpy(gfile->gridfile_name,"default.msh");
    strcpy(sim->visualization_path,"WRK/");
}

void initialize_grid_and_solver(dg4est_t *dg4est,
                                p4est_t **p4est,
                                p4est_connectivity_t **conn,
                                char regrid_now,
                                char partition_now,
                                char visualize_now){
    ctx_t        *ctx = &dg4est->d_ctx;
    simulation_t *sim = &ctx->d_simulation;

    int initial = 1;

    /* set initial flag */
    ctx->initial = initial = 1;

    /* build new grid */
    initialize_grid_constants(ctx);

    /* read gridfile, initialize connectivity, and reorder */
    initialize_grid(ctx,conn);

    /* construct level 0 grid: create forest and read element coordinates */
    *p4est = initialize_p4est(dg4est,ctx,*conn);

    if (ctx->rank == 0 && ctx->log_info < DGLOG_ALL) {
        DGOUT(ctx->log_io,"[ amr ] Done initializing quadrant data\n");
    }

    /* deallocate gridfile element line bookkeeper */
    gridfile_deallocate_data(ctx);

    /* assign boundary conditions */
    p4est_iterate(*p4est,NULL,(void *) ctx,
                         NULL,
                         initialize_quadrant_bc_data,
                   ARG3D(NULL)
                         NULL);

    if (ctx->rank == 0 && ctx->log_info < DGLOG_ALL) {
        DGOUT(ctx->log_io,"[ amr ] Done initializing quadrant bc data\n");
    }

    /* deallocate quadrant faces info */
    if(sim->unstructured_flag) geometry_deallocate_face_data(ctx);

    /* regrid to min_level */
    regrid_2minlevel(ctx,*p4est,0);

    /* calculate domain volume */
    initialize_grid_volume(ctx,*p4est);

    /* display grid information */
    amr_utilities_grid_message(ctx);

    /* adapt grid initial solution */
    if(regrid_now) regrid_solution(ctx,*p4est,initial);

    /* partition grid */
    if(partition_now) regrid_partition(*p4est);

    /* visualize initial solution*/
    if(visualize_now) vtk_write_solution(sim->visualization_path,*p4est,1,
                                         ctx->high_order_viz,ctx->high_order_viz_npts);

    /* reset initial flag after building grid */
    ctx->initial = 0;
}

void initialize_grid_constants(ctx_t *ctx){
    grid_t *grd = &ctx->d_grid;

    Real mesh_scale[3];
    Real base_length;
    Real x_length;
    Real y_length;
    Real z_length;

    /* grid data is taken from the Cartesian grid in the input file. */
    /* determine the coarse level size */
    x_length = grd->xhi[0] - grd->xlo[0];
    y_length = grd->xhi[1] - grd->xlo[1];
    z_length = grd->xhi[2] - grd->xlo[2];

    mesh_scale[0] = x_length / (Real) grd->nelem[0];
    mesh_scale[1] = y_length / (Real) grd->nelem[1];
    mesh_scale[2] = z_length / (Real) grd->nelem[2];

    grd->max_dx[0] = mesh_scale[0];
    grd->max_dx[1] = mesh_scale[1];
    grd->max_dx[2] = mesh_scale[2];

    base_length = (Real) (1 << grd->max_level);
    grd->min_dx[0] = grd->max_dx[0] / base_length;
    grd->min_dx[1] = grd->max_dx[1] / base_length;
    grd->min_dx[2] = grd->max_dx[2] / base_length;
}

void initialize_grid(ctx_t *ctx,p4est_connectivity_t **conn){
    grid_t       *grd   = &ctx->d_grid;
    simulation_t *sim   = &ctx->d_simulation;
    gridfile_t   *gfile = &ctx->d_gridfile;

    int nelem_cLevel_x;
    int nelem_cLevel_y;
    int nelem_cLevel_z;
    int i;

    nelem_cLevel_x = grd->nelem[0];
    nelem_cLevel_y = grd->nelem[1];
    nelem_cLevel_z = grd->nelem[2];

    /* create a new forest */
    if (!sim->unstructured_flag) {
        /* ================ */
        /* Cartesian Domain */
        /* ================ */
#ifndef P4_TO_P8
        *conn = p4est_connectivity_new_brick(nelem_cLevel_x,
                                             nelem_cLevel_y,
                                             grd->periodic[0],
                                             grd->periodic[1]);
#else
        *conn = p8est_connectivity_new_brick(nelem_cLevel_x,
                                             nelem_cLevel_y,
                                             nelem_cLevel_z,
                                             grd->periodic[0],
                                             grd->periodic[1],
                                             grd->periodic[2]);
#endif

        /* scale mesh and translate */
        Real *vertices = (*conn)->vertices;
        for (i = 0; i < (*conn)->num_vertices; ++i) {
            vertices[3*i+0] = vertices[3*i+0]*grd->max_dx[0] + grd->xlo[0];
            vertices[3*i+1] = vertices[3*i+1]*grd->max_dx[1] + grd->xlo[1];
            vertices[3*i+2] = vertices[3*i+2]*grd->max_dx[2] + grd->xlo[2];
        }
    } else {
        /* =================== */
        /* Unstructured Domain */
        /* =================== */
        sc_array_t *newid;
        int *ip;
        int i;

        /* read unstructured grid */
        *conn = gridfile_reader(ctx,gfile->gridfile_name);

        /* Reorder connectivity */
#ifdef P4EST_WITH_METIS
        newid = sc_array_new_count(sizeof(size_t), 10);
        p4est_connectivity_reorder_newid(ctx->comm,0,*conn,P4EST_CONNECT_FACE,newid);

        /* copy newid */
        for (i = 0; i < newid->elem_count; i++) {
            ip = (int *) sc_array_index(newid, i);
            gfile->p4est_part[ip[0]] = i;
        }

        /* destroy newid */
        sc_array_destroy(newid);
#endif /* P4EST_WITH_METIS */
    }
}

p4est_t *initialize_p4est(dg4est_t *dg4est,ctx_t *ctx,p4est_connectivity_t *conn){
    return p4est_new_ext(ctx->comm,                /* mpi communicator */
                         conn,                     /* p4est connectivity */
                         0,                        /* minimum quadrants per MPI process */
                         0,                        /* minimum level of mesh refinement */
                         1,                        /* fill the forest with a uniform mesh */
                         sizeof(quad_data_t),      /* user quadrant data type */
                         initialize_quadrant_data, /* p4est callback for data initialization */
                         dg4est);                  /* user data pointer to be stored in p4est */
}

void initialize_grid_volume(ctx_t *ctx,p4est_t *p4est){
    grid_t *grd = &ctx->d_grid;
    Real local_volume = 0.0;

    /* calculate local grid volume */
    p4est_iterate(p4est,NULL,(void *)&local_volume,
                  domain_volume_callback,
                  NULL,
            ARG3D(NULL)
                  NULL);

    /* accumulate total grid volume */
    MPI_Allreduce(&local_volume,&grd->domain_volume,1,MPI_DGREAL,MPI_SUM,ctx->comm);
}

void initialize_quadrant_data(p4est_t *p4est,
                              p4est_topidx_t which_tree,
                              p4est_quadrant_t *quad){
    Real dx[3] = {0.0};
    int display_quad = 0;
    int analytic_solution = 0;
    int qorder;
    int ibc;
    Uint vol_tag;

    quad_properties_t quad_prop;
    quad_user_data_t user_data;

    /* access our simulation context data */
    dg4est_t  *dg4est = (dg4est_t *) p4est->user_pointer;
    ctx_t        *ctx = &dg4est->d_ctx;
    geometry_t  *geom = &ctx->d_geometry;
    grid_t       *grd = &ctx->d_grid;
    gridfile_t *gfile = &ctx->d_gridfile;
    simulation_t *sim = &ctx->d_simulation;

    /* access our per quadrant data pointer */
    quad_data_t *data = (quad_data_t *) quad->p.user_data;

    /* calculate the length ratio between this level and level 0 */
    const Real base_length =
        (Real) P4EST_QUADRANT_LEN (quad->level) /
        (Real) P4EST_ROOT_LEN;

    /* quadrant length in each direction: dx */
    dx[0] = grd->max_dx[0] * base_length;
    dx[1] = grd->max_dx[1] * base_length;
    dx[2] = grd->max_dx[2] * base_length;

    /* set linear elements (even if high order for now) */
    qorder = 1;

    /* get quad geometry: xyz */
    Real xyz[3*P4EST_CHILDREN];
    p4est_utilities_quad_coordinates(p4est,which_tree,quad,xyz,display_quad);

    /* assign quad properties */
    quad_prop.dim          = sim->dim;
    quad_prop.level        = quad->level;
    quad_prop.qdegree      = qorder;
    quad_prop.quad_dx      = dx;
    quad_prop.unstructured = sim->unstructured_flag;

    /* set data defaults */
    data->tag       = no_tag;
    data->type      = NO_PATCH;
    data->qdegree   = qorder;
    data->iblank    = field_mask;

    data->soln_size = 0;
    data->soln      = NULL;
    data->geom_size = 3*DGVDIM(qorder+1)*sizeof(Real);
    data->geom      = (Real *) malloc(data->geom_size);
    memcpy(data->geom,xyz,data->geom_size);

    /* assign level */
    data->level = quad->level;

    /* set default bcs */
    for(ibc = 0; ibc < 6; ibc++) data->bc[ibc] = NO_PATCH;

    /* assign data pointers from defaults above */
    user_data.qdegree   = data->qdegree;
    user_data.geom_size = 0;
    user_data.quad_geom = NULL;
}

int find_geom_face(geometry_t *geom,p4est_topidx_t *vquad,int iface){
    int nnodes = P4EST_CHILDREN/2;
    int nface_info = nnodes + 2;
    p4est_topidx_t vface[nnodes];
    int i,j,k;

    /* get face node list */
    if (iface == 0) { /* x_lo */
        vface[0] = vquad[0];
        vface[1] = vquad[2];
    X3D(vface[2] = vquad[4];)
    X3D(vface[3] = vquad[6];)
    } else
    if (iface == 1) { /* x_hi */
        vface[0] = vquad[1];
        vface[1] = vquad[3];
    X3D(vface[2] = vquad[5];)
    X3D(vface[3] = vquad[7];)
    } else
    if (iface == 2) { /* y_lo */
        vface[0] = vquad[0];
        vface[1] = vquad[1];
    X3D(vface[2] = vquad[4];)
    X3D(vface[3] = vquad[5];)
    } else
    if (iface == 3) { /* y_hi */
        vface[0] = vquad[2];
        vface[1] = vquad[3];
    X3D(vface[2] = vquad[6];)
    X3D(vface[3] = vquad[7];)
    } else
    if (iface == 4) { /* z_lo */
        vface[0] = vquad[0];
        vface[1] = vquad[1];
        vface[2] = vquad[2];
        vface[3] = vquad[3];
    } else
    if (iface == 5) { /* z_hi */
        vface[0] = vquad[4];
        vface[1] = vquad[5];
        vface[2] = vquad[6];
        vface[3] = vquad[7];
    }

    /* loop through each face in mesh */
    for (i = 0; i < geom->nelem_face; i++) {
        const Uint *this_face = &geom->face_nodes_info[nface_info*i + EFNDS_IND];

        /* loop through each vertex of the current (vface) */
        char nnode_match = 1;
        for (j = 0; j < nnodes; j++) {
            char inode_match = 0;

            /* check if each current face vertex is part of the i-face  */
            for (k = 0; k < nnodes; k++) {
                inode_match |= (vface[j] == this_face[k]);
                if(inode_match) break; // found node -- break from search
            }
            nnode_match &= inode_match;
            if(!inode_match) break;    // inode_match=0 if no nodes were found so stop
        }
        if(nnode_match) return i;
    }

    /* no face match */
    return -1;
}

void initialize_quadrant_bc_data(p4est_iter_face_info_t *info,void *user_data){
    p4est_t                       *p4est = info->p4est;
    p4est_connectivity_t   *connectivity = p4est->connectivity;
    const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
    sc_array_t *sides = &(info->sides);

    ctx_t *ctx = (ctx_t *) user_data;
    geometry_t *geom = &ctx->d_geometry;

    const int nface_info = 2 + P4EST_CHILDREN/2;
    p4est_topidx_t vquad[P4EST_CHILDREN];

    p4est_iter_face_side_t *side[2];
    p4est_topidx_t which_tree;
    int fullside,hanging;

    Uint face_tag,entity_tag;
    int patch_id, ipatch;
    int face_id;
    int iface;
    int i,j,k;

    /* Check if the face is a forest boundary */
    int boundary = (sides->elem_count > 1) ? 0:1;

    if (boundary) {
        /* access our per quadrant data pointer */
        side[0] = p4est_iter_fside_array_index_int(sides,0);
        quad_data_t *data = (quad_data_t *) side[0]->is.full.quad->p.user_data;

        /* get quad side */
        iface = side[0]->face;

        if (ctx->d_simulation.unstructured_flag) {
            /* get quadrant node list (only corner vertices) */
            which_tree = side[0]->treeid;

            for (k = 0; k < P4EST_CHILDREN; ++k) {
                vquad[k] = tree_to_vertex[which_tree*P4EST_CHILDREN + k];
            }

            /* find face ID corresponding to this quadrants iface */
            face_id = find_geom_face(geom,vquad,iface);
            if(face_id >= 0){
                entity_tag = geom->face_nodes_info[face_id*(nface_info) + EFENT_IND];
                face_tag   = geom->face_nodes_info[face_id*(nface_info) + EFTAG_IND];

                /* get entity tag index */
                for (ipatch = 0; ipatch < geom->nentity_face; ipatch++) {
                    if (geom->entity_face_tag[ipatch] == entity_tag) {
                        patch_id = ipatch;
                        break;
                    }
                }

                /* write the physical boundary tag of this face */
                data->bc[iface] = geom->entity_face_patchtag[patch_id];

                /* Sanity check:
                 *  check that patch_tag is a number between 0 and the number of patches.
                 *  This might happen if the grid is not generated properly.
                 */
                if (data->bc[iface] < 0 || data->bc[iface] > (geom->npatch-1)) {
                    DGOUT(stderr,"\x1B[1;31m[ amr ] Boundary face is not associated to any physical patch. Check grid file!\x1B[0m\n");
                    exit(0);
                }
                return;
            }
        } else {
            /* Structured Grids ONLY:
             *   For WAKE3D/tioga bc registration.
             *   (see amr_external_func.c:setup_wake3d_allnodes)
             */
            data->bc[iface] = 1; // not equal to -1
            return;
        }
        /* if nothing is found, raise an error */
        DGOUT(stderr,"\x1B[1;31m[ amr ] Boundary face not found. No boundary condition assigned!\x1B[0m\n");
        exit(0);
    }
}

void domain_volume_callback(p4est_iter_volume_info_t *info,void *user_data){
    Real *local_volume = (Real *) user_data;
    p4est_quadrant_t *quad = info->quad;
    quad_data_t *data = (quad_data_t *) quad->p.user_data;

    Real quad_volume = 0.0;

    /* calculate volume of octant */
    // FIXME

    /* add to local volume */
    *local_volume += quad_volume;
}