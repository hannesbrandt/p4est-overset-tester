/**
 * \file    amr_geometry.cxx
 * \author  mugolott
 *
 * \brief Geometry
 */

/* header files */
#include "amr_geometry.h"

void geometry_allocate_patch_data(geometry_t *geom){
    /* allocate memory for patch names */
    geom->patch_dim.malloc(geom->npatch,NO_PATCH);
    geom->wake3d_bc.malloc(geom->npatch,W3D_NO);

    /* initialize arrays */
    geom->patch_name = (char **) malloc(geom->npatch*sizeof(char *));
    for(int ii = 0; ii < geom->npatch; ii++) geom->patch_name[ii] = (char *) malloc(BUFF_SIZE);
}

void geometry_allocate_entity_data(geometry_t *geom){
    /* allocate memory for patch info */
    geom->entity_face_tag.malloc(geom->nentity_face,NO_TAG);
    geom->entity_face_patchtag.malloc(geom->nentity_face,NO_TAG);
    geom->entity_vol_tag.malloc(geom->nentity_vol,NO_TAG);
    geom->entity_vol_patchtag.malloc(geom->nentity_vol,NO_TAG);
}

void geometry_allocate_facenode_data(geometry_t *geom,size_t num_faces){
    int nface_info = 2*DIM;
    geom->face_nodes_info.malloc(nface_info*num_faces);
}

void geometry_deallocate_data(ctx_t *ctx){
    geometry_t *geom = &ctx->d_geometry;
    int i;

    /* deallocate memory */
    if (ctx->d_simulation.unstructured_flag) {
        for (i = 0; i < geom->npatch; i++) {
            free(geom->patch_name[i]);
                 geom->patch_name[i] = NULL;
        }
        if(geom->patch_name) {free(geom->patch_name);} geom->patch_name = NULL;

        /* set counters to zero */
        geom->npatch = 0;
        geom->nentity_face = 0;
        geom->nentity_vol = 0;
        geom->nelem_face = 0;
    }
}

void geometry_deallocate_face_data(ctx_t *ctx){
    geometry_t *geom = &ctx->d_geometry;

    geom->entity_vol_tag.free();
    geom->face_nodes_info.free();
    geom->entity_face_tag.free();
    geom->entity_vol_patchtag.free();
    geom->entity_face_patchtag.free();
}

void geometry_broadcast_data(ctx_t *ctx){
    geometry_t *geom = &ctx->d_geometry;
    int ii;

    /* send single quantities first: needed for allocation */
    MPI_Bcast(&geom->npatch,      1,MPI_INT,0,ctx->comm);
    MPI_Bcast(&geom->nentity_face,1,MPI_UNSIGNED_LONG_LONG,0,ctx->comm);
    MPI_Bcast(&geom->nentity_vol, 1,MPI_UNSIGNED_LONG_LONG,0,ctx->comm);
    MPI_Bcast(&geom->nelem_face,  1,MPI_UNSIGNED_LONG_LONG,0,ctx->comm);

    /* allocate memory on other ranks */
    if(ctx->rank != 0) geometry_allocate_patch_data(geom);
    if(ctx->rank != 0) geometry_allocate_entity_data(geom);
    if(ctx->rank != 0) geometry_allocate_facenode_data(geom,geom->nelem_face);

    /* broadcast arrays */
    MPI_Bcast(geom->patch_dim.ptr(),           geom->npatch,MPI_INT,0,ctx->comm);
    MPI_Bcast(geom->entity_face_tag.ptr(),     geom->nentity_face,MPI_UNSIGNED_LONG_LONG,0,ctx->comm);
    MPI_Bcast(geom->entity_face_patchtag.ptr(),geom->nentity_face,MPI_INT,0,ctx->comm);
    MPI_Bcast(geom->entity_vol_tag.ptr(),      geom->nentity_vol,MPI_UNSIGNED_LONG_LONG,0,ctx->comm);
    MPI_Bcast(geom->entity_vol_patchtag.ptr(), geom->nentity_vol,MPI_INT,0,ctx->comm);
    MPI_Bcast(geom->face_nodes_info.ptr(),     2*DIM*geom->nelem_face,MPI_UNSIGNED_LONG_LONG,0,ctx->comm);
    MPI_Bcast(geom->wake3d_bc.ptr(),           geom->npatch,MPI_INT,0,ctx->comm);
    MPI_Bcast(geom->wake3d_wbc,BUFF_SIZE,MPI_CHAR,0,ctx->comm);
    MPI_Bcast(geom->wake3d_obc,BUFF_SIZE,MPI_CHAR,0,ctx->comm);

    for(ii = 0; ii < geom->npatch; ii++) MPI_Bcast(geom->patch_name[ii],BUFF_SIZE,MPI_CHAR,0,ctx->comm);
}