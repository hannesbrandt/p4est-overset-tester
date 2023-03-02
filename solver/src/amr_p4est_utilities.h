/**
 * \file   amr_p4est_utilities.h
 * \author akirby
 */

#ifndef AMR_P4EST_UTILITIES_H
#define AMR_P4EST_UTILITIES_H

/* header files */
#include "amr_var_tag.h"
#include "amr_var_quad.h"
#include "amr_containers.h"
#include "dg4est_solver.hxx"
#include "precision_types.h"

#ifndef DOXYGEN_IGNORE
#  include <mpi.h>
#  include <math.h>
#  include <inttypes.h>
#endif

#ifndef P4_TO_P8
#  include <p4est_extended.h>
#  include <p4est_search.h>
#  include <p4est_nodes.h>
#else
#  include <p8est_extended.h>
#  include <p8est_search.h>
#  include <p8est_nodes.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    Real x[2];
}
vec2;

typedef struct {
    Real x[3];
}
vec3;

typedef struct {
    vec3 min;
    vec3 max;
}
box;

typedef struct {
    Real x[P4EST_DIM];  /**< pt coordinates */
    char *quad_hit;     /**< flag list of intersected quadrants */
}
pt_t;

typedef struct {
    Real x[P4EST_DIM];  /**< pt coordinates */
    int quad_id;        /**< intersected quadrant id */
}
pt2_t;

typedef struct {
    Real ray_origin[P4EST_DIM];   /**< ray origin coordinates */
    Real ray_end[P4EST_DIM];      /**< ray end coordinates */
    char *quad_hit;               /**< flag list of intersected quadrants */
}
ray_t;

typedef struct {
    Real bbox_lo[P4EST_DIM];  /**< bounding box lower coordinates */
    Real bbox_hi[P4EST_DIM];  /**< bounding box upper coordinates */
    char *quad_hit;           /**< flag list of intersected quadrants */
}
bbox_t;

typedef struct {
    Real bbox_lo[P4EST_DIM];  /**< bounding box lower coordinates */
    Real bbox_hi[P4EST_DIM];  /**< bounding box upper coordinates */
    char box_hit;             /**< flag list of intersected quadrants */
}
bbox2_t;

typedef struct {
    Real lo;  /**< lower bound */
    Real hi;  /**< upper bound */
}
bound_t;

typedef struct {
    bound_t x;  /**< x bounds */
    bound_t y;  /**< y bounds */
    bound_t z;  /**< z bounds */
}
box_t;

typedef struct {
    char *ghost_data;
    char correct_flag;
    char display_flag;
}
iblank_data_t;

void iblank_ext_to_p4est_callback(p4est_iter_volume_info_t *info,void *user_data);

/** Displays the mesh statistics to the terminal
 *
 * @param [in] p4est        p4est tree structure
 */
void p4est_utilities_mesh_stats(p4est_t *p4est);

/** Calculates amr quadrant geometry coordinates
 *    * NOTE: This is only suitable for linear elements! This assumption implies
 *    *       that the Jacobian is constant for all reference points
 *    *       located inside the element. Thus we take the center point of
 *    *       the quadrant for the Jacobian transformation.
 *
 * @param [in]  p4est           p4est tree structure
 * @param [in]  which_tree      p4est tree in forest of trees
 * @param [in]  quad            p4est quadrant
 * @param [out] xyz             quad geometry coordinates
 * @param [in]  display_quad    flag to display the quad geometry coordinates
 */
void p4est_utilities_quad_coordinates(p4est_t *p4est,p4est_topidx_t which_tree,
                                      p4est_quadrant_t *quad,Real *xyz,
                                      int display_quad);

/** Calculates amr high-order quadrant geometry coordinates
 *    * NOTE: This is only used to initialize the parent quads coordiantes
 *    *       at the get-go when we read in the grid file.
 *
 * @param [in]  p4est           p4est tree structure
 * @param [out] xyz             quad geometry coordinates
 * @param [in]  quad_nodes      Data type with array containing quad's node IDs (z-ordered)
 * @param [in]  display_quad    flag to display the quad geometry coordinates
 */
void p4est_utilities_init_quad_coordinates(p4est_t *p4est,Real *xyz,
                                           ivector_t *quad_nodes,
                                           int display_quad);

/** Calculates amr quadrant surface area from geometry.
 *
 * @param [in] xyz_quad     quad geometry coordinates
 * @return surface area of quadrant
 */
Real p4est_utilities_surface_area_calc(Real *xyz_quad);

/** Calculates amr quadrant volume from geometry.
 *
 * @param [in] xyz_quad     quad geometry coordinates
 * @return volume of quadrant
 */
Real p4est_utilities_volume_calc(Real *xyz_quad);

/** Calculates amr quadrant volume
 *
 * @param [in] p4est        p4est tree structure
 * @param [in] which_tree   p4est tree of this quadrant
 * @param [in] quad         p4est quadrant
 * @return     volume of quadrant
 */
Real p4est_utilities_quad_volume(p4est_t *p4est,p4est_topidx_t which_tree,
                                   p4est_quadrant_t *quad);

/** Calculates AMR unique grid node geometry coordinates
 *
 * @param [inout] xgeom     node geometry to fill out
 * @param [in]    p4est     p4est tree structure
 * @param [in]    indeps    p4est independent nodes
 */
void p4est_utilities_indep_node_coordinates(Real *xgeom,p4est_t *p4est,
                                            sc_array_t *indeps);

/** Calculates amr quadrant/octant lower corner geometry coordinates
 *
 * @param [in]  p4est           p4est tree structure
 * @param [in]  which_tree      p4est tree in forest of trees
 * @param [in]  q               p4est quadrant
 * @param [out] xyz             quad lower geometry coordinates
 */
void p4est_utilities_get_xlo(p4est_t *p4est,p4est_topidx_t which_tree,
                             p4est_quadrant_t *q,Real xyz[3]);

/** Provides local id of quad volume
 *
 * @param [in] info     p4est quad volume information
 * @return local id of quadrant
 */
p4est_locidx_t p4est_utilities_get_local_id_volume(p4est_iter_volume_info_t *info);

/** Provides local id of quad from hanging face
 *
 * @param [in] info         p4est face information
 * @param [in] which_side   p4est which side face is on
 * @param [in] which_hang   p4est which face is hanging
 * @return local id of quadrant
 */
p4est_locidx_t p4est_utilities_get_local_id_face_hang(p4est_iter_face_info_t *info,
                                                      int which_side, int which_hang);

/** Provides local id of quad from full face
 *
 * @param [in] info         p4est face information
 * @param [in] which_side   p4est which side face is on
 * @return local id of quadrant
 */
p4est_locidx_t p4est_utilities_get_local_id_face_full(p4est_iter_face_info_t *info,
                                                      int which_side);

/** Provides the p4est quadrant and tree given a local id
 *
 * @param p4est         p4est tree structure
 * @param local_id      local quadrant id
 * @param quad          p4est quadrant corresponding to the \a local_id
 * @param which_tree    p4est tree containing \a local_id
 */
void p4est_utilities_get_quadrant_from_local_id(p4est_t *p4est,
                                                p4est_locidx_t local_id,
                                                p4est_quadrant_t **quad,
                                                p4est_topidx_t *which_tree);

/** Search utility function for searching for points in 2D unstructured mesh using Newton's method
 *
 * @param [in]  qdegree quadrant geometry order of accuracy
 * @param [in]  geom    quadrant geometry coordinates (z-ordered)
 * @param [in]  xyz     physical point coordinates to search
 * @param [out] rst     natural point coordinates
 * @return  flag indicator for point found in quadrant
 */
int p4est_pt_rst_coordinates_2d(int qdegree, Real *geom,Real *xyz,Real *rst);

/** Search utility function for searching for points in 3D unstructured mesh using Newton's method
 *
 * @param [in]  qdegree quadrant geometry order of accuracy
 * @param [in]  geom    quadrant geometry coordinates (z-ordered)
 * @param [in]  xyz     physical point coordinates to search
 * @param [out] rst     natural point coordinates
 * @return  flag indicator for point found in quadrant
 */
int p4est_pt_rst_coordinates_3d(int qdegree, Real *geom,Real *xyz,Real *rst);

/** Search utility function for checking boundary-box-surface/quadrant intersection.
 *
 * @param [in]  bbox    bounding box coordinates
 * @param [in]  qbox    quadrant coordinates
 * @param [out] mask    bit-mask indicating face overlap
 *                      * 8 bits:
 *                      *  bit 0: xlo face  bit 1: xhi face
 *                      *  bit 2: ylo face  bit 3: yhi face
 *                      *  bit 4: zlo face  bit 5: zhi face
 *                      *  bit 6: free      bit 7: free
 */
void p4est_box_intersection_mask_set(box_t bbox,box_t qbox,char *mask);

/** Search utility function for searching for points in unstructured using Newton's method on linear quads
 *
 * @param [in] p4est        p4est tree structure
 * @param [in] which_tree   p4est tree of this quadrant
 * @param [in] quad         p4est quadrant
 * @param [in] local_num    local number of this quadrant
 * @param [inout] pt_in     point to check intersection with this quadrant
 * @return  flag indicator for box and quadrant intersection
 */
int p4est_utilities_unst_pt_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                        p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                        void *pt_in);

/** Search utility function for searching for points in structured quads
 * 
 * @param [in] p4est        p4est tree structure
 * @param [in] which_tree   p4est tree of this quadrant
 * @param [in] quad         p4est quadrant
 * @param [in] local_num    local number of this quadrant
 * @param [inout] pt_in     point to check intersection with this quadrant
 * @return  flag indicator for box and quadrant intersection 
 */
int p4est_utilities_structured_pt_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                              p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                              void *pt_in);

/** Search utility function for overlapping boxes: uses bbox_t data structure
 *
 * @param [in] p4est        p4est tree structure
 * @param [in] which_tree   p4est tree of this quadrant
 * @param [in] quad         p4est quadrant
 * @param [in] local_num    local number of this quadrant
 * @param [inout] box_in    bounding box to check intersection with this quadrant
 * @return  flag indicator for box and quadrant intersection
 */
int p4est_utilities_unst_bbox_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                          p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                          void *box_in);

/** Search utility function for overlapping boxes: uses bbox2_t data structure
 *
 * @param [in] p4est        p4est tree structure
 * @param [in] which_tree   p4est tree of this quadrant
 * @param [in] quad         p4est quadrant
 * @param [in] local_num    local number of this quadrant
 * @param [inout] box_in    bounding box to check intersection with this quadrant
 * @return  flag indicator for box and quadrant intersection
 */
int p4est_utilities_unst_bbox2_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                           p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                           void *box_in);

/** Search utility function for overlapping boxes
 *
 * @param [in] p4est        p4est tree structure
 * @param [in] which_tree   p4est tree of this quadrant
 * @param [in] quad         p4est quadrant
 * @param [in] local_num    local number of this quadrant
 * @param [inout] box_in    bounding box to check intersection with this quadrant
 * @return  flag indicator for box and quadrant intersection
 */
int p4est_utilities_structured_bbox_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                                p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                                void *box_in);

/** Search utility function for boundary-box surface
 *
 * @param [in] p4est        p4est tree structure
 * @param [in] which_tree   p4est tree of this quadrant
 * @param [in] quad         p4est quadrant
 * @param [in] local_num    local number of this quadrant
 * @param [inout] box_in    bounding box to check intersection with this quadrant
 * @return  flag indicator for box and quadrant intersection
 */
int p4est_utilities_bbox_side_search_func(p4est_t *p4est,p4est_topidx_t which_tree,
                                          p4est_quadrant_t *quad,p4est_locidx_t local_num,
                                          void *box_in);

/** Search a list of points for the quadrants that contain the points
 *
 * @param [in] p4est    p4est tree structure
 * @param [in] npts     number of points to search (globally)
 * @param [in] xyz      point coordinates to search (dim*npts) [pt0_x,pt0_y,pt1_x,pt1_y,...]
 * @return list of the unique quadrants found
 */
quad_list_t * p4est_utilities_points_search(p4est_t *p4est,int npts,Real *xyz);

/** Search for the quadrants that contain a line
 *
 * @param [in] p4est    p4est tree structure
 * @param [in] pt_a     global coordinate of the starting point (x,y,{z})
 * @param [in] pt_b     global coordinate of the ending   point (x,y,{z})
 * @return list of the unique quadrants found
 */
quad_list_t * p4est_utilities_ray_search(p4est_t *p4est,Real *pt_a,Real *pt_b);

/** Search for the quadrants that intersect a bounding box
 *
 * @param [in] p4est    p4est tree structure
 * @param [in] bbox_lo  global coordinate of the lower bounding box point (x,y,{z})
 * @param [in] bbox_hi  global coordinate of the upper bounding box point (x,y,{z})
 * @return list of the unique quadrants found
 */
quad_list_t * p4est_utilities_bbox_search(p4est_t *p4est,Real *bbox_lo,Real *bbox_hi);

/** Search for the quadrants that intersect a bounding box surface
 *
 * @param [in] p4est    p4est tree structure
 * @param [in] bbox_lo  global coordinate of the lower bounding box point (x,y,{z})
 * @param [in] bbox_hi  global coordinate of the upper bounding box point (x,y,{z})
 * @return list of the unique quadrants found
 */
quad_side_list_t * p4est_utilities_bbox_surface_search(p4est_t *p4est,Real *bbox_lo,Real *bbox_hi);

/** Search for the quadrants that belong to the surface patch.
 *
 * @param geom          simulation geometry entity data
 * @param surface_name  surface patch entity name to search
 * @return list of the unique quadrants found
 */
quad_side_list_t * p4est_utilities_surface_name_list(p4est_t *p4est,
                                                     geometry_t *geom,
                                                     char *surface_patch_name);

/** Search for the quadrants that belong to the volume patches (assembled into single list).
 *
 * @param geom                  simulation geometry entity data
 * @param volume_patch_names    volume patch entity name(s) to search
 * @return list of the unique quadrants found
 */
quad_list_t * p4est_utilities_volume_name_list(p4est_t *p4est,
                                               geometry_t *geom,
                                               char *volume_patch_names);

#ifdef __cplusplus
}
#endif
#endif /* AMR_P4EST_UTILITIES_H */
