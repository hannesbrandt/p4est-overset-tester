/**
 * \file    dg4est_solver.hxx
 * \author  akirby
 */

#ifndef DG4EST_SOLVER_H
#define DG4EST_SOLVER_H

/* header files */
#include "precision_types.h"
#include "memory.hxx"

/* system header files */
#ifndef DOXYGEN_IGNORE
#  include <mpi.h>
#  include <stdio.h>
#  include <stdlib.h>
#  include <stddef.h>
#  include <sys/stat.h>
#  include <memory>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define VERSION     301     /**< @brief dg4est_solver.hxx version number for compatibility */

#define BUFF_SIZE   1024    /**< @brief string buffer size */
#define MAX_LEVELS  19      /**< @brief number of maximum levels allowed */
#define MAX_QDEGREE 10      /**< @brief maximum p-degree allowed */

#define BC_TYPE   1
#define FULL_TYPE 2
#define HANG_TYPE 3

/* offset for geometry \face_nodes_info */
#define EFENT_IND  0    /**< @brief Element face entity tag */
#define EFTAG_IND  1    /**< @brief Element face tag */
#define EFNDS_IND  2    /**< @brief Element face node list */

/* patch id or tag not assigned */
#define NO_PATCH   -1
#define NO_TAG     -1

/* wake3d boundary tag */
#define W3D_NO     0
#define W3D_WBC    1
#define W3D_OBC    2

/* max number of bcs to iterate */
#define MAX_NBC    10

/* face indicator bit masks */
#define xlo_mask 0b00000001 /**< xlo face intersection flag */
#define xhi_mask 0b00000010 /**< xhi face intersection flag */
#define ylo_mask 0b00000100 /**< ylo face intersection flag */
#define yhi_mask 0b00001000 /**< yhi face intersection flag */
#define zlo_mask 0b00010000 /**< zlo face intersection flag */
#define zhi_mask 0b00100000 /**< zhi face intersection flag */

/**
 * @brief Iblank values type.
 */
typedef enum {
  receptor_mask = -1, /**< receptor cell */
  hole_mask     =  0, /**< hole cell */
  field_mask    =  1  /**< valid cell */
}
iblank_names;

/**
 * @brief Generic quadrant list type.
 */
typedef struct {
    int nquads;
    int *quad_list;
}
quad_list_t;

/**
 * @brief Generic quadrant list type with side intersection (bit-mask) info.
 */
typedef struct {
    int nsides;
    int nquads;
    int *quad_list;
    char *side_list;
}
quad_side_list_t;

/**
 * @brief Quadrant physical properties type.
 */
typedef struct {
    int dim;
    int level;
    int qdegree;
    int unstructured;

    Real *quad_dx;
    Real *quad_xyz;
}
quad_properties_t;

/**
 * @brief Quadrant user data type.
 */
typedef struct {
    char bc[6];
    char type;
    int qdegree;

    int geom_size;
    Real *quad_geom;
}
quad_user_data_t;

typedef struct {
    int amr_nhboxes;                    /**< number of box regions for level refinement */
    dg::memory<Real> amr_hboxes_lo;     /**< box region lo coordinates for level refinement */
    dg::memory<Real> amr_hboxes_hi;     /**< box region hi coordinates for level refinement */
    dg::memory<int> amr_hboxes_level;   /**< box region refinement level */
    dg::memory<int> amr_hboxes_autotag; /**< flag for autotagging all elements in box region */
}
amr_hbox_t;

/**
 * @brief Contains all WAKE3D data.
 */
class wake3d_t {
  public:
    int nobc;           /**< Number of overset bc nodes */
    int nwbc;           /**< Number of wall bc nodes */
    int nnodes;         /**< Number of nodes */

    Real wake3d_dt;     /**< size of WAKE3D time step */
    Real igbp_scale;    /**< inter-grid boundary point scaling constant */
    Real igbp_scale_cf; /**< coarse level scaling factor (scales igbp_scale for lower levels) */

    dg::memory<int> ndcTri;      /**< Triangle cell connectivities */
    dg::memory<int> ndcQuad;     /**< Quad cell connectivities */
    dg::memory<int> ndcTet;      /**< Tetrahedra cell connectivities */
    dg::memory<int> ndcPyr;      /**< Pyramid cell connectivities */
    dg::memory<int> ndcPrism;    /**< Prism cell connectivities */
    dg::memory<int> ndcHex;      /**< Hexahedron cell connectivities */
    dg::memory<int> iblank;      /**< Node iblank */
    dg::memory<int> iobcnode;    /**< List of outer bc nodes */
    dg::memory<int> iwbcnode;    /**< List of wall bc nodes */
    dg::memory<int> iblank_cell; /**< Cell iblank */
    dg::memory<Real> xgeom;      /**< Unique node geometry coordinates */

    /* constructor */
    wake3d_t() :
        nobc{0},nwbc{0},nnodes{0},
        wake3d_dt{1.0},igbp_scale{1.0},igbp_scale_cf{1.0} {}

    /* destructor */
   ~wake3d_t()=default;
};

/**
 * @brief Contains all grid data.
 */
typedef struct {
    int nlevels;             /**< Number of AMR levels */
    int min_level;           /**< Minumum level of grid refinement */
    int max_level;           /**< Maximum level of grid refinement */
    int max_level_pmax;      /**< Maximum level allowing pdegree_max */
    int max_level_global;    /**< Maximum level override if hbox level > max_level */
    int nelem[3];            /**< Level-0 number of elements in mesh */
    int periodic[3];         /**< Periodic boundary condition indicator */
    int construct_grid;      /**< Initialize grid by building up levels */
    int total_quads;         /**< Total Mesh Quadrants */
    unsigned int total_dofs; /**< Total Degrees of Freedom */
    Real xlo[3];             /**< Lower coordinates of the grid domain */
    Real xhi[3];             /**< Upper coordinates of the grid domain */
    Real min_dx[3];          /**< Size of the finest element */
    Real max_dx[3];          /**< Size of the coarsest element */
    Real domain_volume;      /**< Grid domain volume */
}
grid_t;

/**
 * @brief Contains all simulation data.
 */
typedef struct {
    int dim;                                /**< Spatial dimension of the grid */
    int regrid_hboxes;                      /**< Regrid based on h-amr box regions flag */
    int regrid_coarsen;                     /**< Regrid coarsening flag */
    int regrid_nspread;                     /**< Regrid spread layers count */
    int regrid_buffer;                      /**< Regrid buffer flag */
    int regrid_interval;                    /**< Interval to regrid solution */
    int regrid_features;                    /**< Regrid based on features flag */
    int regrid_strategy;                    /**< Regrid stategy flag */
    int mapped_grid_flag;                   /**< Mapped Carteisian grid indicator */
    int unstructured_flag;                  /**< Unstructured grid indicator */
    int visualization_interval;             /**< Interval to visualize solution */
    char input_file[BUFF_SIZE];             /**< Input file name */
    char visualization_path[BUFF_SIZE];     /**< Visualization output directory */
    time_t input_file_mod_time;             /**< Input file modification time */
}
simulation_t;

/**
 * @brief Contains boundary condition data.
 */
typedef struct {
    int  npatch;                          /**< Number of the physical patches defined */
    Uint nentity_face;                    /**< Number of boundary geometry entities (edges/faces) */
    Uint nentity_vol;                     /**< Number of volume geometry entities (faces/volumes) */
    Uint nelem_face;                      /**< Number of element faces in the grid_file */
    char wake3d_wbc[BUFF_SIZE];           /**< Char containing the patches name of wall boundaries */
    char wake3d_obc[BUFF_SIZE];           /**< Char containing the patches name of outern boundaries */
    char **patch_name;                    /**< List of physical patches' names (size npatch)*/
    dg::memory<int> patch_dim;            /**< List of physical patches's dimension 2: 1: 1D, 2D, 3: 3D (size npatch) */
    dg::memory<int> entity_face_patchtag; /**< List of boundary geometry entities' patch tag (edges/faces) */
    dg::memory<int> entity_vol_patchtag;  /**< List of volume geometry entities' patch tag (faces/volumes) */
    dg::memory<Uint> entity_face_tag;     /**< List of boundary geometry entities' tag (edges/faces) */
    dg::memory<Uint> entity_vol_tag;      /**< List of volume geometry entities' tag (faces/volumes) */
    dg::memory<Uint> face_nodes_info;     /**< List of element face (edges/faces) info */
    dg::memory<int> wake3d_bc;            /**< List of indeces tagging patches as wall[1], outer[2] or null[0] (size npatch)*/
}
geometry_t;

/**
 * @brief Contains gridfile line number of elements (used to retrieve HiO info faster)
 */
typedef struct {
    char gridfile_name[BUFF_SIZE]; /**< Unstructured grid file name */
    int binary;                    /**< Binary file (1: yes, 0: no) */
    int swap;                      /**< Binary swap needed (1:yes, 0: no) */
    int file_format;               /**< File format (41 or 2.0) */
    size_t entities_bytes;         /**< number of bytes for gmsh $ENTITIES section (only for binary) */
    size_t nodes_bytes;            /**< number of bytes for gmsh $NODES section (only for binary) */
    int nvertices;                 /**< number of vertices per high-order volume element (only for binary) */
    dg::memory<Uint> elem_tag;     /**< Element volume tag */
    dg::memory<Uint> p4est_part;   /**< p4est partition mapping */
}
gridfile_t;

/**
 * @brief Contains all context data of DG4EST.
 *        This contains all the other data structures.
 */
class ctx_t {
  public:
    int version;                /**< Header file version flag */
    int rank;                   /**< MPI rank */
    int nranks;                 /**< Number of ranks in this communicator */
    MPI_Comm comm;              /**< MPI communicator */
    FILE *log_io;               /**< Output file stream */
    int log_info;               /**< Logging threshold flag:
                                 *      *Values:
                                 *      * 0     log everything
                                 *      * 1     prefix file and line number
                                 *      * 2     information on the internal state
                                 *      * 3     information on conditions, decisions
                                 *      * 4     main information a function is doing
                                 *      * 5     important consistency/performance information
                                 *      * 6     few lines for a major API function
                                 *      * 7     logs a few lines max per program
                                 *      * 8     logs errors only
                                 *      * 9     never logs anything */
    int initial;                /**< Simulation initialization flag */
    int visualize;              /**< Visualize solution flag */
    int nvisualize;             /**< Visualize count */
    int high_order_viz;         /**< Visualize high-order flag (VTK_LAGRANGE_XXX) */
    int high_order_viz_npts;    /**< Visualize high-order solution points (max 10) */

    grid_t d_grid;              /**< Grid data */
    wake3d_t d_wake3d;          /**< WAKE3D data */
    amr_hbox_t d_amr_hboxes;    /**< AMR level box region refinement data */
    simulation_t d_simulation;  /**< Simulation data */
    geometry_t d_geometry;      /**< Geometry data */
    gridfile_t d_gridfile;      /**< Grifile data */

    /* constructors */
    ctx_t()=default;
   ~ctx_t()=default;
};

typedef struct {
    int p_levels;         /**< Number of multigrid levels */
}
multigrid_t;

/**
 * @brief Contains all data of DG4EST.
 */
class dg4est_t {
  public:
    ctx_t d_ctx;                       /**< Context dg4est data */
    multigrid_t d_multigrid;           /**< Multigrid data */

    /* constructors */
    dg4est_t()=default;
   ~dg4est_t()=default;
};

/** Utility function for deleting quad_list_t data structure.
 *
 * @param [inout] list  quad_list_t data to delete
 */
static inline
void dg4est_destroy_quad_list(quad_list_t **list){
    if((*list)[0].quad_list) free((*list)[0].quad_list);
    free(*list); *list = NULL;
}

/** Utility function for deleting quad_side_list_t data structure.
 *
 * @param [inout] list  quad_side_list_t data to delete
 */
static inline
void dg4est_destroy_quad_side_list(quad_side_list_t **list){
    if((*list)[0].quad_list) free((*list)[0].quad_list);
    if((*list)[0].side_list) free((*list)[0].side_list);
    free(*list); *list = NULL;
}

/** External wrapper function for initializing quadrant solution data
 *
 * @param [in]  quad_prop   AMR quadrant properties
 * @param [out] user_data   AMR quadrant user data
 * @param [in]  init_strategy   initialization solution flag
 * @param [in]  time_init       initialization time
 */
void dg4est_solver_initialize_quadrant(quad_properties_t *quad_prop,
                                       quad_user_data_t *user_data,
                                       int *init_strategy,
                                       Real time_init);

///** External wrapper function for initializing mesh geometry data
// *
// * @param [inout]  ext  AMR external data properties
// */
//void dg4est_solver_initialize_geometry(external_t *ext);
//
///** External wrapper function for initializing boundary condition data
// *
// * @param [inout]  ext  AMR external data properties
// */
//void dg4est_solver_initialize_bc_data(external_t *ext,ctx_t *ctx);

/** External wrapper function of amr tag feature function
 *
 * @param [in]  quad_prop       AMR quadrant properties
 * @param [in]  user_data       AMR quadrant user data
 * @param [out] tag             AMR tag indicator (0=do not tag, 1=tag)
 */
void dg4est_solver_tag_feature(quad_properties_t *quad_prop,
                               quad_user_data_t *user_data,
                               int *tag);

/** Wrapper function for external coarsen operator
 *
 * @param [out] parent_user_data    AMR parent quadrant user data
 * @param [in]  children_user_data  AMR children quadrant user data
 * @param [in]  initial             AMR initial solution
 */
void dg4est_solver_coarsen_operator(quad_user_data_t *parent_user_data,
                                    quad_user_data_t *children_user_data,
                                    ctx_t *ctx);

/** Wrapper function for external refine operator
 *
 * @param [in]  parent_user_data    AMR parent quadrant user data
 * @param [out] children_user_data  AMR children quadrant user data
 */
void dg4est_solver_refine_operator(quad_user_data_t *parent_user_data,
                                   quad_user_data_t *children_user_data,
                                   ctx_t *ctx);

/** Wrapper function for external change q-degree operator
 *
 * @param [in] quad_prop    AMR quadrant properties
 * @param [in] user_data    AMR quadrant user data
 */
void dg4est_solver_change_qdegree(quad_properties_t *quad_prop,
                                  quad_user_data_t *quad_data);

/** External wrapper function for projecting the geometry coordinates to set of reference points
 *
 * @param [in]  npts1D       number of 1D projection point locations
 * @param [in]  ref_pts      1D point locations in reference space
 * @param [in]  geom_xyz_in  AMR quadrant coordinates
 * @param [out] geom_xyz_out AMR projected quadrant coordinates
 */
void dg4est_solver_project_geometry(int *npts1D,
                                    Real *ref_pts,
                                    Real *geom_xyz_in,
                                    Real *geom_xyz_out);

/** External wrapper function for integrating a function f over a quadrant
 *
 * @param [in]  dim         dimension of the simulation
 * @param [in]  nfields     number of fields to integrate per quadrature point
 * @param [in]  pdegree     AMR quadrant p-degree
 * @param [in]  quad_geom   AMR quadrant geometry
 * @param [in]  volume      AMR quadrant volume
 * @param [in]  f           function to integrate over quadrant
 * @param [out] u           integrated value of f over quadrant
 */
void dg4est_solver_integrate_volume(int *dim,int *nfields,int *pdegree,
                                    Real *quad_geom,Real *volume,Real *f,Real *u);

/** External wrapper function for computing quadrant volume
 *
 * @param [in]  pdegree     AMR quadrant p-degree
 * @param [in]  quad_geom   AMR quadrant geometry
 * @param [out] volume      quadrant volume
 */
void dg4est_solver_compute_volume(int *pdegree,Real *geom_xyz,Real *volume);

/** External wrapper function for creating the receptor nodes for Tioga
 *
 * @param [in]  pdegree         AMR quadrant p-degree
 * @param [in]  quad_geom       AMR quadrant geometry
 * @param [out] receptor_nodes  receptor nodes for Tioga
 */
void dg4est_solver_create_receptor_nodes(int *pdegree,Real *quad_geom,
                                         Real *receptor_nodes);

/** External wrapper function for creating the donor fraction weights for Tioga
 *
 * @param [in]  pdegree     AMR quadrant p-degree
 * @param [in]  tm          total modes
 * @param [in]  r           natural coordinate #1
 * @param [in]  s           natural coordinate #2
 * @param [in]  t           natural coordinate #3
 * @param [out] frac        donor fraction weights
 */
void dg4est_solver_create_donor_frac(int *pdegree,int *tm,
                                     Real *r,Real *s,Real *t,
                                     Real *frac);

/** External wrapper function for converting receptor coefficients to solutions for Tioga
 *
 * @param [in]  pdegree         AMR quadrant p-degree
 * @param [in]  nfields         number of fields per solution point
 * @param [in]  nqp             number of solution points 1D
 * @param [in]  tm              total modes
 * @param [in]  quad_geom       quadrant geometry coordinates
 * @param [in]  coefficients    solution coefficients
 * @param [out] soln            solution values
 */
void dg4est_solver_convert_to_receptor_coefficients(int *pdegree,int *nfields,int *nqp,int *tm,
                                                    Real *quad_geom,Real *coefficients,Real *soln);

/** External wrapper function for donor inclusion for Tioga
 *
 * @param [in]  pt_xyz      physical point to test cell inclusion
 * @param [in]  pdegree     AMR quadrant pdegree
 * @param [in]  quad_xyz    AMR quadrant geometry
 * @param [out] pass_flag   flag indicator if point found in quadrant
 * @param [out] rst_out     natural coordinates of point found in quadrant
 */
void dg4est_solver_donor_inclusion_test(Real *pt_xyz,int pdegree,Real *quad_xyz,
                                        int *pass_flag,Real *rst_out);

#ifdef __cplusplus
}
#endif
#endif /* DG4EST_SOLVER_H */