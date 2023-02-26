/**
 * File:   flow_solver.h
 * Author: akirby
 *
 * Created on February 25, 2023, 11:40 AM
 */

#ifndef FLOW_SOLVER_H
#define FLOW_SOLVER_H

/* system header files */
#include <dlfcn.h>
#include <stdio.h>
#include <vector>

/* header files */
#include "defs.h"

#ifndef P4_TO_P8
# include <p4est_extended.h>
#else
# include <p8est_extended.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* node index offset */
#define BASE 1

/* qpoint definitions */
#define XI 0
#define YI 1
#define ZI 2
#define VI 3
#define WALL_NODE_VOL  -2.0
#define OUTER_NODE_VOL -1.0

/* a struct to store flow solver function pointers and mesh */
class flow_solver {
private:
  const int numverts[4][6]={3,3,3,3,0,0,4,3,3,3,3,0,3,4,4,4,3,0,4,4,4,4,4,4};
  const int faceInfo[4][6][4]={1,2,3,0,1,4,2,0,2,4,3,0,1,3,4,0,0,0,0,0,0,0,0,0,
                               1,2,3,4,1,5,2,0,2,5,3,0,4,3,5,0,1,4,5,0,0,0,0,0,
                               1,2,3,0,1,4,5,2,2,5,6,3,1,3,6,4,4,6,5,0,0,0,0,0,
                               1,2,3,4,1,5,6,2,2,6,7,3,3,7,8,4,1,4,8,5,5,8,7,6};

public:
  char solver_so_file[buff_size];
  void *group_handle;
  char high_order_mode_flag;
  MPI_Fint new_group_comm;

  /* mandatory callback functions */
  void (*initialize_group_mpi)(int *group_comm);
  void (*initialize)(void);
  void (*set_pointers)(double**,
                       double**,
                       int**,
                       int**,
                       int**,
                       int**,
                       int**,
                       int**,
                       int**,
                       int**,
                       void**,
                       void**,
                       void**,
                       void**,
                       void**);
  void (*get_data)(int*,int*,int*,int*,int*,int*,int*,int*);
  void (*output_solution)(int *t);
  void (*point_inclusion)(int *npoint, double *x, int *cell_id);
  p4est_t* (*get_p4est)(void);

  /* pointers to data inside flow solvers */
  p4est_t *p4est;
  sc_array_t *qpoints;
  double* soln;
  double* xgeom;
  int* iblank;
  int* iwbcnode;
  int* iobcnode;
  int* ndc4;
  int* ndc5;
  int* ndc6;
  int* ndc8;
  int* iblank_cell;
  double translation[3];
  double *xgeom_translated;

  /* node and cell resolution to force p4est_overset to connect */
  double* node_res;
  double* cell_res;

  /* function pointers */
  void *count_receptor_nodes;
  void *create_receptor_nodes;
  void *donor_inclusion_test;
  void *create_donor_frac;
  void *convert_to_receptor_coefficients;

  /* static data that may need to get updated */
  int nnode;
  int nwbc;
  int nobc;
  int ntet;
  int nprism;
  int npyramid;
  int nhex;
  int body_tag;

  int ncell_types = 4;
  int kstride4 = 4;
  int kstride5 = 5;
  int kstride6 = 6;
  int kstride8 = 8;
  int nfield = 5;

  igbp_t igbp;
  obc_t obc;

  flow_solver() :
    group_handle{NULL},
    high_order_mode_flag{0},
    p4est{NULL},
    qpoints{NULL},
    soln{NULL},
    xgeom{NULL},
    iblank{NULL},
    iwbcnode{NULL},
    iobcnode{NULL},
    ndc4{NULL},
    ndc5{NULL},
    ndc6{NULL},
    ndc8{NULL},
    iblank_cell{NULL},
    xgeom_translated{NULL},
    node_res{NULL},
    cell_res{NULL},
    nnode{0},
    nwbc{0},
    nobc{0},
    ntet{0},
    nprism{0},
    npyramid{0},
    nhex{0},
    body_tag{0} {}

 ~flow_solver(){
    flow_close_dynamic_library();

    if(cell_res) delete [] cell_res;
    if(node_res) delete [] node_res;
    if(xgeom_translated) delete [] xgeom_translated;
 }

  /* flow solver routines */
  void flow_load_dynamic_library();
  void flow_close_dynamic_library();
  void flow_initialize_group_mpi(MPI_Comm group_comm);
  void flow_set_pointers();
  void flow_output_solution(int t);
  void flow_point_inclusion(int *npoint, double *x, int *cell_id);
};

#ifdef __cplusplus
}
#endif
#endif /* FLOW_SOLVER_H */