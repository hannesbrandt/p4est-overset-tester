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

/* element definitions */
#define TRI   0
#define QUAD  1
#define TET   2
#define PYR   3
#define PRISM 4
#define HEX   5

/* a struct to store flow solver function pointers and mesh */
class flow_solver {
private:

  MPI_Comm mpicomm;
  int rank;

  const int numverts[6][6]={2,2,2,0,0,0, // tri
                            2,2,2,2,0,0, // quad
                            3,3,3,3,0,0, // tet
                            4,3,3,3,3,0, // pyramid
                            3,4,4,4,3,0, // prism
                            4,4,4,4,4,4};// hex

                            /* Face 1, Face 2, Face 3, Face 4, Face 5, Face 6 */
  const int faceInfo[6][6][4]={1,2,0,0,2,3,0,0,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // tri
                               1,3,0,0,2,4,0,0,1,2,0,0,3,4,0,0,0,0,0,0,0,0,0,0, // quad
                               1,2,3,0,1,4,2,0,2,4,3,0,1,3,4,0,0,0,0,0,0,0,0,0, // tet
                               1,2,3,4,1,5,2,0,2,5,3,0,4,3,5,0,1,4,5,0,0,0,0,0, // pyramid
                               1,2,3,0,1,4,5,2,2,5,6,3,1,3,6,4,4,6,5,0,0,0,0,0, // prism
                               1,2,4,3,1,5,6,2,2,6,8,4,4,8,7,3,1,3,7,5,5,7,8,6};// hex

  double triangleArea(double ax,double ay,double bx,double by,double cx,double cy){
    return 0.5*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by));
  }
  double scalarProduct(double a1,double a2,double a3,double b1,double b2,double b3,double c1,double c2,double c3){
    return (a1*b2*c3 - a1*b3*c2 + a2*b3*c1 - a2*b1*c3 + a3*b1*c2 - a3*b2*c1);
  }

  double cellArea(double xc[8][3],const int numverts[6],const int fconn[6][4],int nfaces){
    double area = 0.0;
    if (nfaces == 3) {
      double ax = xc[0][0]; double ay = xc[0][1];
      double bx = xc[1][0]; double by = xc[1][1];
      double cx = xc[2][0]; double cy = xc[2][1];

      area = triangleArea(ax,ay,bx,by,cx,cy);
    } else {
      double ax = xc[0][0]; double ay = xc[0][1];
      double bx = xc[1][0]; double by = xc[1][1];
      double cx = xc[2][0]; double cy = xc[2][1];
      double dx = xc[3][0]; double dy = xc[3][1];

      area  = triangleArea(ax,ay, // node 1
                           bx,by, // node 2
                           cx,cy);// node 3
      area += triangleArea(bx,by, // node 2
                           dx,dy, // node 4
                           cx,cy);// node 3
    }
    return area;
  }

  double cellVolume(double xc[8][3],const int numverts[6],const int fconn[6][4],int nfaces){
    double vol = 0.0;
    int X=0,Y=1,Z=2;
    int N0=0,N1=1,N2=2,N3=3;

      for (int iface = 1; iface <= nfaces; iface++) {
        int face = iface-BASE;

        if (numverts[face] == 3) {
          vol-= 0.50*scalarProduct(xc[fconn[face][N0]-BASE][X], xc[fconn[face][N0]-BASE][N1], xc[fconn[face][N0]-BASE][N2],
                                   xc[fconn[face][N1]-BASE][X], xc[fconn[face][N1]-BASE][N1], xc[fconn[face][N1]-BASE][N2],
                                   xc[fconn[face][N2]-BASE][X], xc[fconn[face][N2]-BASE][N1], xc[fconn[face][N2]-BASE][N2]);
        }else{
          vol-= 0.25*scalarProduct(xc[fconn[face][N0]-BASE][X], xc[fconn[face][N0]-BASE][Y], xc[fconn[face][N0]-BASE][Z],
                                   xc[fconn[face][N1]-BASE][X], xc[fconn[face][N1]-BASE][Y], xc[fconn[face][N1]-BASE][Z],
                                   xc[fconn[face][N2]-BASE][X], xc[fconn[face][N2]-BASE][Y], xc[fconn[face][N2]-BASE][Z]);
          vol-= 0.25*scalarProduct(xc[fconn[face][N0]-BASE][X], xc[fconn[face][N0]-BASE][Y], xc[fconn[face][N0]-BASE][Z],
                                   xc[fconn[face][N2]-BASE][X], xc[fconn[face][N2]-BASE][Y], xc[fconn[face][N2]-BASE][Z],
                                   xc[fconn[face][N3]-BASE][X], xc[fconn[face][N3]-BASE][Y], xc[fconn[face][N3]-BASE][Z]);
          vol-= 0.25*scalarProduct(xc[fconn[face][N0]-BASE][X], xc[fconn[face][N0]-BASE][Y], xc[fconn[face][N0]-BASE][Z],
                                   xc[fconn[face][N1]-BASE][X], xc[fconn[face][N1]-BASE][Y], xc[fconn[face][N1]-BASE][Z],
                                   xc[fconn[face][N3]-BASE][X], xc[fconn[face][N3]-BASE][Y], xc[fconn[face][N3]-BASE][Z]);
          vol-= 0.25*scalarProduct(xc[fconn[face][N1]-BASE][X], xc[fconn[face][N1]-BASE][Y], xc[fconn[face][N1]-BASE][Z],
                                   xc[fconn[face][N2]-BASE][X], xc[fconn[face][N2]-BASE][Y], xc[fconn[face][N2]-BASE][Z],
                                   xc[fconn[face][N3]-BASE][X], xc[fconn[face][N3]-BASE][Y], xc[fconn[face][N3]-BASE][Z]);
        }
      }
    vol /= 3.0;
    return vol;
  }

  double computeCellVolume(double xv[8][3],int cellType,int dim){
    double vol;
    int nfaces;

    switch(cellType){
      case TRI:   nfaces=3; break;
      case QUAD:  nfaces=4; break;
      case TET:   nfaces=4; break;
      case PYR:   nfaces=5; break;
      case PRISM: nfaces=5; break;
      case HEX:   nfaces=6; break;
    }
    vol = (dim==2) ? cellArea(xv,numverts[cellType],faceInfo[cellType],nfaces):
                   cellVolume(xv,numverts[cellType],faceInfo[cellType],nfaces);
    return vol;
  }

public:
  char solver_so_file[buff_size];
  void *group_handle;
  char high_order_mode_flag;
  MPI_Fint new_group_comm;

  /* mandatory callback functions */
  void (*initialize_group_mpi)(int *group_comm);
  void (*initialize)(void);
  void (*set_pointers)(double **soln,
                       double **geom,
                       int **iblank,
                       int **wbcnode,
                       int **obcnode,
                       int **ndcTri,
                       int **ndcQuad,
                       int **ndcTet,
                       int **ndcPyr,
                       int **ndcPrism,
                       int **ndcHex,
                       int **iblank_cell,
                       void **count_receptor_nodes_callback,
                       void **create_receptor_nodes_callback,
                       void **donor_inclusion_test_callback,
                       void **create_donor_frac_callback,
                       void **convert_to_receptor_coefficients_callback);
  void (*get_data)(int *body_tag,
                   int *nnode,
                   int *nwbc,
                   int *nobc,
                   int *ntri,
                   int *nquad,
                   int *ntet,
                   int *npyramid,
                   int *nprism,
                   int *nhex);
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
  int* ndcTri;
  int* ndcQuad;
  int* ndcTet;
  int* ndcPyr;
  int* ndcPrism;
  int* ndcHex;
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
  int ntri;
  int nquad;
  int ntet;
  int nprism;
  int npyramid;
  int nhex;
  int body_tag;

  int ncell_types  = 6;
  int kstrideTri   = 3;
  int kstrideQuad  = 4;
  int kstrideTet   = 4;
  int kstridePyr   = 5;
  int kstridePrism = 6;
  int kstrideHex   = 8;
  int nfield = 1;

  igbp_t igbp;
  obc_t obc;

  flow_solver(MPI_Comm _comm,int _rank) :
    mpicomm{_comm},
    rank{_rank},
    group_handle{NULL},
    high_order_mode_flag{0},
    p4est{NULL},
    qpoints{NULL},
    soln{NULL},
    xgeom{NULL},
    iblank{NULL},
    iwbcnode{NULL},
    iobcnode{NULL},
    ndcTri{NULL},
    ndcQuad{NULL},
    ndcTet{NULL},
    ndcPyr{NULL},
    ndcPrism{NULL},
    ndcHex{NULL},
    iblank_cell{NULL},
    xgeom_translated{NULL},
    node_res{NULL},
    cell_res{NULL},
    nnode{0},
    nwbc{0},
    nobc{0},
    ntri{0},
    nquad{0},
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
  void flow_set_pointers(int group_index);
  void flow_output_solution(int t);
  void flow_point_inclusion(int *npoint, double *x, int *cell_id);
};

#ifdef __cplusplus
}
#endif
#endif /* FLOW_SOLVER_H */