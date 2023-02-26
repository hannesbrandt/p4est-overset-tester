/**
 * File:   defs.h
 * Author: akirby
 *
 * Created on February 25, 2023, 12:10 PM
 */

#ifndef defs_h
#define defs_h

#define buff_size 1024
#include <mpi.h>
#include <sys/stat.h>

/* outer boundary condition points for tagging on the off-body */
class obc_t {
public:
  int nobc;   /**< number of obc points */
  double *x;  /**< x location of points */
  double *y;  /**< y location of points */
  double *z;  /**< z location of points */
  double *dx; /**< mesh size at points */

  obc_t() :
    nobc{0},
    x{NULL},
    y{NULL},
    z{NULL},
    dx{NULL} {}

  ~obc_t(){
      if(x) delete [] x;
      if(y) delete [] y;
      if(z) delete [] z;
      if(dx) delete [] dx;
  }
};

/** intergrid boundary points for tagging on the off-body */
class igbp_t {
public:
  /** variables that come from get igbp */
  int n_lcl;       /**< number of igbp on proc */
  int *index;      /**< index of igbp node on proc */
  double *dx_lcl;  /**< mesh size pointer on proc */

  /** variables that get accumulated and copied to all procs */
  int n;           /**< number of igbp on all procs */
  double *x;       /**< x location of igbp copied on all procs */
  double *y;       /**< y location of igbp copied on all procs */
  double *z;       /**< z location of igbp copied on all procs */
  double *dx;      /**< mesh size of igbp copied on all procs */

  igbp_t() :
    n_lcl{0},
    index{NULL},
    dx_lcl{NULL},
    n{0},
    x{NULL},
    y{NULL},
    z{NULL},
    dx{NULL} {}

  ~igbp_t(){
      if(x) delete [] x;
      if(y) delete [] y;
      if(z) delete [] z;
      if(dx) delete [] dx;
      n = 0;
  }
};

#endif /* defs_h */
