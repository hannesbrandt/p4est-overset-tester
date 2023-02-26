/**
 * File:   p4est_overset_interface.h
 * Author: akirby
 *
 * Created on February 25, 2023, 11:46 AM
 */

#ifndef P4EST_OVERSET_INTERFACE_H
#define P4EST_OVERSET_INTERFACE_H

/* header files */
#include "defs.h"
#include "flow_solver.h"

#ifndef P4_TO_P8
# include <p4est_extended.h>
# include <p4est_search.h>
# include <p4est_vtk.h>
# include <p4est_bits.h>
# include "p4est_multi_overset.h"
#else
# include <p8est_extended.h>
# include <p8est_search.h>
# include <p8est_vtk.h>
# include <p8est_bits.h>
# include "p8est_multi_overset.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* p4est_overset to store function pointers and files */
class p4est_overset_interface {
public:
  p4est_overset_interface(){};
 ~p4est_overset_interface();

  /* p4est_overset wrapper routines */
  void p4est_overset_init(MPI_Comm glocomm,
                          MPI_Comm leadercomm,
                          MPI_Comm rolecomm,
                          int mesh_index,
                          int num_meshes,
                          const int *mesh_offsets,
                          p4est_t *bgp4est,
                          sc_array_t *qpoints);
};

#ifdef __cplusplus
}
#endif
#endif /* P4EST_OVERSET_INTERFACE_H */