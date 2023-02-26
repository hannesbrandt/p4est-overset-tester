/**
 * File:   flow_solver.cxx
 * Author: akirby
 *
 * Created on February 25, 2023, 11:40 AM
 */

#include <cmath>

/* header files */
#include "flow_solver.h"

void flow_solver::flow_load_dynamic_library(){
  char *filename  = solver_so_file;
  char *error;

  /* get group handle */
  group_handle = dlopen (filename, RTLD_LAZY);
  if(!group_handle) {fputs (dlerror(), stderr);}

  /* load functions */
  *(void **)(&initialize_group_mpi) = dlsym(group_handle,"driver_interface_initialize_group_mpi");
  if((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&initialize) = dlsym(group_handle, "driver_interface_initialize");
  if((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&get_data) = dlsym(group_handle, "driver_interface_get_data");
  if((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&set_pointers) = dlsym(group_handle, "driver_interface_set_pointers");
  if((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&output_solution) = dlsym(group_handle, "driver_interface_output_solution");
  if((error = dlerror()) != NULL) {fputs(error, stderr);}

  *(void **)(&point_inclusion) = NULL;
  *(void **)(&point_inclusion) = dlsym(group_handle, "driver_interface_point_inclusion");

  *(void **)(&get_p4est) = NULL;
  *(void **)(&get_p4est) = dlsym(group_handle, "driver_interface_get_p4est");
}

void flow_solver::flow_close_dynamic_library(){
  if(group_handle) dlclose(group_handle);
}

void flow_solver::flow_initialize_group_mpi(MPI_Comm group_comm){
  new_group_comm = MPI_Comm_c2f(group_comm);
  initialize_group_mpi(&new_group_comm);
}

void flow_solver::flow_set_pointers(){
  int inode[8];
  double xv[8][3];
  double vol;

  int wbc,obc;
  int i,k,m,n;
  int nvert;
  double *pt;

  set_pointers(&soln,
               &xgeom,
               &iblank,
               &iwbcnode,
               &iobcnode,
               &ndc4,
               &ndc5,
               &ndc6,
               &ndc8,
               &iblank_cell,
               &count_receptor_nodes,
               &create_receptor_nodes,
               &donor_inclusion_test,
               &create_donor_frac,
               &convert_to_receptor_coefficients);

  get_data(&body_tag,
           &nnode,
           &nwbc,
           &nobc,
           &ntet,
           &npyramid,
           &nprism,
           &nhex);

  /* translate the geometry using translation from input file */
  if(xgeom_translated) delete [] xgeom_translated;
  xgeom_translated = new double[3*nnode];

  for(i=0; i<nnode; ++i){
    for(int j=0;j<3;j++){
      xgeom_translated[3*i+j] = xgeom[3*i+j] + translation[j];
    }
  }

  /* calculate approximate node resolutions */
  std::vector<double> nodeRes(nnode,0.0);
  std::vector<int> iflag(nnode,0);

  int nv[ncell_types] = {kstride4,kstride5,kstride6,kstride8};
  int nc[ncell_types] = {ntet,npyramid,nprism,nhex};
  int *vconn[ncell_types] = {ndc4,ndc5,ndc6,ndc8};

  // loop cell types
  for(n=0,k=0; n<ncell_types; n++){
    nvert = nv[n];
    // loop cells of each type
    for(i=0; i<nc[n]; i++){
      // loop vertices of each cell
      for(m=0; m<nvert; m++){
        inode[m] = vconn[n][nvert*i+m]-BASE;
        iflag[inode[m]]++;

        const double *node_xgeom = &xgeom[3*inode[m]];
        xv[m][XI] = node_xgeom[XI];
        xv[m][YI] = node_xgeom[YI];
        xv[m][ZI] = node_xgeom[ZI];
      }

      vol = 1.0; //TODO: computeCellVolume(xv,nvert);
      for(m=0; m<nvert; m++) nodeRes[inode[m]] += vol;
    }
  }
  // nodal resolution: average of all the cells associated with it
  for(i=0; i<nnode; i++) nodeRes[i] /= iflag[i];

  /* set p4est info */
  if(get_p4est != NULL) p4est = get_p4est();

  /* set qpoints */
  if(qpoints) sc_array_destroy(qpoints);
  qpoints = sc_array_new_count(4*sizeof(double),nnode);

  wbc = obc = 0;
  for(i=0; i<nnode; ++i){
    pt = (double *) sc_array_index(qpoints,i);

    /* set point coordinates */
    pt[XI] = xgeom_translated[3*i+XI];
    pt[YI] = xgeom_translated[3*i+YI];
    pt[ZI] = xgeom_translated[3*i+ZI];

    /* set volume associated with node */
    pt[VI] = nodeRes[i];
    if(iwbcnode[wbc] == i) {pt[VI] = WALL_NODE_VOL; wbc++; continue;}
    if(iobcnode[obc] == i) {pt[VI] = OUTER_NODE_VOL; obc++; continue;}
  }
}

void flow_solver::flow_output_solution(int t){
    output_solution(&t);
}

void flow_solver::flow_point_inclusion(int *npoint, double *x, int *cell_id){
  if(point_inclusion != NULL){
    point_inclusion(npoint,x,cell_id);
  } else {
    for(int i=0;i<*npoint;++i) cell_id[i] = -1;
  }
}