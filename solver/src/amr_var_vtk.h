/**
 * \file    amr_var_vtk.h
 * \ingroup amr_group
 * \author  akirby
 *
 */

#ifndef AMR_VAR_VTK_H
#define AMR_VAR_VTK_H

#ifndef P4_TO_P8
#  include <p4est_vtk.h>
#else
#  include <p8est_vtk.h>
#endif

/**
 * vtk_plot_fluid_t contains conservative fluid variable visualization-related data
 */
typedef struct {
    int nplotpts1D;
    int nplotpts;

    int high_order_viz;
    sc_array_t *positions;  /**< plot coordinates [3*nplotpts] */

    sc_array_t *tag_id;     /**< quad AMR Tag */
    sc_array_t *type_id;    /**< quad Type  */
    sc_array_t *quad_id;    /**< quad Number */
    sc_array_t *iblank;     /**< quad iblank */
}
vtk_plot_fluid_t; /**< Visualization Data Type: conservative fluid variables */

typedef struct {
    int ncell_scalars;
    int ncell_vectors;
    int npt_scalars;
    int npt_vectors;
    char **cell_names;
    char **pt_names;
    sc_array_t **cell_scalar_data;
    sc_array_t **cell_vector_data;
    sc_array_t **pt_scalar_data;
    sc_array_t **pt_vector_data;
}
vtk_data_t; /**< Visualization Data Type: conservative fluid variables */

#endif /* AMR_VAR_VTK_H */