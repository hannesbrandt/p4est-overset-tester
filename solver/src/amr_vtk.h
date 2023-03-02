/**
 * \file   amr_vtk.h
 * \author akirby
 */

#ifndef AMR_VTK_H
#define AMR_VTK_H

/* header files */
#include "amr_var_defines.h"
#include "amr_var_quad.h"
#include "amr_var_vtk.h"

#include "amr_p4est_utilities.h"
#include "amr_utilities.h"
#include "dg4est_solver.hxx"
#include "precision_types.h"

#ifndef DOXYGEN_IGNORE
#  include <math.h>
#endif

#ifndef P4_TO_P8
#  include <p4est_vtk.h>
#  include <p4est_extended.h>
#else
#  include <p8est_vtk.h>
#  include <p8est_extended.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Visualize function for just the AMR grid
 *
 * @param [in] p4est        p4est forest structure
 */
void vtk_visualize_grid(p4est_t *p4est);

/** Visualize function for the AMR grid and solution
 *
 * @param [in] filepath     path to output directory
 * @param [in] p4est        p4est forest structure
 * @param [in] timestep     output time step number
 * @param [in] HO_VIZ       flag for VTK high-order visualization (paraview)
 * @param [in] HO_NPTS      order of VTK visualization (max 10)
 * @return visualization wall-clock time
 */
Real vtk_write_solution(char *filepath,p4est_t *p4est,int timestep,
                        int HO_VIZ,int HO_NPTS);

/** Interpolation function of the quadrant data to the corners (external)
 *
 * @param [in] info             p4est volume information related to quadrant
 * @param [in,out] user_data    interpolated data
 */
void interpolate_solution(p4est_iter_volume_info_t *info, void *user_data);

static void
GLL_points(int npts,Real *points){
    switch (npts) {
        case(1):
            points[0]  = 0.000000000000000000;
            break;
        case(2):
            points[0] = -1.000000000000000000;
            points[1] =  1.000000000000000000;
            break;
        case(3):
            points[0] = -1.000000000000000000;
            points[1] =  0.000000000000000000;
            points[2] =  1.000000000000000000;
            break;
        case(4):
            points[0] = -1.000000000000000000;
            points[1] = -0.447213595499957939;
            points[2] =  0.447213595499957939;
            points[3] =  1.000000000000000000;
            break;
        case(5):
            points[0] = -1.000000000000000000;
            points[1] = -0.654653670707977144;
            points[2] =  0.000000000000000000;
            points[3] =  0.654653670707977144;
            points[4] =  1.000000000000000000;
            break;
        case(6):
            points[0] = -1.000000000000000000;
            points[1] = -0.765055323929464690;
            points[2] = -0.285231516480645096;
            points[3] =  0.285231516480645096;
            points[4] =  0.765055323929464690;
            points[5] =  1.000000000000000000;
            break;
        case(7):
            points[0] = -1.000000000000000000;
            points[1] = -0.83022389627856693;
            points[2] = -0.468848793470714214;
            points[3] =  0.000000000000000000;
            points[4] =  0.468848793470714214;
            points[5] =  0.83022389627856693;
            points[6] =  1.000000000000000000;
            break;
        case(8):
            points[0] = -1.000000000000000000;
            points[1] = -0.871740148509606615;
            points[2] = -0.591700181433142302;
            points[3] = -0.209299217902478869;
            points[4] =  0.209299217902478869;
            points[5] =  0.591700181433142302;
            points[6] =  0.871740148509606615;
            points[7] =  1.000000000000000000;
            break;
        case(9):
            points[0] = -1.000000000000000000;
            points[1] = -0.899757995411460157;
            points[2] = -0.677186279510737753;
            points[3] = -0.363117463826178159;
            points[4] =  0.000000000000000000;
            points[5] =  0.363117463826178159;
            points[6] =  0.677186279510737753;
            points[7] =  0.899757995411460157;
            points[8] =  1.000000000000000000;
            break;
        case(10):
            points[0]  = -1.000000000000000000;
            points[1]  = -0.919533908166458814;
            points[2]  = -0.738773865105505075;
            points[3]  = -0.477924949810444496;
            points[4]  = -0.165278957666387025;
            points[5]  =  0.165278957666387025;
            points[6]  =  0.477924949810444496;
            points[7]  =  0.738773865105505075;
            points[8]  =  0.919533908166458814;
            points[9]  =  1.000000000000000000;
            break;
    }
}

#ifdef __cplusplus
}
#endif
#endif /* AMR_VTK_H */