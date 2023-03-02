/**
 * \file   amr_tag.h
 * \author akirby
 */

#ifndef AMR_TAG_H
#define AMR_TAG_H

/* header files */
#include "precision_types.h"
#include "dg4est_solver.hxx"
#include "amr_var_tag.h"
#include "amr_p4est_utilities.h"

#ifndef DOXYGEN_IGNORE
#  include <math.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Resets all refinement tags to no_tag
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] tag_default  tag number to set all quadrant tags
 * @param [in] keep_pt_tag  flag to keep point tags active
 * @return flag for refinement of quadrant (always 0)
 */
int tag_reset(p4est_t *p4est,int tag_default,char keep_pt_tag);

/** Resets all refinement tags on top level to tag_default
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] tag_default  tag number to set all quadrant tags
 * @return flag for refinement of quadrant (always 0)
 */
int tag_reset_max_level(p4est_t *p4est,int tag_default);

/** Tags quadrant for refinement if user-specified point is found in quadrant
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_point(p4est_t *p4est,
              p4est_topidx_t which_tree,
              p4est_quadrant_t *q);

/** Point search in p4est forest for list of points
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @param [in] local_num    leaf id number
 * @param [in] point        user-specified point list
 * @return flag for point found in quadrant
 */
int tag_point_search(p4est_t *p4est,
                     p4est_topidx_t which_tree,
                     p4est_quadrant_t *quadrant,
                     p4est_locidx_t local_num,
                     void *point);

/** Tags quadrant for refinement if user-specified feature found in quadrant
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_feature(p4est_t *p4est,
                p4est_topidx_t which_tree,
                p4est_quadrant_t *q);

/** Tags quadrant for refinement if spread tag is set in quadrant
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_spread(p4est_t *p4est,
               p4est_topidx_t which_tree,
               p4est_quadrant_t *q);

/** Tags quadrant for refinement if buffer tag is set in quadrant, e.g. adjacent hanging sides
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_buffer(p4est_t *p4est,
               p4est_topidx_t which_tree,
               p4est_quadrant_t *q);

/** Tags quadrant for coarsening if no refinement tags found
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] children     quadrant children tagged for coarsening
 * @return flag for coarsening quadrant
 */
int tag_coarsen(p4est_t *p4est,
                p4est_topidx_t which_tree,
                p4est_quadrant_t *children[]);

/** Tags quadrant for refinement if quadrant level < grid->min_level
 *
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_level(p4est_t *p4est,
              p4est_topidx_t which_tree,
              p4est_quadrant_t *q);

#ifdef __cplusplus
}
#endif
#endif /* AMR_TAG_H */