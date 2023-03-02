/**
 * \file    amr_var_tag.h
 * \ingroup amr_group
 * \author  akirby
 *
 * \brief   AMR enum types
 */

#ifndef AMR_VAR_TAG_H
#define AMR_VAR_TAG_H

/**
 * Name assignment for different amr tags
 */
typedef enum {
    no_tag,         /**< No tagging assigned */
    point_tag,      /**< Tagged from point refinement callback */
    feature_tag,    /**< Tagged from feature refinement callback */
    spread_tag,     /**< Tagged from spreading callback */
    level_tag,      /**< Tagged from min level refinement callback */
    buffer_tag      /**< Tagged from buffer refinement callback */
}
tag_names;

#endif /* AMR_VAR_TAG_H */