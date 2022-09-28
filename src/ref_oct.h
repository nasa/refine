
/* Copyright 2006, 2014, 2021 United States Government as represented
 * by the Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine version 3 unstructured grid adaptation platform is
 * licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * https://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#ifndef REF_OCT_H
#define REF_OCT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_OCT_STRUCT REF_OCT_STRUCT;
typedef REF_OCT_STRUCT *REF_OCT;
END_C_DECLORATION

#include "ref_node.h"

BEGIN_C_DECLORATION
struct REF_OCT_STRUCT {
  REF_DBL bbox[6];
  REF_INT n, max, nnode;
  REF_INT *children;
  REF_INT *nodes;
};

#define ref_oct_n(ref_oct) ((ref_oct)->n)
#define ref_oct_max(ref_oct) ((ref_oct)->max)
#define ref_oct_nnode(ref_oct) ((ref_oct)->nnode)
#define ref_oct_child(ref_oct, corner, tree_node) \
  ((ref_oct)->children[(corner) + 8 * (tree_node)])
#define ref_oct_internal_node(ref_oct, tree_node) \
  (REF_EMPTY != ref_oct_child(ref_oct, 0, tree_node))
#define ref_oct_leaf_node(ref_oct, tree_node) \
  (REF_EMPTY == ref_oct_child(ref_oct, 0, tree_node))
#define ref_oct_c2n(ref_oct, corner, tree_node) \
  ((ref_oct)->nodes[(corner) + 27 * (tree_node)])

REF_FCN REF_STATUS ref_oct_create(REF_OCT *ref_oct);

REF_FCN REF_STATUS ref_oct_free(REF_OCT ref_oct);

REF_FCN REF_STATUS ref_oct_child_bbox(REF_DBL *parent_bbox, REF_INT child_index,
                                      REF_DBL *child_bbox);
REF_FCN REF_STATUS ref_oct_bbox_diag(REF_DBL *bbox, REF_DBL *diag);

REF_FCN REF_STATUS ref_oct_split(REF_OCT ref_oct, REF_INT node);
REF_FCN REF_STATUS ref_oct_split_at(REF_OCT ref_oct, REF_DBL *xyz, REF_DBL h);
REF_FCN REF_STATUS ref_oct_split_touching(REF_OCT ref_oct, REF_DBL *bbox,
                                          REF_DBL h);
REF_FCN REF_STATUS ref_oct_gradation(REF_OCT ref_oct);
REF_FCN REF_STATUS ref_oct_unique_nodes(REF_OCT ref_oct, REF_NODE ref_node);
REF_FCN REF_STATUS ref_oct_set_node_at(REF_OCT ref_oct, REF_INT insert_node,
                                       REF_DBL *xyz);

REF_FCN REF_STATUS ref_oct_contains(REF_OCT ref_oct, REF_DBL *xyz,
                                    REF_INT *node, REF_DBL *bbox);

REF_FCN REF_STATUS ref_oct_bbox_overlap(REF_DBL *bbox0, REF_DBL *bbox1,
                                        REF_BOOL *overlap);
REF_FCN REF_STATUS ref_oct_bbox_scale(REF_DBL *bbox0, REF_DBL factor,
                                      REF_DBL *bbox1);
REF_FCN REF_STATUS ref_oct_bbox_corner(REF_DBL *bbox, REF_INT corner,
                                       REF_DBL *xyz);

REF_FCN REF_STATUS ref_oct_tec(REF_OCT ref_oct, const char *filename);

REF_FCN REF_STATUS ref_oct_nleaf(REF_OCT ref_oct, REF_INT *nleaf);

END_C_DECLORATION

#endif /* REF_OCT_H */
