
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

#include "ref_oct.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_malloc.h"

REF_FCN REF_STATUS ref_oct_create(REF_OCT *ref_oct_ptr) {
  REF_OCT ref_oct;

  ref_malloc(*ref_oct_ptr, 1, REF_OCT_STRUCT);
  ref_oct = (*ref_oct_ptr);

  ref_oct->bbox[0] = 0.0;
  ref_oct->bbox[1] = 1.0;
  ref_oct->bbox[2] = 0.0;
  ref_oct->bbox[3] = 1.0;
  ref_oct->bbox[4] = 0.0;
  ref_oct->bbox[5] = 1.0;

  ref_oct->n = 1;
  ref_oct->max = 100 * 1024;
  ref_oct->nnode = 0;
  ref_oct->children = NULL;
  ref_malloc_init(ref_oct->children, 8 * ref_oct->max, REF_INT, REF_EMPTY);
  ref_malloc_init(ref_oct->nodes, 27 * ref_oct->max, REF_INT, REF_EMPTY);

  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_free(REF_OCT ref_oct) {
  if (NULL == (void *)ref_oct) return REF_NULL;
  ref_free(ref_oct->nodes);
  ref_free(ref_oct->children);
  ref_free(ref_oct);
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_child_bbox(REF_DBL *bbox, REF_INT child_index,
                                      REF_DBL *box) {
  switch (child_index) {
    case 0:
      box[0] = bbox[0];
      box[1] = 0.5 * (bbox[0] + bbox[1]);
      box[2] = bbox[2];
      box[3] = 0.5 * (bbox[2] + bbox[3]);
      box[4] = bbox[4];
      box[5] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 1:
      box[0] = 0.5 * (bbox[0] + bbox[1]);
      box[1] = bbox[1];
      box[2] = bbox[2];
      box[3] = 0.5 * (bbox[2] + bbox[3]);
      box[4] = bbox[4];
      box[5] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 2:
      box[0] = 0.5 * (bbox[0] + bbox[1]);
      box[1] = bbox[1];
      box[2] = 0.5 * (bbox[2] + bbox[3]);
      box[3] = bbox[3];
      box[4] = bbox[4];
      box[5] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 3:
      box[0] = bbox[0];
      box[1] = 0.5 * (bbox[0] + bbox[1]);
      box[2] = 0.5 * (bbox[2] + bbox[3]);
      box[3] = bbox[3];
      box[4] = bbox[4];
      box[5] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 4:
      box[0] = bbox[0];
      box[1] = 0.5 * (bbox[0] + bbox[1]);
      box[2] = bbox[2];
      box[3] = 0.5 * (bbox[2] + bbox[3]);
      box[4] = 0.5 * (bbox[4] + bbox[5]);
      box[5] = bbox[5];
      break;
    case 5:
      box[0] = 0.5 * (bbox[0] + bbox[1]);
      box[1] = bbox[1];
      box[2] = bbox[2];
      box[3] = 0.5 * (bbox[2] + bbox[3]);
      box[4] = 0.5 * (bbox[4] + bbox[5]);
      box[5] = bbox[5];
      break;
    case 6:
      box[0] = 0.5 * (bbox[0] + bbox[1]);
      box[1] = bbox[1];
      box[2] = 0.5 * (bbox[2] + bbox[3]);
      box[3] = bbox[3];
      box[4] = 0.5 * (bbox[4] + bbox[5]);
      box[5] = bbox[5];
      break;
    case 7:
      box[0] = bbox[0];
      box[1] = 0.5 * (bbox[0] + bbox[1]);
      box[2] = 0.5 * (bbox[2] + bbox[3]);
      box[3] = bbox[3];
      box[4] = 0.5 * (bbox[4] + bbox[5]);
      box[5] = bbox[5];
      break;
    default:
      THROW("not 2^3");
  }

  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_bbox_diag(REF_DBL *bbox, REF_DBL *diag) {
  *diag = sqrt(pow(bbox[1] - bbox[0], 2) + pow(bbox[3] - bbox[2], 2) +
               pow(bbox[5] - bbox[4], 2));
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_tattle(REF_OCT ref_oct, REF_INT node) {
  printf("corners %d %d %d %d  %d %d %d %d\n", ref_oct_c2n(ref_oct, 0, node),
         ref_oct_c2n(ref_oct, 1, node), ref_oct_c2n(ref_oct, 2, node),
         ref_oct_c2n(ref_oct, 3, node), ref_oct_c2n(ref_oct, 4, node),
         ref_oct_c2n(ref_oct, 5, node), ref_oct_c2n(ref_oct, 6, node),
         ref_oct_c2n(ref_oct, 7, node));
  printf("edges %d %d %d %d  %d %d %d %d  %d %d %d %d\n",
         ref_oct_c2n(ref_oct, 8, node), ref_oct_c2n(ref_oct, 9, node),
         ref_oct_c2n(ref_oct, 10, node), ref_oct_c2n(ref_oct, 11, node),
         ref_oct_c2n(ref_oct, 12, node), ref_oct_c2n(ref_oct, 13, node),
         ref_oct_c2n(ref_oct, 14, node), ref_oct_c2n(ref_oct, 15, node),
         ref_oct_c2n(ref_oct, 16, node), ref_oct_c2n(ref_oct, 17, node),
         ref_oct_c2n(ref_oct, 18, node), ref_oct_c2n(ref_oct, 19, node));
  printf("faces %d  %d %d %d %d  %d center %d\n",
         ref_oct_c2n(ref_oct, 20, node), ref_oct_c2n(ref_oct, 21, node),
         ref_oct_c2n(ref_oct, 22, node), ref_oct_c2n(ref_oct, 23, node),
         ref_oct_c2n(ref_oct, 24, node), ref_oct_c2n(ref_oct, 25, node),
         ref_oct_c2n(ref_oct, 26, node));
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_split(REF_OCT ref_oct, REF_INT node) {
  REF_INT i;
  RAS(node >= 0, "node negative");
  RAS(node < ref_oct_n(ref_oct), "node larger than n");
  RAS(node < ref_oct_max(ref_oct), "node larger than max");
  for (i = 0; i < 8; i++) {
    REIS(REF_EMPTY, ref_oct_child(ref_oct, i, node), "child not empty");
  }
  if (ref_oct_n(ref_oct) + 8 > ref_oct_max(ref_oct)) {
    REF_INT orig, chunk;
    REF_INT max_limit = REF_INT_MAX / 27;
    orig = ref_oct_max(ref_oct);
    /* geometric growth for efficiency */
    chunk = (REF_INT)(1.5 * (REF_DBL)orig);
    RAS(max_limit - orig > 0, "chunk limit at max");
    chunk = MIN(chunk, max_limit - orig);
    ref_oct_max(ref_oct) = orig + chunk;
    ref_realloc(ref_oct->children, 8 * ref_oct->max, REF_INT);
    ref_realloc(ref_oct->nodes, 27 * ref_oct->max, REF_INT);
    for (i = 8 * orig; i < 8 * ref_oct_max(ref_oct); i++)
      ref_oct->children[i] = REF_EMPTY;
    for (i = 27 * orig; i < 27 * ref_oct_max(ref_oct); i++)
      ref_oct->nodes[i] = REF_EMPTY;
  }
  for (i = 0; i < 8; i++) {
    ref_oct_child(ref_oct, i, node) = ref_oct_n(ref_oct);
    (ref_oct_n(ref_oct))++;
  }
  return REF_SUCCESS;
}

REF_STATUS ref_oct_split_at(REF_OCT ref_oct, REF_DBL *xyz, REF_DBL h) {
  REF_INT node;
  REF_DBL bbox[6], diag;
  RSS(ref_oct_contains(ref_oct, xyz, &node, bbox), "contains oct");
  RAS(node >= 0, "not found");
  RSS(ref_oct_bbox_diag(bbox, &diag), "bbox diag");
  if (diag > h) {
    RSS(ref_oct_split(ref_oct, node), "split");
    RSS(ref_oct_split_at(ref_oct, xyz, h), "again");
  }
  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_split_touching_node(
    REF_OCT ref_oct, REF_INT node, REF_DBL *my_bbox, REF_DBL *bbox, REF_DBL h) {
  REF_BOOL overlap;
  REF_DBL diag;
  REF_INT child_index;
  RSS(ref_oct_bbox_overlap(my_bbox, bbox, &overlap), "overlap");
  if (!overlap) return REF_SUCCESS;
  if (ref_oct_leaf_node(ref_oct, node)) {
    RSS(ref_oct_bbox_diag(my_bbox, &diag), "bbox diag");
    if (diag > h) {
      RSS(ref_oct_split(ref_oct, node), "split");
    } else {
      return REF_SUCCESS;
    }
  }
  for (child_index = 0; child_index < 8; child_index++) {
    REF_DBL box[6];
    RSS(ref_oct_child_bbox(my_bbox, child_index, box), "bbox");
    RSS(ref_oct_split_touching_node(
            ref_oct, ref_oct->children[child_index + 8 * node], box, bbox, h),
        "recurse");
  }
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_split_touching(REF_OCT ref_oct, REF_DBL *bbox,
                                          REF_DBL h) {
  RSS(ref_oct_split_touching_node(ref_oct, 0, ref_oct->bbox, bbox, h),
      "descend");
  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_gradation_node(REF_OCT ref_oct, REF_INT node,
                                                 REF_DBL *bbox) {
  if (ref_oct_leaf_node(ref_oct, node)) {
    REF_DBL tool[6], factor = 1.1, diag, h;
    RSS(ref_oct_bbox_scale(bbox, factor, tool), "scale");
    RSS(ref_oct_bbox_diag(bbox, &diag), "scale");
    h = factor * 2.0 * diag;
    RSS(ref_oct_split_touching(ref_oct, tool, h), "split region");
  } else {
    REF_INT child_index;
    for (child_index = 0; child_index < 8; child_index++) {
      REF_DBL box[6];
      RSS(ref_oct_child_bbox(bbox, child_index, box), "bbox");
      RSS(ref_oct_gradation_node(
              ref_oct, ref_oct_child(ref_oct, child_index, node), box),
          "recurse");
    }
  }
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_gradation(REF_OCT ref_oct) {
  REF_INT n, last_n, nleaf;
  last_n = REF_EMPTY;
  n = ref_oct_n(ref_oct);
  while (n != last_n) {
    RSS(ref_oct_nleaf(ref_oct, &nleaf), "count leaves");
    printf("ncell %d nleaf %d %5.1f %%\n", n, nleaf,
           ((REF_DBL)nleaf / (REF_DBL)n) * 100.0);
    RSS(ref_oct_gradation_node(ref_oct, 0, ref_oct->bbox), "descend");
    last_n = n;
    n = ref_oct_n(ref_oct);
  }
  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_unique_corner_nodes_node(REF_OCT ref_oct,
                                                           REF_INT node,
                                                           REF_DBL *bbox,
                                                           REF_NODE ref_node) {
  if (ref_oct_leaf_node(ref_oct, node)) {
    REF_INT corner;
    for (corner = 0; corner < 8; corner++) {
      if (REF_EMPTY == ref_oct_c2n(ref_oct, corner, node)) {
        REF_INT new_node;
        REF_DBL xyz[3];
        REF_INT insert_node;
        insert_node = ref_oct_nnode(ref_oct);
        ref_oct_nnode(ref_oct)++;
        RSS(ref_node_add(ref_node, insert_node, &new_node), "add node");
        REIS(insert_node, new_node, "expects to match");
        RSS(ref_oct_bbox_corner(bbox, corner, xyz), "corner xyz");
        ref_node_xyz(ref_node, 0, new_node) = xyz[0];
        ref_node_xyz(ref_node, 1, new_node) = xyz[1];
        ref_node_xyz(ref_node, 2, new_node) = xyz[2];
        RSS(ref_oct_set_node_at(ref_oct, insert_node, xyz), "set node");
      }
      RUS(REF_EMPTY, ref_oct_c2n(ref_oct, corner, node), "self not set");
    }
  } else {
    REF_INT child_index;
    for (child_index = 0; child_index < 8; child_index++) {
      REF_DBL box[6];
      RSS(ref_oct_child_bbox(bbox, child_index, box), "bbox");
      RSS(ref_oct_unique_corner_nodes_node(
              ref_oct, ref_oct_child(ref_oct, child_index, node), box,
              ref_node),
          "recurse");
    }
  }
  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_unique_face_nodes_node(REF_OCT ref_oct,
                                                         REF_INT node,
                                                         REF_DBL *bbox,
                                                         REF_NODE ref_node) {
  if (ref_oct_leaf_node(ref_oct, node)) {
    REF_INT cell_face;
    for (cell_face = 0; cell_face < 6; cell_face++) {
      REF_DBL xyz[3], cxyz[3];
      REF_INT insert_node, i, j;
      REF_INT nodes[9];
      REF_BOOL has_edge;
      RSS(ref_oct_he2_qu2(cell_face, nodes), "he2qu2");
      has_edge = REF_EMPTY != ref_oct_c2n(ref_oct, nodes[4], node) ||
                 REF_EMPTY != ref_oct_c2n(ref_oct, nodes[5], node) ||
                 REF_EMPTY != ref_oct_c2n(ref_oct, nodes[6], node) ||
                 REF_EMPTY != ref_oct_c2n(ref_oct, nodes[7], node);
      if (has_edge && REF_EMPTY == ref_oct_c2n(ref_oct, nodes[8], node)) {
        REF_INT new_node;
        insert_node = ref_oct_nnode(ref_oct);
        ref_oct_nnode(ref_oct)++;
        RSS(ref_node_add(ref_node, insert_node, &new_node), "add node");
        REIS(insert_node, new_node, "expects to match");
        for (i = 0; i < 3; i++) xyz[i] = 0.0;
        for (j = 0; j < 4; j++) {
          RSS(ref_oct_bbox_corner(bbox, nodes[0], cxyz), "corner xyz");
          for (i = 0; i < 3; i++) xyz[i] += 0.25 * cxyz[i];
        }
        ref_node_xyz(ref_node, 0, new_node) = xyz[0];
        ref_node_xyz(ref_node, 1, new_node) = xyz[1];
        ref_node_xyz(ref_node, 2, new_node) = xyz[2];
        RSS(ref_oct_set_node_at(ref_oct, insert_node, xyz), "set node");
        RUS(REF_EMPTY, ref_oct_c2n(ref_oct, nodes[8], node), "self not set");
      }
    }
  } else {
    REF_INT child_index;
    for (child_index = 0; child_index < 8; child_index++) {
      REF_DBL box[6];
      RSS(ref_oct_child_bbox(bbox, child_index, box), "bbox");
      RSS(ref_oct_unique_face_nodes_node(
              ref_oct, ref_oct_child(ref_oct, child_index, node), box,
              ref_node),
          "recurse");
    }
  }
  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_unique_center_nodes_node(REF_OCT ref_oct,
                                                           REF_INT node,
                                                           REF_DBL *bbox,
                                                           REF_NODE ref_node) {
  if (ref_oct_leaf_node(ref_oct, node)) {
    REF_BOOL has_edge;
    REF_INT edge_node;
    has_edge = REF_FALSE;
    for (edge_node = 8; edge_node < 20; edge_node++)
      has_edge =
          has_edge || (REF_EMPTY != ref_oct_c2n(ref_oct, edge_node, node));
    if (has_edge) {
      REF_DBL xyz[3], cxyz[3];
      REF_INT insert_node, i, j;
      REF_INT new_node;
      insert_node = ref_oct_nnode(ref_oct);
      ref_oct_nnode(ref_oct)++;
      RSS(ref_node_add(ref_node, insert_node, &new_node), "add node");
      REIS(insert_node, new_node, "expects to match");
      for (i = 0; i < 3; i++) xyz[i] = 0.0;
      for (j = 0; j < 8; j++) {
        RSS(ref_oct_bbox_corner(bbox, j, cxyz), "corner xyz");
        for (i = 0; i < 3; i++) xyz[i] += 0.125 * cxyz[i];
      }
      ref_node_xyz(ref_node, 0, new_node) = xyz[0];
      ref_node_xyz(ref_node, 1, new_node) = xyz[1];
      ref_node_xyz(ref_node, 2, new_node) = xyz[2];
      RSS(ref_oct_set_node_at(ref_oct, insert_node, xyz), "set node");
      RUS(REF_EMPTY, ref_oct_c2n(ref_oct, 26, node), "self not set");
    }
  } else {
    REF_INT child_index;
    for (child_index = 0; child_index < 8; child_index++) {
      REF_DBL box[6];
      RSS(ref_oct_child_bbox(bbox, child_index, box), "bbox");
      RSS(ref_oct_unique_face_nodes_node(
              ref_oct, ref_oct_child(ref_oct, child_index, node), box,
              ref_node),
          "recurse");
    }
  }
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_unique_nodes(REF_OCT ref_oct, REF_NODE ref_node) {
  REIS(0, ref_oct_nnode(ref_oct), "expected zero oct nodes");
  REIS(0, ref_node_n(ref_node), "expected zero grid nodes");
  RSS(ref_oct_unique_corner_nodes_node(ref_oct, 0, ref_oct->bbox, ref_node),
      "corner node descend");
  RSS(ref_oct_unique_face_nodes_node(ref_oct, 0, ref_oct->bbox, ref_node),
      "face node descend");
  RSS(ref_oct_unique_center_nodes_node(ref_oct, 0, ref_oct->bbox, ref_node),
      "center node descend");
  RSS(ref_node_initialize_n_global(ref_node, ref_oct_nnode(ref_oct)),
      "init glob");
  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_set_node_at_node(REF_OCT ref_oct,
                                                   REF_INT node, REF_DBL *bbox,
                                                   REF_INT insert_node,
                                                   REF_DBL *xyz) {
  REF_DBL fuzz[6], point[6], expand = 1.01;
  REF_BOOL overlap;

  point[0] = xyz[0];
  point[1] = xyz[0];
  point[2] = xyz[1];
  point[3] = xyz[1];
  point[4] = xyz[2];
  point[5] = xyz[2];
  RSS(ref_oct_bbox_scale(bbox, expand, fuzz), "scale bbox");
  RSS(ref_oct_bbox_overlap(point, fuzz, &overlap), "overlap");
  if (!overlap) return REF_SUCCESS;

  if (ref_oct_leaf_node(ref_oct, node)) {
    REF_INT corner;
    REF_DBL h;
    RSS(ref_oct_bbox_diag(bbox, &h), "diag");
    for (corner = 0; corner < 27; corner++) {
      REF_DBL my_xyz[3], dist, tol = 0.1;
      RSS(ref_oct_bbox_corner(bbox, corner, my_xyz), "corner xyz");
      dist = sqrt(pow(xyz[0] - my_xyz[0], 2) + pow(xyz[1] - my_xyz[1], 2) +
                  pow(xyz[2] - my_xyz[2], 2));
      if (dist < tol * h) {
        ref_oct_c2n(ref_oct, corner, node) = insert_node;
      }
    }
  } else {
    REF_INT child_index;
    for (child_index = 0; child_index < 8; child_index++) {
      REF_DBL box[6];
      RSS(ref_oct_child_bbox(bbox, child_index, box), "bbox");
      RUS(node, ref_oct_child(ref_oct, child_index, node), "same");
      RSS(ref_oct_set_node_at_node(ref_oct,
                                   ref_oct_child(ref_oct, child_index, node),
                                   box, insert_node, xyz),
          "recurse");
    }
  }
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_set_node_at(REF_OCT ref_oct, REF_INT insert_node,
                                       REF_DBL *xyz) {
  RSS(ref_oct_set_node_at_node(ref_oct, 0, ref_oct->bbox, insert_node, xyz),
      "descend");
  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_contains_node(REF_OCT ref_oct, REF_DBL *xyz,
                                                REF_DBL *bbox, REF_INT current,
                                                REF_INT *node,
                                                REF_DBL *node_bbox) {
  REF_INT child_index;
  *node = REF_EMPTY;
  if (xyz[0] < bbox[0] || bbox[1] < xyz[0] || xyz[1] < bbox[2] ||
      bbox[3] < xyz[1] || xyz[2] < bbox[4] || bbox[5] < xyz[2]) {
    *node = REF_EMPTY;
    return REF_SUCCESS;
  }
  if (ref_oct_leaf_node(ref_oct, current)) {
    REF_INT i;
    for (i = 0; i < 6; i++) node_bbox[i] = bbox[i];
    *node = current;
    return REF_SUCCESS;
  }
  for (child_index = 0; child_index < 8; child_index++) {
    REF_DBL box[6];
    RSS(ref_oct_child_bbox(bbox, child_index, box), "bbox");
    RSS(ref_oct_contains_node(ref_oct, xyz, box,
                              ref_oct->children[child_index + 8 * current],
                              node, node_bbox),
        "recurse");
    if (*node != REF_EMPTY) {
      return REF_SUCCESS;
    }
  }
  return REF_SUCCESS;
}
REF_FCN REF_STATUS ref_oct_contains(REF_OCT ref_oct, REF_DBL *xyz,
                                    REF_INT *node, REF_DBL *bbox) {
  RSS(ref_oct_contains_node(ref_oct, xyz, ref_oct->bbox, 0, node, bbox),
      "wrapper");
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_bbox_corner(REF_DBL *bbox, REF_INT corner,
                                       REF_DBL *xyz) {
  switch (corner) {
    case 0:
      xyz[0] = bbox[0];
      xyz[1] = bbox[2];
      xyz[2] = bbox[4];
      break;
    case 1:
      xyz[0] = bbox[1];
      xyz[1] = bbox[2];
      xyz[2] = bbox[4];
      break;
    case 2:
      xyz[0] = bbox[1];
      xyz[1] = bbox[3];
      xyz[2] = bbox[4];
      break;
    case 3:
      xyz[0] = bbox[0];
      xyz[1] = bbox[3];
      xyz[2] = bbox[4];
      break;
    case 4:
      xyz[0] = bbox[0];
      xyz[1] = bbox[2];
      xyz[2] = bbox[5];
      break;
    case 5:
      xyz[0] = bbox[1];
      xyz[1] = bbox[2];
      xyz[2] = bbox[5];
      break;
    case 6:
      xyz[0] = bbox[1];
      xyz[1] = bbox[3];
      xyz[2] = bbox[5];
      break;
    case 7:
      xyz[0] = bbox[0];
      xyz[1] = bbox[3];
      xyz[2] = bbox[5];
      break;

    case 8:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = bbox[2];
      xyz[2] = bbox[4];
      break;
    case 9:
      xyz[0] = bbox[1];
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = bbox[4];
      break;
    case 10:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = bbox[3];
      xyz[2] = bbox[4];
      break;
    case 11:
      xyz[0] = bbox[0];
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = bbox[4];
      break;

    case 12:
      xyz[0] = bbox[0];
      xyz[1] = bbox[2];
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 13:
      xyz[0] = bbox[1];
      xyz[1] = bbox[2];
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 14:
      xyz[0] = bbox[1];
      xyz[1] = bbox[3];
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 15:
      xyz[0] = bbox[0];
      xyz[1] = bbox[3];
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;

    case 16:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = bbox[2];
      xyz[2] = bbox[5];
      break;
    case 17:
      xyz[0] = bbox[1];
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = bbox[5];
      break;
    case 18:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = bbox[3];
      xyz[2] = bbox[5];
      break;
    case 19:
      xyz[0] = bbox[0];
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = bbox[5];
      break;

    case 20:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = bbox[4];
      break;

    case 21:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = bbox[2];
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 22:
      xyz[0] = bbox[1];
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 23:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = bbox[3];
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;
    case 24:
      xyz[0] = bbox[0];
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;

    case 25:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = bbox[5];
      break;

    case 26:
      xyz[0] = 0.5 * (bbox[0] + bbox[1]);
      xyz[1] = 0.5 * (bbox[2] + bbox[3]);
      xyz[2] = 0.5 * (bbox[4] + bbox[5]);
      break;

    default:
      THROW("not Q2 index");
  }
  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_tec_node(REF_OCT ref_oct, REF_INT node,
                                           REF_DBL *bbox, FILE *f) {
  REF_DBL box[6];
  REF_INT i;

  RAS(0 <= node && node < ref_oct->n, "out of range node");
  if (ref_oct_leaf_node(ref_oct, node)) {
    for (i = 0; i < 8; i++) {
      REF_DBL xyz[3];
      RSS(ref_oct_bbox_corner(bbox, i, xyz), "corner xyz");
      fprintf(f, "%f %f %f\n", xyz[0], xyz[1], xyz[2]);
    }
  } else {
    for (i = 0; i < 8; i++) {
      RSS(ref_oct_child_bbox(bbox, i, box), "bbox");
      if (ref_oct_internal_node(ref_oct, node))
        RSS(ref_oct_tec_node(ref_oct, ref_oct_child(ref_oct, i, node), box, f),
            "recurse");
    }
  }

  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_bbox_overlap(REF_DBL *bbox0, REF_DBL *bbox1,
                                        REF_BOOL *overlap) {
  *overlap = REF_FALSE;
  if (bbox0[1] < bbox1[0] || bbox1[1] < bbox0[0]) return REF_SUCCESS;
  if (bbox0[3] < bbox1[2] || bbox1[3] < bbox0[2]) return REF_SUCCESS;
  if (bbox0[5] < bbox1[4] || bbox1[5] < bbox0[4]) return REF_SUCCESS;
  *overlap = REF_TRUE;
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_bbox_scale(REF_DBL *bbox0, REF_DBL factor,
                                      REF_DBL *bbox1) {
  bbox1[0] = bbox0[0] - 0.5 * (factor - 1.0) * (bbox0[1] - bbox0[0]);
  bbox1[1] = bbox0[1] + 0.5 * (factor - 1.0) * (bbox0[1] - bbox0[0]);
  bbox1[2] = bbox0[2] - 0.5 * (factor - 1.0) * (bbox0[3] - bbox0[2]);
  bbox1[3] = bbox0[3] + 0.5 * (factor - 1.0) * (bbox0[3] - bbox0[2]);
  bbox1[4] = bbox0[4] - 0.5 * (factor - 1.0) * (bbox0[5] - bbox0[4]);
  bbox1[5] = bbox0[5] + 0.5 * (factor - 1.0) * (bbox0[5] - bbox0[4]);
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_tec(REF_OCT ref_oct, const char *filename) {
  FILE *f;
  REF_INT i, nleaf;
  const char *zonetype = "febrick";
  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  RSS(ref_oct_nleaf(ref_oct, &nleaf), "count leaves");

  fprintf(f, "title=\"tecplot refine octree\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  fprintf(
      f,
      "zone t=\"octree\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      8 * nleaf, nleaf, "point", zonetype);

  RSS(ref_oct_tec_node(ref_oct, 0, ref_oct->bbox, f), "draw root box");

  for (i = 0; i < nleaf; i++) {
    fprintf(f, "%d %d %d %d %d %d %d %d\n", 1 + 8 * i, 2 + 8 * i, 3 + 8 * i,
            4 + 8 * i, 5 + 8 * i, 6 + 8 * i, 7 + 8 * i, 8 + 8 * i);
  }
  fclose(f);
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_nleaf(REF_OCT ref_oct, REF_INT *nleaf) {
  REF_INT i;
  *nleaf = 0;
  for (i = 0; i < ref_oct_n(ref_oct); i++) {
    if (ref_oct_leaf_node(ref_oct, i)) (*nleaf)++;
  }
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_he2_qu2(REF_INT cell_face, REF_INT *qu2) {
  switch (cell_face) {
    case 0:
      qu2[0] = 0;
      qu2[1] = 3;
      qu2[2] = 7;
      qu2[3] = 4;
      qu2[4] = 11;
      qu2[5] = 15;
      qu2[6] = 19;
      qu2[7] = 12;
      qu2[8] = 24;
      break;
    case 1:
      qu2[0] = 1;
      qu2[1] = 5;
      qu2[2] = 6;
      qu2[3] = 2;
      qu2[4] = 13;
      qu2[5] = 17;
      qu2[6] = 14;
      qu2[7] = 9;
      qu2[8] = 22;
      break;
    case 2:
      qu2[0] = 0;
      qu2[1] = 4;
      qu2[2] = 5;
      qu2[3] = 1;
      qu2[4] = 12;
      qu2[5] = 16;
      qu2[6] = 13;
      qu2[7] = 8;
      qu2[8] = 21;
      break;
    case 3:
      qu2[0] = 2;
      qu2[1] = 6;
      qu2[2] = 7;
      qu2[3] = 3;
      qu2[4] = 14;
      qu2[5] = 18;
      qu2[6] = 15;
      qu2[7] = 10;
      qu2[8] = 23;
      break;
    case 4:
      qu2[0] = 0;
      qu2[1] = 1;
      qu2[2] = 2;
      qu2[3] = 3;
      qu2[4] = 8;
      qu2[5] = 9;
      qu2[6] = 10;
      qu2[7] = 11;
      qu2[8] = 20;
      break;
    case 5:
      qu2[0] = 4;
      qu2[1] = 7;
      qu2[2] = 6;
      qu2[3] = 5;
      qu2[4] = 19;
      qu2[5] = 18;
      qu2[6] = 17;
      qu2[7] = 16;
      qu2[8] = 25;
      break;
    default:
      THROW("not 0-5 cell_face");
  }

  return REF_SUCCESS;
}

REF_FCN static REF_STATUS ref_oct_export_node(REF_OCT ref_oct, REF_INT node,
                                              REF_DBL *bbox,
                                              REF_GRID ref_grid) {
  REF_INT i;
  if (ref_oct_leaf_node(ref_oct, node)) {
    REF_CELL ref_cell;
    REF_INT nodes[REF_CELL_MAX_SIZE_PER], new_cell;
    REF_DBL diag, tol = 0.1;
    REF_INT cell_face, face_node;
    for (i = 0; i < 8; i++)
      RAS(REF_EMPTY != ref_oct_c2n(ref_oct, i, node), "expects to be set");
    for (i = 8; i < 27; i++)
      REIB(REF_EMPTY, ref_oct_c2n(ref_oct, i, node), "implement 2-1",
           { ref_oct_tattle(ref_oct, node); });
    for (i = 0; i < 8; i++) nodes[i] = ref_oct_c2n(ref_oct, i, node);
    ref_cell = ref_grid_hex(ref_grid);
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add hex");
    ref_cell = ref_grid_qua(ref_grid);
    RSS(ref_oct_bbox_diag(bbox, &diag), "diag");

    for (cell_face = 0; cell_face < 6; cell_face++) {
      if (ABS(bbox[cell_face] - ref_oct->bbox[cell_face]) < tol * diag) {
        for (face_node = 0; face_node < 4; face_node++) {
          nodes[face_node] =
              ref_cell_f2n_gen(ref_grid_hex(ref_grid), face_node, cell_face);
          nodes[face_node] = ref_oct_c2n(ref_oct, nodes[face_node], node);
        }
        nodes[ref_cell_id_index(ref_cell)] = cell_face + 1;
        RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add quad");
      }
    }

  } else {
    REF_INT child_index;
    for (i = 0; i < 27; i++)
      REIS(REF_EMPTY, ref_oct_c2n(ref_oct, i, node),
           "interor node with grid node");
    for (child_index = 0; child_index < 8; child_index++) {
      REF_DBL box[6];
      RSS(ref_oct_child_bbox(bbox, child_index, box), "bbox");
      RSS(ref_oct_export_node(ref_oct,
                              ref_oct_child(ref_oct, child_index, node), box,
                              ref_grid),
          "recurse");
    }
  }
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_export(REF_OCT ref_oct, REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  RSS(ref_oct_unique_nodes(ref_oct, ref_node), "make nodes");
  RSS(ref_oct_export_node(ref_oct, 0, ref_oct->bbox, ref_grid), "descend");
  return REF_SUCCESS;
}
