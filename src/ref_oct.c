
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
  ref_oct->children = NULL;
  ref_malloc_init(ref_oct->children, 8 * ref_oct->max, REF_INT, REF_EMPTY);

  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_free(REF_OCT ref_oct) {
  if (NULL == (void *)ref_oct) return REF_NULL;
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

REF_FCN REF_STATUS ref_oct_split(REF_OCT ref_oct, REF_INT node) {
  REF_INT i;
  RAS(node >= 0, "node negative");
  RAS(node < ref_oct_n(ref_oct), "node larger than n");
  RAS(node < ref_oct_max(ref_oct), "node larger than max");
  for (i = 0; i < 8; i++) {
    REIS(REF_EMPTY, ref_oct_child(ref_oct, i, node), "child not empty");
  }
  if (ref_oct_n(ref_oct) + 8 > ref_oct_max(ref_oct)) {
    THROW("out of children, implememnt realloc");
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
              ref_oct, ref_oct->children[child_index + 8 * node], box),
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
    printf("ncell %d nleaf %d\n", n, nleaf);
    RSS(ref_oct_gradation_node(ref_oct, 0, ref_oct->bbox), "descend");
    last_n = n;
    n = ref_oct_n(ref_oct);
  }
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

REF_FCN static REF_STATUS ref_oct_tec_node(REF_OCT ref_oct, REF_INT node,
                                           REF_DBL *bbox, FILE *f) {
  REF_DBL box[6];
  REF_INT i;

  RAS(0 <= node && node < ref_oct->n, "out of range node");
  fprintf(f, "%f %f %f\n", bbox[0], bbox[2], bbox[4]);
  fprintf(f, "%f %f %f\n", bbox[1], bbox[2], bbox[4]);
  fprintf(f, "%f %f %f\n", bbox[1], bbox[3], bbox[4]);
  fprintf(f, "%f %f %f\n", bbox[0], bbox[3], bbox[4]);
  fprintf(f, "%f %f %f\n", bbox[0], bbox[2], bbox[5]);
  fprintf(f, "%f %f %f\n", bbox[1], bbox[2], bbox[5]);
  fprintf(f, "%f %f %f\n", bbox[1], bbox[3], bbox[5]);
  fprintf(f, "%f %f %f\n", bbox[0], bbox[3], bbox[5]);

  for (i = 0; i < 8; i++) {
    RSS(ref_oct_child_bbox(bbox, i, box), "bbox");
    if (ref_oct_internal_node(ref_oct, node))
      RSS(ref_oct_tec_node(ref_oct, ref_oct_child(ref_oct, i, node), box, f),
          "recurse");
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
  REF_INT i;
  const char *zonetype = "febrick";
  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "title=\"tecplot refine octree\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  fprintf(f,
          "zone t=\"octree\", nodes=%d, elements=%d, datapacking=%s, "
          "zonetype=%s\n",
          8 * ref_oct->n, ref_oct->n, "point", zonetype);

  RSS(ref_oct_tec_node(ref_oct, 0, ref_oct->bbox, f), "draw root box");

  for (i = 0; i < ref_oct_n(ref_oct); i++) {
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
