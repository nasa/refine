

/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#include "ref_geom.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_egads.h"
#include "ref_export.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_meshlink.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

REF_STATUS ref_geom_initialize(REF_GEOM ref_geom) {
  REF_INT geom;
  ref_geom_n(ref_geom) = 0;
  for (geom = 0; geom < ref_geom_max(ref_geom); geom++) {
    ref_geom_type(ref_geom, geom) = REF_EMPTY;
    ref_geom_id(ref_geom, geom) = geom + 1;
  }
  ref_geom_id(ref_geom, ref_geom_max(ref_geom) - 1) = REF_EMPTY;
  ref_geom_blank(ref_geom) = 0;
  if (NULL != (void *)(ref_geom->ref_adj))
    RSS(ref_adj_free(ref_geom->ref_adj), "free to prevent leak");
  RSS(ref_adj_create(&(ref_geom->ref_adj)), "create ref_adj for ref_geom");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_create(REF_GEOM *ref_geom_ptr) {
  REF_GEOM ref_geom;

  (*ref_geom_ptr) = NULL;

  ref_malloc(*ref_geom_ptr, 1, REF_GEOM_STRUCT);

  ref_geom = (*ref_geom_ptr);

  ref_geom_max(ref_geom) = 10;

  ref_malloc(ref_geom->descr, REF_GEOM_DESCR_SIZE * ref_geom_max(ref_geom),
             REF_INT);
  ref_malloc(ref_geom->param, 2 * ref_geom_max(ref_geom), REF_DBL);
  ref_geom->ref_adj = (REF_ADJ)NULL;
  RSS(ref_geom_initialize(ref_geom), "init geom list");

  ref_geom->uv_area_sign = NULL;
  ref_geom->initial_cell_height = NULL;
  ref_geom->face_min_length = NULL;
  ref_geom->face_seg_per_rad = NULL;
  ref_geom->segments_per_radian_of_curvature = 2.0;
  ref_geom->segments_per_bounding_box_diagonal = 10.0;
  ref_geom->tolerance_protection = 100.0;
  ref_geom->gap_protection = 10.0;

  ref_geom->nnode = REF_EMPTY;
  ref_geom->nedge = REF_EMPTY;
  ref_geom->nface = REF_EMPTY;
  ref_geom->zip_pcurve = REF_FALSE;
  ref_geom->effective = REF_FALSE;
  ref_geom->effective_curvature = REF_TRUE;
  ref_geom->manifold = REF_TRUE;
  ref_geom->contex_owned = REF_TRUE;
  ref_geom->context = NULL;
  RSS(ref_egads_open(ref_geom), "open egads");
  ref_geom->model = NULL;
  ref_geom->body = NULL;
  ref_geom->faces = NULL;
  ref_geom->edges = NULL;
  ref_geom->nodes = NULL;
  ref_geom->pcurves = NULL;
  ref_geom->e2f = NULL;

  ref_geom->cad_data_size = 0;
  ref_geom->cad_data = (REF_BYTE *)NULL;

  ref_geom->meshlink = NULL;
  ref_geom->meshlink_projection = NULL;

  ref_geom->ref_facelift = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_free(REF_GEOM ref_geom) {
  if (NULL == (void *)ref_geom) return REF_NULL;
  ref_facelift_free(ref_geom_facelift(ref_geom));
  ref_free(ref_geom->cad_data);
  ref_free(ref_geom->e2f);
  if (ref_geom->contex_owned)
    RSS(ref_egads_close(ref_geom), "close egads contex");
  RSS(ref_adj_free(ref_geom->ref_adj), "adj free");
  ref_free(ref_geom->face_seg_per_rad);
  ref_free(ref_geom->face_min_length);
  ref_free(ref_geom->initial_cell_height);
  ref_free(ref_geom->uv_area_sign);
  ref_free(ref_geom->param);
  ref_free(ref_geom->descr);
  ref_free(ref_geom);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_deep_copy(REF_GEOM *ref_geom_ptr, REF_GEOM original) {
  REF_GEOM ref_geom;
  REF_INT geom, i;
  (*ref_geom_ptr) = NULL;

  ref_malloc(*ref_geom_ptr, 1, REF_GEOM_STRUCT);

  ref_geom = (*ref_geom_ptr);

  ref_geom_n(ref_geom) = ref_geom_n(original);
  ref_geom_max(ref_geom) = ref_geom_max(original);

  ref_malloc(ref_geom->descr, REF_GEOM_DESCR_SIZE * ref_geom_max(ref_geom),
             REF_INT);
  ref_malloc(ref_geom->param, 2 * ref_geom_max(ref_geom), REF_DBL);
  ref_geom->uv_area_sign = NULL;
  ref_geom->initial_cell_height = NULL;
  ref_geom->face_min_length = NULL;
  ref_geom->face_seg_per_rad = NULL;
  ref_geom->segments_per_radian_of_curvature =
      original->segments_per_radian_of_curvature;
  ref_geom->segments_per_bounding_box_diagonal =
      original->segments_per_bounding_box_diagonal;
  ref_geom->tolerance_protection = original->tolerance_protection;
  ref_geom->gap_protection = original->gap_protection;

  for (geom = 0; geom < ref_geom_max(ref_geom); geom++)
    for (i = 0; i < REF_GEOM_DESCR_SIZE; i++)
      ref_geom_descr(ref_geom, i, geom) = ref_geom_descr(original, i, geom);
  ref_geom_blank(ref_geom) = ref_geom_blank(original);
  for (geom = 0; geom < ref_geom_max(ref_geom); geom++)
    for (i = 0; i < 2; i++)
      ref_geom_param(ref_geom, i, geom) = ref_geom_param(original, i, geom);

  RSS(ref_adj_deep_copy(&(ref_geom->ref_adj), original->ref_adj),
      "deep copy ref_adj for ref_geom");

  ref_geom->zip_pcurve = original->zip_pcurve;

  ref_geom->contex_owned = REF_FALSE;
  RSS(ref_geom_share_context(ref_geom, original), "share egads");

  ref_geom->cad_data_size = 0;
  ref_geom->cad_data = (REF_BYTE *)NULL;

  ref_geom->meshlink = NULL;
  ref_geom->meshlink_projection = NULL;

  ref_geom->ref_facelift = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_share_context(REF_GEOM ref_geom_recipient,
                                  REF_GEOM ref_geom_donor) {
  if (ref_geom_recipient->contex_owned)
    RSS(ref_egads_close(ref_geom_recipient), "close egads contex");

  ref_geom_recipient->contex_owned = REF_FALSE;
  ref_geom_recipient->nnode = ref_geom_donor->nnode;
  ref_geom_recipient->nedge = ref_geom_donor->nedge;
  ref_geom_recipient->nface = ref_geom_donor->nface;
  ref_geom_recipient->effective = ref_geom_donor->effective;
  ref_geom_recipient->effective_curvature = ref_geom_donor->effective_curvature;
  ref_geom_recipient->manifold = ref_geom_donor->manifold;
  ref_geom_recipient->context = ref_geom_donor->context;
  ref_geom_recipient->model = ref_geom_donor->model;
  ref_geom_recipient->body = ref_geom_donor->body;
  ref_geom_recipient->faces = ref_geom_donor->faces;
  ref_geom_recipient->edges = ref_geom_donor->edges;
  ref_geom_recipient->nodes = ref_geom_donor->nodes;
  ref_geom_recipient->pcurves = ref_geom_donor->pcurves;

  if (NULL == ref_geom_donor->e2f) {
    ref_geom_recipient->e2f = NULL;
  } else {
    REF_INT i;
    ref_malloc(ref_geom_recipient->e2f, 2 * ref_geom_recipient->nedge, REF_INT);
    for (i = 0; i < 2 * ref_geom_recipient->nedge; i++)
      ref_geom_recipient->e2f[i] = ref_geom_donor->e2f[i];
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_pack(REF_GEOM ref_geom, REF_INT *o2n) {
  REF_INT geom, compact, i;
  compact = 0;
  each_ref_geom(ref_geom, geom) {
    for (i = 0; i < REF_GEOM_DESCR_SIZE; i++)
      ref_geom_descr(ref_geom, i, compact) = ref_geom_descr(ref_geom, i, geom);
    ref_geom_node(ref_geom, compact) = o2n[ref_geom_node(ref_geom, geom)];
    for (i = 0; i < 2; i++)
      ref_geom_param(ref_geom, i, compact) = ref_geom_param(ref_geom, i, geom);
    compact++;
  }
  REIS(compact, ref_geom_n(ref_geom), "count mismatch");
  if (ref_geom_n(ref_geom) < ref_geom_max(ref_geom)) {
    for (geom = ref_geom_n(ref_geom); geom < ref_geom_max(ref_geom); geom++) {
      ref_geom_type(ref_geom, geom) = REF_EMPTY;
      ref_geom_id(ref_geom, geom) = geom + 1;
    }
    ref_geom_id(ref_geom, ref_geom_max(ref_geom) - 1) = REF_EMPTY;
    ref_geom_blank(ref_geom) = ref_geom_n(ref_geom);
  } else {
    ref_geom_blank(ref_geom) = REF_EMPTY;
  }
  RSS(ref_adj_free(ref_geom->ref_adj), "free to prevent leak");
  RSS(ref_adj_create(&(ref_geom->ref_adj)), "create ref_adj for ref_geom");

  each_ref_geom(ref_geom, geom) {
    RSS(ref_adj_add(ref_geom->ref_adj, ref_geom_node(ref_geom, geom), geom),
        "register geom");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_uv_area(REF_GEOM ref_geom, REF_INT *nodes,
                            REF_DBL *uv_area) {
  REF_DBL uv0[2], uv1[2], uv2[2];
  REF_INT sens;
  double a, b, c, d;
  RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_FACE, uv0, &sens),
      "uv0");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_FACE, uv1, &sens),
      "uv1");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[2], nodes, REF_GEOM_FACE, uv2, &sens),
      "uv2");
  a = uv0[0] - uv2[0];
  b = uv0[1] - uv2[1];
  c = uv1[0] - uv2[0];
  d = uv1[1] - uv2[1];
  *uv_area = 0.5 * (a * d - b * c);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_uv_area_sign(REF_GRID ref_grid, REF_INT id, REF_DBL *sign) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  if (NULL == ((ref_geom)->uv_area_sign)) {
    REF_CELL ref_cell = ref_grid_tri(ref_grid);
    REF_INT face;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_DBL uv_area;
    if (REF_EMPTY == ref_geom->nface)
      RSS(ref_geom_infer_nedge_nface(ref_grid), "infer counts");
    ref_malloc_init(ref_geom->uv_area_sign, ref_geom->nface, REF_DBL, 0.0);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      face = nodes[3];
      if (face < 1 || ref_geom->nface < face) continue;
      RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area), "uv area");
      if (uv_area < 0.0) {
        ((ref_geom)->uv_area_sign)[face - 1] -= 1.0;
      } else {
        ((ref_geom)->uv_area_sign)[face - 1] += 1.0;
      }
    }
    for (face = 0; face < ref_geom->nface; face++) {
      if (((ref_geom)->uv_area_sign)[face] < 0.0) {
        ((ref_geom)->uv_area_sign)[face] = -1.0;
      } else {
        ((ref_geom)->uv_area_sign)[face] = 1.0;
      }
    }
  }

  if (id < 1 || id > ref_geom->nface) return REF_INVALID;
  *sign = ((ref_geom)->uv_area_sign)[id - 1];

  return REF_SUCCESS;
}

REF_STATUS ref_geom_uv_area_report(REF_GRID ref_grid) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT geom, id, min_id, max_id;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL no_cell;
  REF_DBL uv_area, total_uv_area, min_uv_area, max_uv_area, sign_uv_area;
  REF_INT n_neg, n_pos;

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++) {
    no_cell = REF_TRUE;
    total_uv_area = 0.0;
    min_uv_area = 0.0;
    max_uv_area = 0.0;
    n_neg = 0;
    n_pos = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (id == nodes[3]) {
        RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area), "uv area");
        total_uv_area += uv_area;
        if (no_cell) {
          min_uv_area = uv_area;
          max_uv_area = uv_area;
          no_cell = REF_FALSE;
        } else {
          min_uv_area = MIN(min_uv_area, uv_area);
          max_uv_area = MAX(max_uv_area, uv_area);
        }
        if (uv_area < 0.0) {
          n_neg++;
        } else {
          n_pos++;
        }
      }
    }
    if (!no_cell) {
      RSS(ref_geom_uv_area_sign(ref_grid, id, &sign_uv_area), "sign");
      printf("face%5d: %4.1f %9.2e total (%10.3e,%10.3e) %d + %d -\n", id,
             sign_uv_area, total_uv_area, min_uv_area, max_uv_area, n_pos,
             n_neg);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_inspect(REF_GEOM ref_geom) {
  REF_INT geom;
  printf("ref_geom = %p\n", (void *)ref_geom);
  printf(" n = %d, max = %d\n", ref_geom_n(ref_geom), ref_geom_max(ref_geom));
  for (geom = 0; geom < ref_geom_max(ref_geom); geom++) {
    switch (ref_geom_type(ref_geom, geom)) {
      case REF_GEOM_NODE:
        printf("%d node: %d id, %d jump, %d degen, %d global\n", geom,
               ref_geom_id(ref_geom, geom), ref_geom_jump(ref_geom, geom),
               ref_geom_degen(ref_geom, geom), ref_geom_node(ref_geom, geom));
        break;
      case REF_GEOM_EDGE:
        printf("%d edge: %d id, %d jump, %d degen, %d global, t=%e\n", geom,
               ref_geom_id(ref_geom, geom), ref_geom_jump(ref_geom, geom),
               ref_geom_node(ref_geom, geom), ref_geom_degen(ref_geom, geom),
               ref_geom_param(ref_geom, 0, geom));
        break;
      case REF_GEOM_FACE:
        printf("%d face: %d id, %d jump, %d degen, %d global, uv= %e %e\n",
               geom, ref_geom_id(ref_geom, geom), ref_geom_jump(ref_geom, geom),
               ref_geom_node(ref_geom, geom), ref_geom_degen(ref_geom, geom),
               ref_geom_param(ref_geom, 0, geom),
               ref_geom_param(ref_geom, 1, geom));
        break;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tattle(REF_GEOM ref_geom, REF_INT node) {
  REF_INT item, geom;

  printf(" tattle on node = %d\n", node);
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    switch (ref_geom_type(ref_geom, geom)) {
      case REF_GEOM_NODE:
        printf("%d node: %d id, %d jump, %d degen, %d global\n", geom,
               ref_geom_id(ref_geom, geom), ref_geom_jump(ref_geom, geom),
               ref_geom_degen(ref_geom, geom), ref_geom_node(ref_geom, geom));
        break;
      case REF_GEOM_EDGE:
        printf("%d edge: %d id, %d jump, %d degen, %d global, t=%e\n", geom,
               ref_geom_id(ref_geom, geom), ref_geom_jump(ref_geom, geom),
               ref_geom_degen(ref_geom, geom), ref_geom_node(ref_geom, geom),
               ref_geom_param(ref_geom, 0, geom));
        break;
      case REF_GEOM_FACE:
        printf("%d face: %d id, %d jump, %d degen, %d global, uv= %e %e\n",
               geom, ref_geom_id(ref_geom, geom), ref_geom_jump(ref_geom, geom),
               ref_geom_degen(ref_geom, geom), ref_geom_node(ref_geom, geom),
               ref_geom_param(ref_geom, 0, geom),
               ref_geom_param(ref_geom, 1, geom));
        break;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_supported(REF_GEOM ref_geom, REF_INT node,
                              REF_BOOL *has_support) {
  *has_support = !ref_adj_empty(ref_geom_adj(ref_geom), node);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_tri_supported(REF_GEOM ref_geom, REF_INT *nodes,
                                  REF_BOOL *has_support) {
  REF_INT node, id, geom;
  REF_STATUS status;
  *has_support = REF_FALSE;

  node = nodes[0];
  id = nodes[3];
  status = ref_geom_find(ref_geom, node, REF_GEOM_FACE, id, &geom);
  if (REF_NOT_FOUND == status) { /* no geom support */
    *has_support = REF_FALSE;
    return REF_SUCCESS;
  }
  RSS(status, "error testing geom support");
  *has_support = REF_TRUE;
  return REF_SUCCESS;
}

REF_STATUS ref_geom_id_supported(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                                 REF_INT id, REF_BOOL *has_support) {
  REF_INT geom;
  REF_STATUS status;
  *has_support = REF_FALSE;

  status = ref_geom_find(ref_geom, node, type, id, &geom);
  if (REF_NOT_FOUND == status) { /* no geom support */
    *has_support = REF_FALSE;
    return REF_SUCCESS;
  }
  RSS(status, "error testing geom support");
  *has_support = REF_TRUE;
  return REF_SUCCESS;
}

static REF_STATUS ref_geom_grow(REF_GEOM ref_geom) {
  REF_INT geom;
  REF_INT orig, chunk;
  REF_INT max_limit = REF_INT_MAX / 3;

  if (REF_EMPTY != ref_geom_blank(ref_geom)) {
    return REF_SUCCESS;
  }

  RAS(ref_geom_max(ref_geom) != max_limit,
      "the number of geoms is too large for integers, cannot grow");
  orig = ref_geom_max(ref_geom);
  /* geometric growth for efficiency */
  chunk = MAX(1000, (REF_INT)(1.5 * (REF_DBL)orig));

  /* try to keep under 32-bit limit */
  RAS(max_limit - orig > 0, "chunk limit at max");
  chunk = MIN(chunk, max_limit - orig);

  ref_geom_max(ref_geom) = orig + chunk;

  ref_realloc(ref_geom->descr, REF_GEOM_DESCR_SIZE * ref_geom_max(ref_geom),
              REF_INT);
  ref_realloc(ref_geom->param, 2 * ref_geom_max(ref_geom), REF_DBL);

  for (geom = orig; geom < ref_geom_max(ref_geom); geom++) {
    ref_geom_type(ref_geom, geom) = REF_EMPTY;
    ref_geom_id(ref_geom, geom) = geom + 1;
  }
  ref_geom_id(ref_geom, ref_geom_max(ref_geom) - 1) = REF_EMPTY;
  ref_geom_blank(ref_geom) = orig;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_add_with_descr(REF_GEOM ref_geom, REF_INT *descr,
                                   REF_DBL *param) {
  REF_INT type, id, gref, jump, degen, node, geom;
  type = descr[REF_GEOM_DESCR_TYPE];
  id = descr[REF_GEOM_DESCR_ID];
  gref = descr[REF_GEOM_DESCR_GREF];
  jump = descr[REF_GEOM_DESCR_JUMP];
  degen = descr[REF_GEOM_DESCR_DEGEN];
  node = descr[REF_GEOM_DESCR_NODE];
  RSS(ref_geom_add(ref_geom, node, type, id, param), "geom add");
  RSS(ref_geom_find(ref_geom, node, type, id, &geom), "geom find");
  ref_geom_gref(ref_geom, geom) = gref;
  ref_geom_degen(ref_geom, geom) = degen;
  ref_geom_jump(ref_geom, geom) = jump;
  return REF_SUCCESS;
}
REF_STATUS ref_geom_add(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                        REF_INT id, REF_DBL *param) {
  REF_INT geom;
  REF_STATUS status;

  if (type < 0 || 2 < type) return REF_INVALID;

  status = ref_geom_find(ref_geom, node, type, id, &geom);
  RXS(status, REF_NOT_FOUND, "find failed");

  if (REF_SUCCESS == status) {
    if (type > 0) ref_geom_param(ref_geom, 0, geom) = param[0];
    if (type > 1) ref_geom_param(ref_geom, 1, geom) = param[1];
    return REF_SUCCESS;
  }

  if (REF_EMPTY == ref_geom_blank(ref_geom)) {
    RSS(ref_geom_grow(ref_geom), "grow add");
  }

  geom = ref_geom_blank(ref_geom);
  ref_geom_blank(ref_geom) = ref_geom_id(ref_geom, geom);

  ref_geom_type(ref_geom, geom) = type;
  ref_geom_id(ref_geom, geom) = id;
  ref_geom_gref(ref_geom, geom) = id; /* assume same until set */
  ref_geom_jump(ref_geom, geom) = 0;
  ref_geom_degen(ref_geom, geom) = 0;
  ref_geom_node(ref_geom, geom) = node;

  ref_geom_param(ref_geom, 0, geom) = 0.0;
  ref_geom_param(ref_geom, 1, geom) = 0.0;
  if (type > 0) ref_geom_param(ref_geom, 0, geom) = param[0];
  if (type > 1) ref_geom_param(ref_geom, 1, geom) = param[1];

  RSB(ref_adj_add(ref_geom->ref_adj, node, geom), "register geom", {
    printf("register node %d geom %d type %d id %d\n",
           ref_geom_node(ref_geom, geom), geom, ref_geom_type(ref_geom, geom),
           ref_geom_id(ref_geom, geom));
  });

  ref_geom_n(ref_geom)++;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_remove(REF_GEOM ref_geom, REF_INT geom) {
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT node;

  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;
  if (REF_EMPTY == ref_geom_type(ref_geom, geom)) return REF_INVALID;

  node = ref_geom_node(ref_geom, geom);
  RSS(ref_adj_remove(ref_adj, node, geom), "unregister geom");

  ref_geom_type(ref_geom, geom) = REF_EMPTY;
  ref_geom_id(ref_geom, geom) = ref_geom_blank(ref_geom);
  ref_geom_blank(ref_geom) = geom;
  ref_geom_n(ref_geom)--;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_remove_all(REF_GEOM ref_geom, REF_INT node) {
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT item, geom;

  item = ref_adj_first(ref_adj, node);
  while (ref_adj_valid(item)) {
    geom = ref_adj_item_ref(ref_adj, item);

    RSS(ref_geom_remove(ref_geom, geom), "remove");

    item = ref_adj_first(ref_adj, node);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_remove_without_cell(REF_GRID ref_grid, REF_INT node) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_CELL ref_cell;
  REF_INT type, item, geom, group;
  REF_BOOL supported_by_cell;

  item = ref_adj_first(ref_adj, node);
  while (ref_adj_valid(item)) {
    geom = ref_adj_item_ref(ref_adj, item);
    type = ref_geom_type(ref_geom, geom);
    switch (type) {
      case REF_GEOM_NODE:
        item = ref_adj_item_next(ref_adj, item);
        break;
      case REF_GEOM_EDGE:
        supported_by_cell = REF_FALSE;
        each_ref_grid_edge_ref_cell(ref_grid, group, ref_cell) {
          supported_by_cell =
              (supported_by_cell || !ref_cell_node_empty(ref_cell, node));
        }
        if (supported_by_cell) {
          item = ref_adj_item_next(ref_adj, item);
        } else {
          RSS(ref_geom_remove(ref_geom, geom), "remove");
          item = ref_adj_first(ref_adj, node);
        }
        break;
      case REF_GEOM_FACE:
        supported_by_cell = REF_FALSE;
        each_ref_grid_face_ref_cell(ref_grid, group, ref_cell) {
          supported_by_cell =
              (supported_by_cell || !ref_cell_node_empty(ref_cell, node));
        }
        if (supported_by_cell) {
          item = ref_adj_item_next(ref_adj, item);
        } else {
          RSS(ref_geom_remove(ref_geom, geom), "remove");
          item = ref_adj_first(ref_adj, node);
        }
        break;
      default:
        RSS(REF_IMPLEMENT, "can't to geom type yet");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_is_a(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                         REF_BOOL *it_is) {
  REF_INT item, geom;
  *it_is = REF_FALSE;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    if (type == ref_geom_type(ref_geom, geom)) {
      *it_is = REF_TRUE;
      return REF_SUCCESS;
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_geom_unique_id(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                              REF_INT *id) {
  REF_INT item, geom;
  REF_BOOL found_one;
  found_one = REF_FALSE;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    if (type == ref_geom_type(ref_geom, geom)) {
      if (found_one) return REF_INVALID; /* second one makes invalid */
      found_one = REF_TRUE;
      *id = ref_geom_id(ref_geom, geom);
    }
  }
  if (found_one) return REF_SUCCESS;
  return REF_NOT_FOUND;
}

REF_STATUS ref_geom_find(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                         REF_INT id, REF_INT *found) {
  REF_INT item, geom;
  *found = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    if (type == ref_geom_type(ref_geom, geom) &&
        id == ref_geom_id(ref_geom, geom)) {
      *found = geom;
      return REF_SUCCESS;
    }
  }
  return REF_NOT_FOUND;
}

REF_STATUS ref_geom_tuv(REF_GEOM ref_geom, REF_INT node, REF_INT type,
                        REF_INT id, REF_DBL *param) {
  REF_INT geom;

  RSS(ref_geom_find(ref_geom, node, type, id, &geom), "not found");

  REIS(0, ref_geom_jump(ref_geom, geom), "use ref_geom_cell_tuv for jumps");
  REIS(0, ref_geom_degen(ref_geom, geom), "use ref_geom_cell_tuv for degen");

  if (type > 0) param[0] = ref_geom_param(ref_geom, 0, geom);
  if (type > 1) param[1] = ref_geom_param(ref_geom, 1, geom);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_cell_tuv_supported(REF_GEOM ref_geom, REF_INT *nodes,
                                       REF_INT type, REF_BOOL *supported) {
  REF_INT node_per;
  REF_INT id, geom0, geom1, geom2;
  REF_BOOL tri_supported;

  *supported = REF_TRUE;

  RAS(1 <= type && type <= 2, "type not allowed");
  node_per = type + 1;
  id = nodes[node_per];

  /* protects unsupported meshlink tri */
  if (REF_GEOM_FACE == type) {
    RSS(ref_geom_tri_supported(ref_geom, nodes, &tri_supported), "tri support");
    if (!tri_supported) { /* no geom support */
      *supported = 1.0;
      return REF_SUCCESS;
    }
  }

  switch (type) {
    case REF_GEOM_EDGE:
      RSS(ref_geom_find(ref_geom, nodes[0], type, id, &geom0), "not found");
      RSS(ref_geom_find(ref_geom, nodes[1], type, id, &geom1), "not found");

      if ((0 != ref_geom_jump(ref_geom, geom0) ||
           0 != ref_geom_degen(ref_geom, geom0)) &&
          (0 != ref_geom_jump(ref_geom, geom1) ||
           0 != ref_geom_degen(ref_geom, geom1))) {
        *supported = REF_FALSE;
      }
      break;
    case REF_GEOM_FACE:
      RSS(ref_geom_find(ref_geom, nodes[0], type, id, &geom0), "not found");
      RSS(ref_geom_find(ref_geom, nodes[1], type, id, &geom1), "not found");
      RSS(ref_geom_find(ref_geom, nodes[2], type, id, &geom2), "not found");
      if ((0 != ref_geom_jump(ref_geom, geom0) ||
           0 != ref_geom_degen(ref_geom, geom0)) &&
          (0 != ref_geom_jump(ref_geom, geom1) ||
           0 != ref_geom_degen(ref_geom, geom1)) &&
          (0 != ref_geom_jump(ref_geom, geom2) ||
           0 != ref_geom_degen(ref_geom, geom2))) {
        *supported = REF_FALSE;
      }
      break;
    default:
      RSS(REF_IMPLEMENT, "can't to geom type yet");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_cell_tuv(REF_GEOM ref_geom, REF_INT node, REF_INT *nodes,
                             REF_INT type, REF_DBL *param, REF_INT *sens) {
  REF_INT node_per;
  REF_INT id, edgeid, geom, from, from_geom;
  REF_INT node_index, cell_node;
  double trange[2], uv[2], uv0[2], uv1[2], uvtmin[2], uvtmax[2];
  REF_DBL from_param[2], t;
  REF_DBL dist0, dist1;
  REF_INT hits;

  RAS(1 <= type && type <= 2, "type not allowed");
  node_per = type + 1;
  id = nodes[node_per];
  node_index = REF_EMPTY;
  for (cell_node = 0; cell_node < node_per; cell_node++) {
    if (node == nodes[cell_node]) {
      REIS(REF_EMPTY, node_index, "node found twice in nodes");
      node_index = cell_node;
    }
  }
  RAS(REF_EMPTY != node_index, "node not found in nodes");

  RSB(ref_geom_find(ref_geom, node, type, id, &geom), "not found", {
    printf(" %d type %d id\n", type, id);
    ref_geom_tattle(ref_geom, node);
  });

  if (0 == ref_geom_jump(ref_geom, geom) &&
      0 == ref_geom_degen(ref_geom, geom)) {
    if (type > 0) param[0] = ref_geom_param(ref_geom, 0, geom);
    if (type > 1) param[1] = ref_geom_param(ref_geom, 1, geom);
    *sens = 0;
    return REF_SUCCESS;
  }

  switch (type) {
    case REF_GEOM_EDGE:
      RSS(ref_egads_edge_trange(ref_geom, id, trange), "trange");
      from = nodes[1 - node_index];
      RSS(ref_geom_tuv(ref_geom, from, type, id, from_param), "from tuv");
      dist0 = from_param[0] - trange[0];
      dist1 = trange[1] - from_param[0];
      if (dist0 < 0.0 || dist1 < 0.0) {
        printf(" from t = %e %e %e, dist = %e %e\n", trange[0], from_param[0],
               trange[1], dist0, dist1);
        THROW("from node not in trange");
      }
      if (dist0 < dist1) {
        *sens = 1;
        param[0] = trange[0];
      } else {
        *sens = -1;
        param[0] = trange[1];
      }
      break;
    case REF_GEOM_FACE:
      if (0 == ref_geom_degen(ref_geom, geom)) {
        from = REF_EMPTY;
        for (cell_node = 0; cell_node < node_per; cell_node++) {
          RSS(ref_geom_find(ref_geom, nodes[cell_node], type, id, &from_geom),
              "not found");
          if (node_index != cell_node &&
              0 == ref_geom_jump(ref_geom, from_geom) &&
              0 == ref_geom_degen(ref_geom, from_geom)) {
            from = nodes[cell_node];
          }
        }
        RAB(REF_EMPTY != from, "can't find from tuv in tri cell", {
          ref_geom_tattle(ref_geom, nodes[0]);
          ref_geom_tattle(ref_geom, nodes[1]);
          ref_geom_tattle(ref_geom, nodes[2]);
          printf("faceid %d node %d node_index %d\n", id, node, node_index);
        });
        edgeid = ref_geom_jump(ref_geom, geom);
        RSS(ref_geom_tuv(ref_geom, from, REF_GEOM_FACE, id, uv), "from uv");
        RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_EDGE, edgeid, &t), "edge t0");
        RSS(ref_egads_edge_face_uv(ref_geom, edgeid, id, 1, t, uv0), "uv 1");
        RSS(ref_egads_edge_face_uv(ref_geom, edgeid, id, -1, t, uv1), "uv -1");
        dist0 = sqrt(pow(uv0[0] - uv[0], 2) + pow(uv0[1] - uv[1], 2));
        dist1 = sqrt(pow(uv1[0] - uv[0], 2) + pow(uv1[1] - uv[1], 2));
        if (dist0 < dist1) {
          *sens = 1;
          param[0] = uv0[0];
          param[1] = uv0[1];
        } else {
          *sens = -1;
          param[0] = uv1[0];
          param[1] = uv1[1];
        }
      } else {
        uv0[0] = 0.0;
        uv0[1] = 0.0;
        hits = 0;
        for (cell_node = 0; cell_node < node_per; cell_node++) {
          RSS(ref_geom_find(ref_geom, nodes[cell_node], type, id, &from_geom),
              "not found");
          if (0 == ref_geom_jump(ref_geom, from_geom) &&
              0 == ref_geom_degen(ref_geom, from_geom)) {
            RSS(ref_geom_tuv(ref_geom, nodes[cell_node], REF_GEOM_FACE, id, uv),
                "from uv");
            uv0[0] += uv[0];
            uv0[1] += uv[1];
            hits++;
          }
        }
        RAS(0 < hits, "no seed uv found for DEGEN");
        uv0[0] /= (REF_DBL)hits;
        uv0[1] /= (REF_DBL)hits;

        *sens = 0;
        edgeid = ABS(ref_geom_degen(ref_geom, geom));
        RSS(ref_egads_edge_trange(ref_geom, edgeid, trange), "trange");
        RSS(ref_egads_edge_face_uv(ref_geom, edgeid,
                                   ref_geom_id(ref_geom, geom), *sens,
                                   trange[0], uvtmin),
            "uv t min");
        RSS(ref_egads_edge_face_uv(ref_geom, edgeid,
                                   ref_geom_id(ref_geom, geom), *sens,
                                   trange[1], uvtmax),
            "uv t max");
        /* edgeid sign convention defined in ref_geom_mark_jump_degen */
        if (0 < ref_geom_degen(ref_geom, geom)) {
          param[0] = ref_geom_param(ref_geom, 0, geom);
          param[1] = uv0[1];
          param[1] = MAX(param[1], MIN(uvtmin[1], uvtmax[1]));
          param[1] = MIN(param[1], MAX(uvtmin[1], uvtmax[1]));
        } else {
          param[0] = uv0[0];
          param[0] = MAX(param[0], MIN(uvtmin[0], uvtmax[0]));
          param[0] = MIN(param[0], MAX(uvtmin[0], uvtmax[0]));
          param[1] = ref_geom_param(ref_geom, 1, geom);
        }
      }
      break;
    default:
      RSS(REF_IMPLEMENT, "can't to geom type yet");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_geom_eval_edge_face_uv(REF_GRID ref_grid,
                                             REF_INT edge_geom) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT node, cell_item, geom_item, cell, face_geom;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  double t;
  double uv[2], edgeuv[2];
  int sense;
  REF_INT faceid;
  REF_BOOL have_jump;

  if (edge_geom < 0 || ref_geom_max(ref_geom) <= edge_geom) return REF_INVALID;
  if (REF_GEOM_EDGE != ref_geom_type(ref_geom, edge_geom)) return REF_INVALID;

  t = ref_geom_param(ref_geom, 0, edge_geom);
  node = ref_geom_node(ref_geom, edge_geom);

  have_jump = REF_FALSE;
  each_ref_adj_node_item_with_ref(ref_adj, node, geom_item, face_geom) {
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, face_geom)) {
      have_jump = have_jump || (0 != ref_geom_jump(ref_geom, face_geom));
    }
  }

  if (have_jump) {
    /* uv update at jump not needed, should always depend on cell_c2n */
    /* keeping for consistency with non-jump */
    each_ref_cell_having_node(ref_cell, node, cell_item, cell) {
      RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
      faceid = nodes[3];
      RSS(ref_geom_cell_tuv(ref_geom, node, nodes, REF_GEOM_FACE, uv, &sense),
          "cell uv");
      if (1 == sense) { /* sense to use is arbitrary */
        each_ref_adj_node_item_with_ref(ref_adj, node, geom_item, face_geom) {
          if (REF_GEOM_FACE == ref_geom_type(ref_geom, face_geom) &&
              faceid == ref_geom_id(ref_geom, face_geom)) {
            ref_geom_param(ref_geom, 0, face_geom) = uv[0];
            ref_geom_param(ref_geom, 1, face_geom) = uv[1];
          }
        }
      }
    }
  } else {
    each_ref_adj_node_item_with_ref(ref_adj, node, geom_item, face_geom) {
      if (REF_GEOM_FACE == ref_geom_type(ref_geom, face_geom)) {
        faceid = ref_geom_id(ref_geom, face_geom);
        sense = 0;
        RSS(ref_egads_edge_face_uv(ref_geom, ref_geom_id(ref_geom, edge_geom),
                                   faceid, sense, t, edgeuv),
            "edge uv");
        ref_geom_param(ref_geom, 0, face_geom) = edgeuv[0];
        ref_geom_param(ref_geom, 1, face_geom) = edgeuv[1];
      }
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_geom_add_constrain_inside_midnode(REF_GRID ref_grid,
                                                        REF_INT *nodes,
                                                        REF_INT new_node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT type, id;
  REF_DBL uv[2], uv0[2], uv1[2], uv2[2];
  REF_INT sense;

  type = REF_GEOM_FACE;
  id = nodes[3];

  RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, type, uv0, &sense),
      "cell uv0");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, type, uv1, &sense),
      "cell uv1");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[2], nodes, type, uv2, &sense),
      "cell uv2");

  uv[0] = (uv0[0] + uv1[0] + uv2[0]) / 3.0;
  uv[1] = (uv0[1] + uv1[1] + uv2[1]) / 3.0;
  RSS(ref_geom_add(ref_geom, new_node, type, id, uv), "new geom");

  RSS(ref_egads_eval_at(ref_geom, type, id, uv,
                        ref_node_xyz_ptr(ref_node, new_node), NULL),
      "eval");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_between_face_area(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_INT new_node,
                                      const char *msg) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT nodes2[REF_CELL_MAX_SIZE_PER];
  REF_INT nodes3[REF_CELL_MAX_SIZE_PER];
  REF_INT type, id, node;
  REF_INT ncell, cells[2];
  REF_DBL uv_area2, uv_area3, uv_area_sign;
  REF_DBL uv_area20, uv_area21;
  REF_DBL uv_area30, uv_area31;

  type = REF_GEOM_FACE;
  RSS(ref_geom_unique_id(ref_geom, new_node, type, &id), "unique face id");

  RSS(ref_geom_uv_area_sign(ref_grid, id, &uv_area_sign), "uv area sign");

  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cells), "list");
  REIS(2, ncell, "expected two tri for box2 nodes");
  RSS(ref_cell_nodes(ref_cell, cells[0], nodes2), "cell nodes");
  RSS(ref_cell_nodes(ref_cell, cells[1], nodes3), "cell nodes");

  RSS(ref_geom_uv_area(ref_geom, nodes2, &uv_area2), "uv area 2");
  uv_area2 *= uv_area_sign;
  RSS(ref_geom_uv_area(ref_geom, nodes3, &uv_area3), "uv area 3");
  uv_area3 *= uv_area_sign;

  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    if (node0 == nodes2[node]) nodes2[node] = new_node;
  RSS(ref_geom_uv_area(ref_geom, nodes2, &uv_area20), "uv area 2 node 0");
  uv_area20 *= uv_area_sign;
  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    if (new_node == nodes2[node]) nodes2[node] = node0;
  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    if (node1 == nodes2[node]) nodes2[node] = new_node;
  RSS(ref_geom_uv_area(ref_geom, nodes2, &uv_area21), "uv area 2 node 1");
  uv_area21 *= uv_area_sign;
  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    if (new_node == nodes2[node]) nodes2[node] = node1;

  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    if (node0 == nodes3[node]) nodes3[node] = new_node;
  RSS(ref_geom_uv_area(ref_geom, nodes3, &uv_area30), "uv area 3 node 0");
  uv_area30 *= uv_area_sign;
  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    if (new_node == nodes3[node]) nodes3[node] = node0;
  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    if (node1 == nodes3[node]) nodes3[node] = new_node;
  RSS(ref_geom_uv_area(ref_geom, nodes3, &uv_area31), "uv area 3 node 1");
  uv_area31 *= uv_area_sign;
  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    if (new_node == nodes3[node]) nodes3[node] = node1;

  if (uv_area2 <= 1e-20 || uv_area3 <= 1e-20 || uv_area20 <= 0 ||
      uv_area21 <= 0 || uv_area30 <= 0 || uv_area31 <= 0) {
    REF_DBL uv[2];
    REF_DBL u0 = 0.0, v0 = 0.0;
    REF_INT sense;
    printf("%s orig uv area %e %e new %e %e %e %e\n", msg, uv_area2, uv_area3,
           uv_area20, uv_area21, uv_area30, uv_area31);
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes2[node]) nodes2[node] = new_node;
    RSS(ref_geom_cell_tuv(ref_geom, nodes2[0], nodes2, type, uv, &sense),
        "cell uv0");
    printf("20uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    RSS(ref_geom_cell_tuv(ref_geom, nodes2[1], nodes2, type, uv, &sense),
        "cell uv0");
    printf("20uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    RSS(ref_geom_cell_tuv(ref_geom, nodes2[2], nodes2, type, uv, &sense),
        "cell uv0");
    printf("20uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes2[node]) nodes2[node] = node0;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes2[node]) nodes2[node] = new_node;
    RSS(ref_geom_cell_tuv(ref_geom, nodes2[0], nodes2, type, uv, &sense),
        "cell uv0");
    printf("21uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    RSS(ref_geom_cell_tuv(ref_geom, nodes2[1], nodes2, type, uv, &sense),
        "cell uv0");
    printf("21uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    RSS(ref_geom_cell_tuv(ref_geom, nodes2[2], nodes2, type, uv, &sense),
        "cell uv0");
    printf("21uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes2[node]) nodes2[node] = node1;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes3[node]) nodes3[node] = new_node;
    RSS(ref_geom_cell_tuv(ref_geom, nodes3[0], nodes3, type, uv, &sense),
        "cell uv0");
    printf("30uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    RSS(ref_geom_cell_tuv(ref_geom, nodes3[1], nodes3, type, uv, &sense),
        "cell uv0");
    printf("30uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    RSS(ref_geom_cell_tuv(ref_geom, nodes3[2], nodes3, type, uv, &sense),
        "cell uv0");
    printf("30uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes3[node]) nodes3[node] = node0;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes3[node]) nodes3[node] = new_node;
    RSS(ref_geom_cell_tuv(ref_geom, nodes3[0], nodes3, type, uv, &sense),
        "cell uv0");
    printf("31uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    RSS(ref_geom_cell_tuv(ref_geom, nodes3[1], nodes3, type, uv, &sense),
        "cell uv0");
    printf("31uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    RSS(ref_geom_cell_tuv(ref_geom, nodes3[2], nodes3, type, uv, &sense),
        "cell uv0");
    printf("31uv %.18e %.18e\n", uv[0] - u0, uv[1] - v0);
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes3[node]) nodes3[node] = node1;
  }

  return REF_SUCCESS;
}
static REF_STATUS ref_geom_add_between_face_interior(REF_GRID ref_grid,
                                                     REF_INT node0,
                                                     REF_INT node1,
                                                     REF_DBL node1_weight,
                                                     REF_INT new_node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT type, id;
  REF_DBL uv[2], uv0[2], uv1[2];
  REF_INT sense, ncell, cells[2];
  REF_INT nodes2[REF_CELL_MAX_SIZE_PER];
  REF_INT nodes3[REF_CELL_MAX_SIZE_PER];
  REF_DBL r0, r1;
  REF_DBL weight, actual, error, relax, last_error;
  REF_BOOL verbose = REF_FALSE;
  REF_INT i;

  type = REF_GEOM_FACE;
  RSS(ref_geom_unique_id(ref_geom, new_node, type, &id), "unique face id");

  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cells), "list");
  REIS(2, ncell, "expected two tri for box2 nodes");
  RSS(ref_cell_nodes(ref_cell, cells[0], nodes2), "cell nodes");
  RSS(ref_cell_nodes(ref_cell, cells[1], nodes3), "cell nodes");

  RSS(ref_geom_cell_tuv(ref_geom, node0, nodes2, type, uv0, &sense),
      "cell uv0");
  RSS(ref_geom_cell_tuv(ref_geom, node1, nodes2, type, uv1, &sense),
      "cell uv1");

  {
    REF_DBL temp_uv[2], temp_xyz[3];
    temp_uv[0] = (1.0 - node1_weight) * uv0[0] + node1_weight * uv1[0];
    temp_uv[1] = (1.0 - node1_weight) * uv0[1] + node1_weight * uv1[1];
    temp_xyz[0] = ref_node_xyz(ref_node, 0, new_node);
    temp_xyz[1] = ref_node_xyz(ref_node, 1, new_node);
    temp_xyz[2] = ref_node_xyz(ref_node, 2, new_node);
    if (REF_SUCCESS ==
        ref_egads_invert(ref_geom, type, id, temp_xyz, temp_uv)) {
      RSS(ref_geom_add(ref_geom, new_node, type, id, temp_uv), "new geom");
    };
  }

  weight = node1_weight;
  relax = 1.0;
  error = 0.0;
  for (i = 0; i < 10; i++) {
    uv[0] = (1.0 - weight) * uv0[0] + weight * uv1[0];
    uv[1] = (1.0 - weight) * uv0[1] + weight * uv1[1];
    RSS(ref_geom_add(ref_geom, new_node, type, id, uv), "new geom");
    RSS(ref_egads_eval_at(ref_geom, type, id, uv,
                          ref_node_xyz_ptr(ref_node, new_node), NULL),
        "eval");
    RSS(ref_node_ratio(ref_node, node0, new_node, &r0), "get r0");
    RSS(ref_node_ratio(ref_node, node1, new_node, &r1), "get r1");
    if (!ref_math_divisible(r0, (r0 + r1))) break;
    actual = r0 / (r0 + r1);
    last_error = error;
    error = actual - node1_weight;
    if (i > 0 &&
        ((error < 0 && last_error > 0) || (error > 0 && last_error < 0)))
      relax *= 0.5;
    if (verbose)
      printf("target %f actual %f adjust %f errro %f\n", node1_weight, actual,
             weight, error);
    error = MAX(-0.05, MIN(0.05, error));
    weight = weight - relax * error;
    weight = MAX(0.01, MIN(0.99, weight));
    if (ABS(error) < 0.005) break;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_add_between(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_DBL node1_weight, REF_INT new_node) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT type, id;
  REF_DBL param[2], param0[2], param1[2];
  REF_BOOL has_edge_support, supported;
  REF_INT edge_geom;
  REF_INT sense, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_STATUS status;
  REF_INT i, ncell, cells[2];
  REF_INT face_geom, geom0, geom1;
  REF_BOOL support0, support1;
  REF_DBL node0_weight = 1.0 - node1_weight;

  RSS(ref_geom_supported(ref_geom, node0, &support0), "node0 supported");
  RSS(ref_geom_supported(ref_geom, node1, &support1), "node1 supported");
  if (!support0 || !support1) {
    return REF_SUCCESS;
  }

  /* insert edge geom on edge cell if present */
  nodes[0] = node0;
  nodes[1] = node1;
  ref_cell = ref_grid_edg(ref_grid);
  status = ref_cell_with(ref_cell, nodes, &cell);
  if (REF_NOT_FOUND == status) {
    has_edge_support = REF_FALSE;
    edge_geom = REF_EMPTY;
  } else {
    RSS(status, "search for edg");
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "get id");
    id = nodes[ref_cell_node_per(ref_cell)];
    type = REF_GEOM_EDGE;
    RSS(ref_geom_cell_tuv(ref_geom, node0, nodes, type, param0, &sense),
        "cell uv");
    RSS(ref_geom_cell_tuv(ref_geom, node1, nodes, type, param1, &sense),
        "cell uv");
    param[0] = node0_weight * param0[0] + node1_weight * param1[0];
    if (ref_geom_model_loaded(ref_geom))
      RSB(ref_egads_inverse_eval(ref_geom, type, id,
                                 ref_node_xyz_ptr(ref_node, new_node), param),
          "inv eval edge", ref_geom_tec(ref_grid, "ref_geom_split_edge.tec"));
    /* enforce bounding box and use midpoint as full-back */
    if (param[0] < MIN(param0[0], param1[0]) ||
        MAX(param0[0], param1[0]) < param[0])
      param[0] = node0_weight * param0[0] + node1_weight * param1[0];
    if (ref_geom_model_loaded(ref_geom)) { /* check weight and distance ratio */
      REF_DBL xyz[3], dx0[3], dx1[3], d0, d1, total, actual_weight;
      REF_DBL mid_t, mid_weight;
      REF_INT ii;
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_EDGE, id, param, xyz, NULL),
          "eval");
      for (ii = 0; ii < 3; ii++)
        dx0[ii] = ref_node_xyz(ref_node, ii, node0) - xyz[ii];
      for (ii = 0; ii < 3; ii++)
        dx1[ii] = ref_node_xyz(ref_node, ii, node1) - xyz[ii];
      d0 = sqrt(ref_math_dot(dx0, dx0));
      d1 = sqrt(ref_math_dot(dx1, dx1));
      total = d0 + d1;
      actual_weight = -1.0;
      if (ref_math_divisible(d0, total)) {
        actual_weight = d0 / total;
      }
      mid_t = node0_weight * param0[0] + node1_weight * param1[0];
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_EDGE, id, &mid_t, xyz, NULL),
          "eval");
      for (ii = 0; ii < 3; ii++)
        dx0[ii] = ref_node_xyz(ref_node, ii, node0) - xyz[ii];
      for (ii = 0; ii < 3; ii++)
        dx1[ii] = ref_node_xyz(ref_node, ii, node1) - xyz[ii];
      d0 = sqrt(ref_math_dot(dx0, dx0));
      d1 = sqrt(ref_math_dot(dx1, dx1));
      total = d0 + d1;
      mid_weight = -1.0;
      if (ref_math_divisible(d0, total)) {
        mid_weight = d0 / total;
      }
      if (ABS(mid_weight - node1_weight) < ABS(actual_weight - node1_weight)) {
        /* printf("request %f actual %f mid %f\n",
                  node1_weight,actual_weight,mid_weight); */
        param[0] = mid_t;
      }
    }
    if (param[0] < MIN(param0[0], param1[0]) ||
        MAX(param0[0], param1[0]) < param[0])
      param[0] = node0_weight * param0[0] + node1_weight * param1[0];
    RSS(ref_geom_add(ref_geom, new_node, type, id, param), "new geom");
    has_edge_support = REF_TRUE;
    RSS(ref_geom_find(ref_geom, new_node, type, id, &edge_geom),
        "find the new edge for later face uv evaluation");
  }

  /* insert face between */
  ref_cell = ref_grid_tri(ref_grid);
  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cells), "list");
  if (0 == ncell) { /* volume edge */
    return REF_SUCCESS;
  }
  if (ref_geom_manifold(ref_geom)) {
    REIB(2, ncell, "expected two tri for between", {
      ref_geom_tattle(ref_geom, node0);
      ref_geom_tattle(ref_geom, node1);
      ref_node_location(ref_node, node0);
      ref_node_location(ref_node, node1);
    });
  }
  for (i = 0; i < ncell; i++) {
    cell = cells[i];
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "get id");

    RSS(ref_geom_tri_supported(ref_geom, nodes, &supported), "tri support");
    if (!supported) continue; /* no geom support, skip */

    id = nodes[ref_cell_node_per(ref_cell)];
    type = REF_GEOM_FACE;
    RSS(ref_geom_cell_tuv(ref_geom, node0, nodes, type, param0, &sense),
        "cell uv");
    RSS(ref_geom_cell_tuv(ref_geom, node1, nodes, type, param1, &sense),
        "cell uv");
    param[0] = node0_weight * param0[0] + node1_weight * param1[0];
    param[1] = node0_weight * param0[1] + node1_weight * param1[1];

    RSS(ref_geom_add(ref_geom, new_node, type, id, param), "new geom");
    RSS(ref_geom_find(ref_geom, new_node, type, id, &face_geom),
        "new face geom");

    RSS(ref_geom_find(ref_geom, node0, type, id, &geom0), "face geom");
    RSS(ref_geom_find(ref_geom, node1, type, id, &geom1), "face geom");
    if (0 != ref_geom_jump(ref_geom, geom0) &&
        0 != ref_geom_jump(ref_geom, geom1) &&
        ref_geom_jump(ref_geom, geom0) == ref_geom_jump(ref_geom, geom1)) {
      ref_geom_jump(ref_geom, face_geom) = ref_geom_jump(ref_geom, geom0);
    }

    /* if there is an edge between, set the face uv based on edge t */
    if (ref_geom_model_loaded(ref_geom) && has_edge_support) {
      REF_INT faceid, edgeid;
      REF_DBL t;
      edgeid = ref_geom_id(ref_geom, edge_geom);
      faceid = ref_geom_id(ref_geom, face_geom);
      t = ref_geom_param(ref_geom, 0, edge_geom);
      sense = 0;
      if (0 != ref_geom_jump(ref_geom, face_geom)) sense = 1;
      RSS(ref_egads_edge_face_uv(ref_geom, edgeid, faceid, sense, t, param),
          "edge t param");
      ref_geom_param(ref_geom, 0, face_geom) = param[0];
      ref_geom_param(ref_geom, 1, face_geom) = param[1];
    }
  }

  if (ref_geom_model_loaded(ref_geom) && !has_edge_support) {
    RSS(ref_geom_add_between_face_interior(ref_grid, node0, node1, node1_weight,
                                           new_node),
        "position new node in uv");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_geom_add_constrain_midnode(REF_GRID ref_grid,
                                                 REF_INT node0, REF_INT node1,
                                                 REF_DBL node1_weight,
                                                 REF_INT new_node) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT type, id;
  REF_DBL param[2], param0[2], param1[2];
  REF_BOOL has_edge_support, supported;
  REF_INT sense, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_STATUS status;
  REF_INT i, ncell, cells[2];
  REF_INT face_geom, geom0, geom1;
  REF_BOOL support0, support1;
  REF_DBL node0_weight = 1.0 - node1_weight;

  RSS(ref_geom_supported(ref_geom, node0, &support0), "node0 supported");
  RSS(ref_geom_supported(ref_geom, node1, &support1), "node1 supported");
  if (!support0 || !support1) {
    return REF_SUCCESS;
  }

  /* insert edge geom on edge cell if present */
  nodes[0] = node0;
  nodes[1] = node1;
  ref_cell = ref_grid_edg(ref_grid);
  status = ref_cell_with(ref_cell, nodes, &cell);
  if (REF_NOT_FOUND == status) {
    has_edge_support = REF_FALSE;
  } else {
    RSS(status, "search for edg");
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "get id");
    id = nodes[ref_cell_node_per(ref_cell)];
    type = REF_GEOM_EDGE;
    RSS(ref_geom_cell_tuv(ref_geom, node0, nodes, type, param0, &sense),
        "cell uv");
    RSS(ref_geom_cell_tuv(ref_geom, node1, nodes, type, param1, &sense),
        "cell uv");
    param[0] = node0_weight * param0[0] + node1_weight * param1[0];
    RSS(ref_egads_eval_at(ref_geom, type, id, param,
                          ref_node_xyz_ptr(ref_node, new_node), NULL),
        "eval");
    RSS(ref_geom_add(ref_geom, new_node, type, id, param), "new geom");
    has_edge_support = REF_TRUE;
  }

  /* insert face between */
  ref_cell = ref_grid_tri(ref_grid);
  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cells), "list");
  if (0 == ncell) { /* volume edge */
    return REF_SUCCESS;
  }
  if (ref_geom_manifold(ref_geom)) {
    REIB(2, ncell, "expected two tri for between", {
      ref_geom_tattle(ref_geom, node0);
      ref_geom_tattle(ref_geom, node1);
      ref_node_location(ref_node, node0);
      ref_node_location(ref_node, node1);
    });
  }
  for (i = 0; i < ncell; i++) {
    cell = cells[i];
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "get id");

    RSS(ref_geom_tri_supported(ref_geom, nodes, &supported), "tri support");
    if (!supported) continue; /* no geom support, skip */

    id = nodes[ref_cell_node_per(ref_cell)];
    type = REF_GEOM_FACE;
    RSS(ref_geom_cell_tuv(ref_geom, node0, nodes, type, param0, &sense),
        "cell uv");
    RSS(ref_geom_cell_tuv(ref_geom, node1, nodes, type, param1, &sense),
        "cell uv");
    param[0] = node0_weight * param0[0] + node1_weight * param1[0];
    param[1] = node0_weight * param0[1] + node1_weight * param1[1];

    RSS(ref_geom_add(ref_geom, new_node, type, id, param), "new geom");
    RSS(ref_geom_find(ref_geom, new_node, type, id, &face_geom),
        "new face geom");

    RSS(ref_geom_find(ref_geom, node0, type, id, &geom0), "face geom");
    RSS(ref_geom_find(ref_geom, node1, type, id, &geom1), "face geom");
    if (0 != ref_geom_jump(ref_geom, geom0) &&
        0 != ref_geom_jump(ref_geom, geom1) &&
        ref_geom_jump(ref_geom, geom0) == ref_geom_jump(ref_geom, geom1)) {
      ref_geom_jump(ref_geom, face_geom) = ref_geom_jump(ref_geom, geom0);
    }
    if (!has_edge_support) {
      RSS(ref_egads_eval_at(ref_geom, type, id, param,
                            ref_node_xyz_ptr(ref_node, new_node), NULL),
          "eval");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_support_between(REF_GRID ref_grid, REF_INT node0,
                                    REF_INT node1, REF_BOOL *needs_support) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT item0, item1;
  REF_INT geom0, geom1;
  REF_INT type, id;
  REF_BOOL has_id;

  *needs_support = REF_FALSE;
  /* assume face check is sufficient */
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node0, item0, geom0)
      each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node1, item1,
                                      geom1) {
    type = REF_GEOM_FACE;
    if (ref_geom_type(ref_geom, geom0) == type &&
        ref_geom_type(ref_geom, geom1) == type &&
        ref_geom_id(ref_geom, geom0) == ref_geom_id(ref_geom, geom1)) {
      id = ref_geom_id(ref_geom, geom0);
      RSS(ref_cell_side_has_id(ref_grid_tri(ref_grid), node0, node1, id,
                               &has_id),
          "has edge id");
      if (has_id) {
        *needs_support = REF_TRUE;
        return REF_SUCCESS;
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tri_uv_bounding_box(REF_GRID ref_grid, REF_INT node,
                                        REF_DBL *uv_min, REF_DBL *uv_max) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_node, id, iuv;
  REF_DBL uv[2];
  REF_INT sense;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  /* get face id and initialize min and max */
  RSS(ref_geom_unique_id(ref_geom, node, REF_GEOM_FACE, &id), "id");
  RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, id, uv_min), "uv_min");
  RSS(ref_geom_tuv(ref_geom, node, REF_GEOM_FACE, id, uv_max), "uv_max");

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_geom_cell_tuv(ref_geom, nodes[cell_node], nodes, REF_GEOM_FACE,
                            uv, &sense),
          "cell uv");
      for (iuv = 0; iuv < 2; iuv++) uv_min[iuv] = MIN(uv_min[iuv], uv[iuv]);
      for (iuv = 0; iuv < 2; iuv++) uv_max[iuv] = MAX(uv_max[iuv], uv[iuv]);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tri_uv_bounding_box2(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_DBL *uv_min,
                                         REF_DBL *uv_max) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell, cell_node, iuv;
  REF_DBL uv[2];
  REF_INT i, ncell, cells[2];
  REF_INT sense;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  RSS(ref_cell_list_with2(ref_cell, node0, node1, 2, &ncell, cells), "list");
  REIS(2, ncell, "expected two tri for box2 nodes");

  for (iuv = 0; iuv < 2; iuv++) uv_min[iuv] = REF_DBL_MAX;
  for (iuv = 0; iuv < 2; iuv++) uv_max[iuv] = REF_DBL_MIN;
  for (i = 0; i < ncell; i++) {
    cell = cells[i];
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_geom_cell_tuv(ref_geom, nodes[cell_node], nodes, REF_GEOM_FACE,
                            uv, &sense),
          "cell uv");
      for (iuv = 0; iuv < 2; iuv++) uv_min[iuv] = MIN(uv_min[iuv], uv[iuv]);
      for (iuv = 0; iuv < 2; iuv++) uv_max[iuv] = MAX(uv_max[iuv], uv[iuv]);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_constrain_all(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_geom_constrain(ref_grid, node), "constrain node");
  }
  return REF_SUCCESS;
}

REF_STATUS ref_geom_constrain(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT item, geom;
  REF_BOOL have_geom_node;
  REF_BOOL have_geom_edge;
  REF_BOOL have_geom_face;
  REF_INT node_geom;
  REF_INT edge_geom;
  REF_INT face_geom;
  REF_DBL xyz[3];

  /* no geom, do nothing */
  if (ref_adj_empty(ref_adj, node)) return REF_SUCCESS;

  if (ref_geom_meshlinked(ref_geom)) {
    RSS(ref_meshlink_constrain(ref_grid, node), "meshlink");
    return REF_SUCCESS;
  }

  have_geom_node = REF_FALSE;
  node_geom = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_adj, node, item, geom) {
    if (REF_GEOM_NODE == ref_geom_type(ref_geom, geom)) {
      have_geom_node = REF_TRUE;
      node_geom = geom;
      break;
    }
  }

  if (have_geom_node) { /* update T of edges? update UV of (degen) faces? */
    RSS(ref_egads_eval(ref_geom, node_geom, xyz, NULL), "eval edge");
    node = ref_geom_node(ref_geom, node_geom);
    ref_node_xyz(ref_node, 0, node) = xyz[0];
    ref_node_xyz(ref_node, 1, node) = xyz[1];
    ref_node_xyz(ref_node, 2, node) = xyz[2];
    return REF_SUCCESS;
  }

  have_geom_edge = REF_FALSE;
  edge_geom = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_adj, node, item, geom) {
    if (REF_GEOM_EDGE == ref_geom_type(ref_geom, geom)) {
      have_geom_edge = REF_TRUE;
      edge_geom = geom;
      break;
    }
  }

  /* edge geom, evaluate edge and update face uv */
  if (have_geom_edge) {
    RSS(ref_egads_eval(ref_geom, edge_geom, xyz, NULL), "eval edge");
    node = ref_geom_node(ref_geom, edge_geom);
    ref_node_xyz(ref_node, 0, node) = xyz[0];
    ref_node_xyz(ref_node, 1, node) = xyz[1];
    ref_node_xyz(ref_node, 2, node) = xyz[2];
    RSS(ref_geom_eval_edge_face_uv(ref_grid, edge_geom), "resol edge uv");
    return REF_SUCCESS;
  }

  /* look for face geom */
  have_geom_face = REF_FALSE;
  face_geom = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_adj, node, item, geom) {
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom)) {
      have_geom_face = REF_TRUE;
      face_geom = geom;
      break;
    }
  }

  /* face geom, evaluate on face uv */
  if (have_geom_face) {
    RSS(ref_egads_eval(ref_geom, face_geom, xyz, NULL), "eval face");
    node = ref_geom_node(ref_geom, face_geom);
    ref_node_xyz(ref_node, 0, node) = xyz[0];
    ref_node_xyz(ref_node, 1, node) = xyz[1];
    ref_node_xyz(ref_node, 2, node) = xyz[2];
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_radian_request(REF_GEOM ref_geom, REF_INT geom,
                                   REF_DBL *delta_radian) {
  REF_INT node, item, face_geom, face;
  REF_DBL segments, face_segments;
  REF_DBL face_set = -998.0;
  REF_DBL face_active = 0.01;
  REF_BOOL turn_off, use_face;
  use_face = REF_FALSE;
  turn_off = REF_FALSE;

  segments = 0.0;

  if (NULL != (ref_geom)->face_seg_per_rad) {
    node = ref_geom_node(ref_geom, geom);
    each_ref_geom_having_node(ref_geom, node, item, face_geom) {
      if (REF_GEOM_FACE == ref_geom_type(ref_geom, face_geom)) {
        face = ref_geom_id(ref_geom, face_geom) - 1;
        face_segments = ref_geom->face_seg_per_rad[face];
        if (face_segments < face_set) continue;
        if (face_segments > face_active) {
          segments = MAX(segments, face_segments);
          use_face = REF_TRUE;
        } else {
          turn_off = REF_TRUE;
          break;
        }
      }
    }
  }

  if (turn_off) {
    *delta_radian = 100.0;
    return REF_SUCCESS;
  }

  if (!use_face) {
    segments = ref_geom_segments_per_radian_of_curvature(ref_geom);
  }

  if (segments > 0.1 && ref_math_divisible(1.0, segments)) {
    *delta_radian = 1.0 / segments;
  } else {
    *delta_radian = 100.0;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_uv_rsn(REF_DBL *uv, REF_DBL *r, REF_DBL *s, REF_DBL *n,
                           REF_DBL *drsduv) {
  REF_INT i;
  REF_DBL dot;
  REF_DBL len;

  for (i = 0; i < 3; i++) r[i] = uv[i];
  drsduv[0] = 1.0;
  drsduv[1] = 0.0;
  for (i = 0; i < 3; i++) s[i] = uv[i + 3];
  drsduv[2] = 0.0;
  drsduv[3] = 1.0;
  len = sqrt(ref_math_dot(r, r));
  drsduv[0] /= len;
  drsduv[1] /= len;
  RAISE(ref_math_normalize(r));
  len = sqrt(ref_math_dot(s, s));
  drsduv[2] /= len;
  drsduv[3] /= len;
  RAISE(ref_math_normalize(s));

  dot = ref_math_dot(r, s);
  for (i = 0; i < 3; i++) s[i] -= dot * r[i];
  drsduv[2] -= dot * drsduv[0];
  drsduv[3] -= dot * drsduv[1];

  len = sqrt(ref_math_dot(s, s));
  drsduv[2] /= len;
  drsduv[3] /= len;
  RAISE(ref_math_normalize(s));

  ref_math_cross_product(r, s, n);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_face_rsn(REF_GEOM ref_geom, REF_INT faceid, REF_DBL *uv,
                             REF_DBL *r, REF_DBL *s, REF_DBL *n) {
  REF_DBL xyz[3];
  REF_DBL dxyz_dtuv[15];
  REF_DBL drsduv[4];
  RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, faceid, uv, xyz, dxyz_dtuv),
      "eval");
  RAISE(ref_geom_uv_rsn(dxyz_dtuv, r, s, n, drsduv));
  return REF_SUCCESS;
}

REF_STATUS ref_geom_tri_centroid(REF_GRID ref_grid, REF_INT *nodes,
                                 REF_DBL *uv) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell_node;
  REF_DBL node_uv[2];
  REF_INT sens;
  uv[0] = 0.0;
  uv[1] = 0.0;
  each_ref_cell_cell_node(ref_cell, cell_node) {
    RSB(ref_geom_cell_tuv(ref_geom, nodes[cell_node], nodes, REF_GEOM_FACE,
                          node_uv, &sens),
        "cell node uv",
        { ref_geom_tec(ref_grid, "ref_geom_tri_centroid_error.tec"); });
    uv[0] += (1.0 / 3.0) * node_uv[0];
    uv[1] += (1.0 / 3.0) * node_uv[1];
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tri_norm_deviation(REF_GRID ref_grid, REF_INT *nodes,
                                       REF_DBL *dot_product) {
  REF_DBL uv[2];
  REF_DBL tri_normal[3];
  REF_DBL r[3], s[3], n[3], area_sign;
  REF_INT id;
  REF_STATUS status;
  *dot_product = -2.0;

  if (ref_geom_meshlinked(ref_grid_geom(ref_grid))) {
    RSS(ref_meshlink_tri_norm_deviation(ref_grid, nodes, dot_product),
        "meshlink");
    return REF_SUCCESS;
  }

  id = nodes[ref_cell_node_per(ref_grid_tri(ref_grid))];
  RSS(ref_node_tri_normal(ref_grid_node(ref_grid), nodes, tri_normal),
      "tri normal");
  /* collapse attempts could create zero area, reject the step with -2.0 */
  status = ref_math_normalize(tri_normal);
  if (REF_DIV_ZERO == status) return REF_SUCCESS;
  RSS(status, "normalize");

  RSS(ref_geom_tri_centroid(ref_grid, nodes, uv), "tri cent");
  RAISE(ref_geom_face_rsn(ref_grid_geom(ref_grid), id, uv, r, s, n));
  RSS(ref_geom_uv_area_sign(ref_grid, id, &area_sign), "a sign");

  *dot_product = area_sign * ref_math_dot(n, tri_normal);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_crease(REF_GRID ref_grid, REF_INT node, REF_DBL *dot_prod) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item0, item1, cell0, cell1, id;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL uv[2];
  REF_DBL r[3], s[3];
  REF_DBL n0[3], area_sign0;
  REF_DBL n1[3], area_sign1;

  *dot_prod = 1.0;

  if (!ref_node_valid(ref_node, node)) {
    return REF_SUCCESS;
  }

  each_ref_cell_having_node(ref_cell, node, item0, cell0) {
    RSS(ref_cell_nodes(ref_grid_tri(ref_grid), cell0, nodes),
        "tri list for edge");
    RSS(ref_geom_tri_centroid(ref_grid, nodes, uv), "tri cent");
    id = nodes[ref_cell_node_per(ref_cell)];
    RSS(ref_geom_face_rsn(ref_grid_geom(ref_grid), id, uv, r, s, n0), "rsn");
    RSS(ref_geom_uv_area_sign(ref_grid, id, &area_sign0), "a sign");
    each_ref_cell_having_node(ref_cell, node, item1, cell1) {
      RSS(ref_cell_nodes(ref_grid_tri(ref_grid), cell1, nodes),
          "tri list for edge");
      RSS(ref_geom_tri_centroid(ref_grid, nodes, uv), "tri cent");
      id = nodes[ref_cell_node_per(ref_cell)];
      RSS(ref_geom_face_rsn(ref_grid_geom(ref_grid), id, uv, r, s, n1), "rsn");
      RSS(ref_geom_uv_area_sign(ref_grid, id, &area_sign1), "a sign");

      *dot_prod =
          MIN(*dot_prod, area_sign0 * area_sign1 * ref_math_dot(n0, n1));
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_max_gap(REF_GRID ref_grid, REF_DBL *max_gap) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT geom;
  REF_INT node;
  REF_DBL xyz[3];
  REF_DBL dist, max, global_max;

  *max_gap = 0.0;

  if (!ref_geom_model_loaded(ref_geom)) return REF_SUCCESS;

  max = 0.0;
  each_ref_geom_node(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_egads_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    max = MAX(max, dist);
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max node");
  max = global_max;
  *max_gap = MAX(*max_gap, max);

  max = 0.0;
  each_ref_geom_edge(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_egads_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    max = MAX(max, dist);
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max edge");
  max = global_max;
  *max_gap = MAX(*max_gap, max);

  max = 0.0;
  each_ref_geom_face(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_egads_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    max = MAX(max, dist);
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max face");
  max = global_max;
  *max_gap = MAX(*max_gap, max);
  RSS(ref_mpi_bcast(ref_mpi, max_gap, 1, REF_DBL_TYPE), "mpi bcast gap");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_verify_param(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT geom;
  REF_INT node;
  REF_DBL xyz[3];
  REF_DBL dist, max, max_node, max_edge, global_max;
  REF_BOOL node_constraint, edge_constraint;

  if (!ref_geom_model_loaded(ref_geom)) return REF_SUCCESS;

  max = 0.0;
  each_ref_geom_node(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_egads_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    max = MAX(max, dist);
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max node");
  max = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo node max eval dist %e\n", max);

  max = 0.0;
  max_node = 0.0;
  each_ref_geom_edge(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_egads_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &node_constraint), "n");
    if (node_constraint) {
      max_node = MAX(max_node, dist);
    } else {
      max = MAX(max, dist);
    }
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max edge");
  max = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo edge max eval dist %e\n", max);
  RSS(ref_mpi_max(ref_mpi, &max_node, &global_max, REF_DBL_TYPE),
      "mpi max edge");
  max_node = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo edge node tol %e\n", max_node);

  max = 0.0;
  max_node = 0.0;
  max_edge = 0.0;
  each_ref_geom_face(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) continue;
    RSS(ref_egads_eval(ref_geom, geom, xyz, NULL), "eval xyz");
    dist = sqrt(pow(xyz[0] - ref_node_xyz(ref_node, 0, node), 2) +
                pow(xyz[1] - ref_node_xyz(ref_node, 1, node), 2) +
                pow(xyz[2] - ref_node_xyz(ref_node, 2, node), 2));
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &node_constraint), "n");
    RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &edge_constraint), "n");
    if (node_constraint) {
      max_node = MAX(max_node, dist);
    } else {
      if (edge_constraint) {
        max_edge = MAX(max_edge, dist);
      } else {
        max = MAX(max, dist);
      }
    }
  }
  RSS(ref_mpi_max(ref_mpi, &max, &global_max, REF_DBL_TYPE), "mpi max face");
  max = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo face max eval dist %e\n", max);
  RSS(ref_mpi_max(ref_mpi, &max_edge, &global_max, REF_DBL_TYPE),
      "mpi max edge");
  max_edge = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo face edge tol %e\n", max_edge);
  RSS(ref_mpi_max(ref_mpi, &max_node, &global_max, REF_DBL_TYPE),
      "mpi max edge");
  max_node = global_max;
  if (ref_grid_once(ref_grid)) printf("CAD topo face node tol %e\n", max_node);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_verify_topo(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT item, geom;
  REF_BOOL geom_node, geom_edge, geom_face;
  REF_BOOL no_face, no_edge;
  REF_BOOL found_one;
  REF_BOOL found_too_many;
  REF_INT cell, ncell, cell_list[2];

  for (node = 0; node < ref_node_max(ref_node); node++) {
    if (ref_node_valid(ref_node, node)) {
      RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &geom_node), "node");
      RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &geom_edge), "edge");
      RSS(ref_geom_is_a(ref_geom, node, REF_GEOM_FACE, &geom_face), "face");
      no_face = ref_cell_node_empty(ref_grid_tri(ref_grid), node) &&
                ref_cell_node_empty(ref_grid_tr2(ref_grid), node) &&
                ref_cell_node_empty(ref_grid_tr3(ref_grid), node) &&
                ref_cell_node_empty(ref_grid_qua(ref_grid), node);
      no_edge = ref_cell_node_empty(ref_grid_edg(ref_grid), node) &&
                ref_cell_node_empty(ref_grid_ed2(ref_grid), node) &&
                ref_cell_node_empty(ref_grid_ed3(ref_grid), node);
      if (geom_node) {
        if (no_edge && ref_node_owned(ref_node, node)) {
          THROW("geom node missing edge");
        }
        if (no_face && ref_node_owned(ref_node, node)) {
          THROW("geom node missing tri or qua");
        }
      }
      if (geom_edge) {
        if (no_edge && ref_node_owned(ref_node, node)) {
          RSS(ref_grid_tattle(ref_grid, node), "tattle");
          RSS(ref_geom_tec_para_shard(ref_grid, "ref_geom_topo_error"),
              "geom tec");
          THROW("geom edge missing edge");
        }
        if (no_face && ref_node_owned(ref_node, node)) {
          RSS(ref_grid_tattle(ref_grid, node), "tattle");
          RSS(ref_geom_tec_para_shard(ref_grid, "ref_geom_topo_error"),
              "geom tec");
          THROW("geom edge missing tri or qua");
        }
      }
      if (geom_face) {
        if (no_face && ref_node_owned(ref_node, node)) {
          printf("no face for geom\n");
          RSS(ref_grid_tattle(ref_grid, node), "tattle");
          RSS(ref_geom_tec_para_shard(ref_grid, "ref_geom_topo_error"),
              "geom tec");
          THROW("geom face missing tri or qua");
        }
      }
      if (!no_edge) {
        if (!geom_edge) {
          printf("no geom for edge\n");
          RSS(ref_grid_tattle(ref_grid, node), "tattle");
          RSS(ref_geom_tec_para_shard(ref_grid, "ref_geom_topo_error"),
              "geom tec");
          THROW("geom edge missing for edg");
        }
      }
      if (!no_face) {
        if (!geom_face && !ref_geom_meshlinked(ref_geom)) {
          printf("no geom for face\n");
          RSS(ref_grid_tattle(ref_grid, node), "tattle");
          RSS(ref_geom_tec_para_shard(ref_grid, "ref_geom_topo_error"),
              "geom tec");
          THROW("geom face missing tri or qua");
        }
      }
      if (geom_edge && !geom_node) {
        found_one = REF_FALSE;
        found_too_many = REF_FALSE;
        each_ref_geom_having_node(ref_geom, node, item, geom) {
          if (REF_GEOM_EDGE == ref_geom_type(ref_geom, geom)) {
            if (found_one) found_too_many = REF_TRUE;
            found_one = REF_TRUE;
          }
        }
        if (!found_one || found_too_many) {
          if (!found_one) printf("none found\n");
          if (found_too_many) printf("found too many\n");
          RSS(ref_grid_tattle(ref_grid, node), "tatt");
          RSS(ref_geom_tec_para_shard(ref_grid, "ref_geom_topo_error"),
              "geom tec");
          THROW("multiple geom edge away from geom node");
        }
      }
      if (geom_face && !geom_edge) {
        found_one = REF_FALSE;
        found_too_many = REF_FALSE;
        each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item,
                                        geom) {
          if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom)) {
            if (found_one) found_too_many = REF_TRUE;
            found_one = REF_TRUE;
          }
        }
        if (!found_one || found_too_many) {
          if (!found_one) printf("none found\n");
          if (found_too_many) printf("found too many\n");
          RSS(ref_grid_tattle(ref_grid, node), "tattle");
          RSS(ref_geom_tec_para_shard(ref_grid, "ref_geom_topo_error"),
              "geom tec");
          THROW("multiple geom face away from geom edge");
        }
      }
    } else {
      if (!ref_adj_empty(ref_geom_adj(ref_geom), node))
        THROW("invalid node has geom");
    }
  }

  ref_cell = ref_grid_edg(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    RSS(ref_cell_list_with2(ref_cell, ref_cell_c2n(ref_cell, 0, cell),
                            ref_cell_c2n(ref_cell, 1, cell), 2, &ncell,
                            cell_list),
        "edge list for edge");
    if (2 == ncell) {
      printf("error: two edg found with same nodes\n");
      printf("edg %d n %d %d id %d\n", cell_list[0],
             ref_cell_c2n(ref_cell, 0, cell_list[0]),
             ref_cell_c2n(ref_cell, 1, cell_list[0]),
             ref_cell_c2n(ref_cell, 2, cell_list[0]));
      printf("edg %d n %d %d id %d\n", cell_list[1],
             ref_cell_c2n(ref_cell, 0, cell_list[1]),
             ref_cell_c2n(ref_cell, 1, cell_list[1]),
             ref_cell_c2n(ref_cell, 2, cell_list[1]));
      RSS(ref_grid_tattle(ref_grid, ref_cell_c2n(ref_cell, 0, cell)), "tattle");
      RSS(ref_grid_tattle(ref_grid, ref_cell_c2n(ref_cell, 1, cell)), "tattle");
      RSS(ref_geom_tec_para_shard(ref_grid, "ref_geom_topo_error"), "geom tec");
    }
    REIS(1, ncell, "expect only one edge cell for two nodes");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_report_topo_at(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell;
  REF_INT geom_item, geom;
  REF_INT cell_item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT have_match;

  each_ref_geom_having_node(ref_geom, node, geom_item, geom) {
    if (REF_GEOM_EDGE == ref_geom_type(ref_geom, geom)) {
      have_match = REF_FALSE;
      ref_cell = ref_grid_edg(ref_grid);
      each_ref_cell_having_node(ref_cell, node, cell_item, cell) {
        RSS(ref_cell_nodes(ref_cell, cell, nodes), "edg nodes");
        if (ref_geom_id(ref_geom, geom) == nodes[ref_cell_id_index(ref_cell)]) {
          have_match = REF_TRUE;
        }
      }
      if (!have_match) {
        RSS(ref_node_location(ref_node, node), "loc");
        RSS(ref_geom_tattle(ref_geom, node), "geom tatt");
        RSS(ref_cell_tattle_about(ref_cell, node), "cell tatt");
      }
    }
    if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom)) {
      have_match = REF_FALSE;
      ref_cell = ref_grid_tri(ref_grid);
      each_ref_cell_having_node(ref_cell, node, cell_item, cell) {
        RSS(ref_cell_nodes(ref_cell, cell, nodes), "tri nodes");
        if (ref_geom_id(ref_geom, geom) == nodes[ref_cell_id_index(ref_cell)]) {
          have_match = REF_TRUE;
        }
      }
      if (!have_match) {
        RSS(ref_node_location(ref_node, node), "loc");
        RSS(ref_geom_tattle(ref_geom, node), "geom tatt");
        RSS(ref_cell_tattle_about(ref_cell, node), "cell tatt");
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tetgen_volume(REF_GRID ref_grid, const char *project,
                                  const char *options) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  char filename[896];
  char command[1024];
  FILE *file;
  REF_INT nnode, ndim, attr, mark;
  REF_GLOB global;
  REF_INT ntri, ntet, node_per, id;
  REF_INT node, nnode_surface, item, new_node;
  REF_DBL xyz[3], dist;
  REF_INT cell, new_cell, nodes[REF_CELL_MAX_SIZE_PER];
  int system_status;
  REF_BOOL delete_temp_files = REF_TRUE;
  REF_BOOL problem;
  REF_BOOL *position;

  printf("%d surface nodes %d triangles\n", ref_node_n(ref_node),
         ref_cell_n(ref_grid_tri(ref_grid)));

  snprintf(filename, 896, "%s-tetgen.poly", project);
  RSS(ref_export_by_extension(ref_grid, filename), "poly");

  printf("  The 'S1000' argument can be added to tetgen options\n");
  printf("    to cap inserted nodes at 1000 and reduce run time.\n");
  printf("  The 'q20/10' arguments (radius-edge-ratio/dihedral-angle)\n");
  printf("    can be increased for faster initial volume refinement.\n");
  printf("  The 'O7/7' arguments (optimization iterations/operation)\n");
  printf("    can be decreased for faster mesh optimization.\n");
  printf("  See 'ref bootstrap -h' for '--mesher-options' description\n");
  printf("    and the TetGen user manual for details.\n");

  if (NULL == options) {
    snprintf(command, 1024,
             "tetgen -pMYq20/10O7/7zVT1e-12 %s < /dev/null > %s-tetgen.txt",
             filename, project);
  } else {
    snprintf(command, 1024, "tetgen %s %s < /dev/null > %s-tetgen.txt", options,
             filename, project);
  }
  printf("%s\n", command);
  fflush(stdout);
  REIS(0, sleep(2), "sleep failed");
  system_status = system(command);
  REIS(0, sleep(2), "sleep failed");
  REIB(0, system_status, "tetgen failed", {
    printf("tec360 ref_geom_test_tetgen_geom.tec\n");
    ref_geom_tec(ref_grid, "ref_geom_test_tetgen_geom.tec");
    printf("tec360 ref_geom_test_tetgen_surf.tec\n");
    ref_export_tec_surf(ref_grid, "ref_geom_test_tetgen_surf.tec");
  });

  snprintf(filename, 896, "%s-tetgen.1.node", project);
  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  REIS(1, fscanf(file, "%d", &nnode), "node header nnode");
  REIS(1, fscanf(file, "%d", &ndim), "node header ndim");
  REIS(3, ndim, "not 3D");
  REIS(1, fscanf(file, "%d", &attr), "node header attr");
  REIS(0, attr, "nodes have attribute 3D");
  REIS(1, fscanf(file, "%d", &mark), "node header mark");
  REIS(0, mark, "nodes have mark");

  /* verify surface nodes */
  nnode_surface = ref_node_n(ref_node);

  printf("%d interior nodes\n", nnode - nnode_surface);

  for (node = 0; node < nnode_surface; node++) {
    REIS(1, fscanf(file, "%d", &item), "node item");
    RES(node, item, "node index");
    RES(1, fscanf(file, "%lf", &(xyz[0])), "x");
    RES(1, fscanf(file, "%lf", &(xyz[1])), "y");
    RES(1, fscanf(file, "%lf", &(xyz[2])), "z");
    dist = sqrt((xyz[0] - ref_node_xyz(ref_node, 0, node)) *
                    (xyz[0] - ref_node_xyz(ref_node, 0, node)) +
                (xyz[1] - ref_node_xyz(ref_node, 1, node)) *
                    (xyz[1] - ref_node_xyz(ref_node, 1, node)) +
                (xyz[2] - ref_node_xyz(ref_node, 2, node)) *
                    (xyz[2] - ref_node_xyz(ref_node, 2, node)));
    if (dist > 1.0e-12) {
      printf("node %d off by %e\n", node, dist);
      THROW("tetgen moved node");
    }
  }

  /* interior nodes */
  for (node = nnode_surface; node < nnode; node++) {
    REIS(1, fscanf(file, "%d", &item), "node item");
    REIS(node, item, "file node index");
    RSS(ref_node_next_global(ref_node, &global), "next global");
    REIS(node, global, "global node index");
    RSS(ref_node_add(ref_node, global, &new_node), "new_node");
    RES(node, new_node, "node index");
    RES(1, fscanf(file, "%lf", &(xyz[0])), "x");
    RES(1, fscanf(file, "%lf", &(xyz[1])), "y");
    RES(1, fscanf(file, "%lf", &(xyz[2])), "z");
    ref_node_xyz(ref_node, 0, new_node) = xyz[0];
    ref_node_xyz(ref_node, 1, new_node) = xyz[1];
    ref_node_xyz(ref_node, 2, new_node) = xyz[2];
  }

  fclose(file);

  /* check .1.face when paranoid, but tetgen -z should not mess with them */
  problem = REF_FALSE;
  snprintf(filename, 896, "%s-tetgen.1.face", project);
  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");
  REIS(1, fscanf(file, "%d", &ntri), "face header ntri");
  REIS(1, fscanf(file, "%d", &mark), "face header mark");
  REIS(1, mark, "face have mark");

  ref_cell = ref_grid_tri(ref_grid);
  ref_malloc_init(position, ref_cell_max(ref_cell), REF_INT, REF_EMPTY);
  for (cell = 0; cell < ntri; cell++) {
    REIS(1, fscanf(file, "%d", &item), "tri item");
    RES(cell, item, "tri index");
    for (node = 0; node < 3; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tri");
    if (1 == mark) REIS(1, fscanf(file, "%d", &id), "tri mark id");
    if (REF_SUCCESS == ref_cell_with(ref_cell, nodes, &new_cell)) {
      position[new_cell] = cell;
    } else {
      problem = REF_TRUE;
      ref_node_location(ref_node, nodes[0]);
      ref_node_location(ref_node, nodes[1]);
      ref_node_location(ref_node, nodes[2]);
      REF_WHERE("tetgen face tri not found in ref_grid");
    }
  }
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (REF_EMPTY == position[cell]) {
      problem = REF_TRUE;
      ref_node_location(ref_node, nodes[0]);
      ref_node_location(ref_node, nodes[1]);
      ref_node_location(ref_node, nodes[2]);
      printf("face id %d\n", nodes[3]);
      REF_WHERE("ref_grid tri not found in tetgen face");
    }
  }
  ref_free(position);
  RAS(!problem, "problem detected in tetgen triangles");

  fclose(file);

  snprintf(filename, 896, "%s-tetgen.1.ele", project);
  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  REIS(1, fscanf(file, "%d", &ntet), "ele header ntet");
  REIS(1, fscanf(file, "%d", &node_per), "ele header node_per");
  REIS(4, node_per, "expected tets");
  REIS(1, fscanf(file, "%d", &mark), "ele header mark");
  REIS(0, mark, "ele have mark");

  ref_cell = ref_grid_tet(ref_grid);
  for (cell = 0; cell < ntet; cell++) {
    REIS(1, fscanf(file, "%d", &item), "tet item");
    RES(cell, item, "tet index");
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tet");
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new tet");
    RES(cell, new_cell, "tet index");
  }

  fclose(file);

  ref_grid_surf(ref_grid) = REF_FALSE;

  if (delete_temp_files) {
    snprintf(filename, 896, "%s-tetgen.1.edge", project);
    REIS(0, remove(filename), "rm .1.edge tetgen output file");
    snprintf(filename, 896, "%s-tetgen.1.face", project);
    REIS(0, remove(filename), "rm .1.face tetgen output file");
    snprintf(filename, 896, "%s-tetgen.1.node", project);
    REIS(0, remove(filename), "rm .1.node tetgen output file");
    snprintf(filename, 896, "%s-tetgen.1.ele", project);
    REIS(0, remove(filename), "rm .1.ele tetgen output file");
    snprintf(filename, 896, "%s-tetgen.poly", project);
    REIS(0, remove(filename), "rm .poly tetgen input file");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_import_ugrid_tets(REF_GRID ref_grid,
                                        const char *filename) {
  REF_CELL ref_cell;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_DBL xyz[3];
  REF_INT new_node, orig_nnode, node, tri;
  REF_GLOB global;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT qua;
  REF_INT face_id;
  REF_INT cell, new_cell;

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RES(1, fscanf(file, "%d", &nnode), "nnode");
  RES(1, fscanf(file, "%d", &ntri), "ntri");
  RES(1, fscanf(file, "%d", &nqua), "nqua");
  RES(1, fscanf(file, "%d", &ntet), "ntet");
  RES(1, fscanf(file, "%d", &npyr), "npyr");
  RES(1, fscanf(file, "%d", &npri), "npri");
  RES(1, fscanf(file, "%d", &nhex), "nhex");

  orig_nnode = ref_node_n(ref_node);

  for (node = 0; node < nnode; node++) {
    REIS(1, fscanf(file, "%lf", &(xyz[0])), "x");
    REIS(1, fscanf(file, "%lf", &(xyz[1])), "y");
    REIS(1, fscanf(file, "%lf", &(xyz[2])), "z");
    if (node >= orig_nnode) {
      RSS(ref_node_next_global(ref_node, &global), "next global");
      REIS(node, global, "global node index");
      RSS(ref_node_add(ref_node, global, &new_node), "new_node");
      REIS(node, new_node, "node index");
      ref_node_xyz(ref_node, 0, new_node) = xyz[0];
      ref_node_xyz(ref_node, 1, new_node) = xyz[1];
      ref_node_xyz(ref_node, 2, new_node) = xyz[2];
    }
  }

  for (tri = 0; tri < ntri; tri++) {
    for (node = 0; node < 3; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tri");
  }
  for (qua = 0; qua < nqua; qua++) {
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "qua");
  }

  for (tri = 0; tri < ntri; tri++) {
    RES(1, fscanf(file, "%d", &face_id), "tri id");
  }

  for (qua = 0; qua < nqua; qua++) {
    RES(1, fscanf(file, "%d", &face_id), "qua id");
  }

  ref_cell = ref_grid_tet(ref_grid);
  for (cell = 0; cell < ntet; cell++) {
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tet");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new tet");
    RES(cell, new_cell, "tet index");
  }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_aflr_volume(REF_GRID ref_grid, const char *project,
                                const char *options) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  char filename[1024];
  char command[1024];
  int system_status;
  REF_INT nnode_surface;

  printf("%d surface nodes %d triangles\n", ref_node_n(ref_node),
         ref_cell_n(ref_grid_tri(ref_grid)));

  snprintf(filename, 1024, "%s-aflr-surface.lb8.ugrid", project);
  RSS(ref_export_by_extension(ref_grid, filename), "ugrid");
  if (NULL == options) {
    sprintf(
        command,
        "aflr3 -igrid %s-aflr-surface.lb8.ugrid -ogrid %s-aflr-volume.ugrid "
        "-mrecrbf=0 -angqbf=179.9 -angqbfmin=0.1 "
        "< /dev/null > %s-aflr.txt",
        project, project, project);
  } else {
    sprintf(
        command,
        "aflr3 -igrid %s-aflr-surface.lb8.ugrid -ogrid %s-aflr-volume.ugrid "
        "%s "
        "< /dev/null > %s-aflr.txt",
        project, project, options, project);
  }
  printf("%s\n", command);
  fflush(stdout);
  system_status = system(command);
  REIS(0, system_status, "aflr failed");

  nnode_surface = ref_node_n(ref_node);

  snprintf(filename, 1024, "%s-aflr-volume.ugrid", project);
  RSS(ref_import_ugrid_tets(ref_grid, filename), "tets only");

  printf("%d interior nodes\n", ref_node_n(ref_node) - nnode_surface);

  ref_grid_surf(ref_grid) = REF_FALSE;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_infer_nedge_nface(REF_GRID ref_grid) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT min_id, max_id;
  RSS(ref_cell_id_range(ref_grid_tri(ref_grid), ref_grid_mpi(ref_grid), &min_id,
                        &max_id),
      "face range");
  REIS(1, min_id, "first face id not 1");
  ref_geom->nface = max_id;
  RSS(ref_cell_id_range(ref_grid_edg(ref_grid), ref_grid_mpi(ref_grid), &min_id,
                        &max_id),
      "edge range");
  REIS(1, min_id, "first edge id not 1");
  ref_geom->nedge = max_id;
  return REF_SUCCESS;
}

REF_STATUS ref_geom_usable(REF_GEOM ref_geom, REF_INT geom, REF_BOOL *usable) {
  REF_DBL kr0, r0[3], ks0, s0[3];
  REF_DBL curvature_is_ok = 1.0;
  REF_DBL uv0[2], xyz0[3], dxyz_dtuv[15];
  REF_DBL h, delta_radian = 0.1; /* 10 segment per radian of curvature */
  REF_DBL drsduv[4], duvdrs[4];
  REF_DBL duv[2], drs[2];
  REF_DBL uv[2];
  REF_DBL kr, r[3], ks, s[3];
  REF_DBL kmin = 0.1, kmax = 10.0;
  REF_BOOL verbose = REF_FALSE;
  *usable = REF_FALSE;
  if (geom < 0 || ref_geom_max(ref_geom) <= geom) return REF_INVALID;

  if (REF_GEOM_FACE == ref_geom_type(ref_geom, geom)) {
    uv0[0] = ref_geom_param(ref_geom, 0, geom); /* ignores periodic */
    uv0[1] = ref_geom_param(ref_geom, 1, geom);

    if (REF_SUCCESS !=
        ref_egads_face_curvature(ref_geom, geom, &kr0, r0, &ks0, s0)) {
      *usable = REF_FALSE;
      return REF_SUCCESS;
    }

    if (ABS(kr0) < curvature_is_ok && ABS(ks0) < curvature_is_ok) {
      *usable = REF_TRUE;
      return REF_SUCCESS;
    }

    RSS(ref_egads_eval(ref_geom, geom, xyz0, dxyz_dtuv), "eval edge");
    /* [x_u,y_u,z_u] [x_v,y_v,z_v] */
    /*  drsduv = [r s] * dxyz_dtuv */
    drsduv[0] =
        dxyz_dtuv[0] * r0[0] + dxyz_dtuv[1] * r0[1] + dxyz_dtuv[2] * r0[2];
    drsduv[1] =
        dxyz_dtuv[0] * s0[0] + dxyz_dtuv[1] * s0[1] + dxyz_dtuv[2] * s0[2];
    drsduv[2] =
        dxyz_dtuv[3] * r0[0] + dxyz_dtuv[4] * r0[1] + dxyz_dtuv[5] * r0[2];
    drsduv[3] =
        dxyz_dtuv[3] * s0[0] + dxyz_dtuv[4] * s0[1] + dxyz_dtuv[5] * s0[2];
    if (REF_SUCCESS != ref_matrix_inv_gen(2, drsduv, duvdrs)) {
      *usable = REF_FALSE;
      return REF_SUCCESS;
    }

    if (ABS(kr0) > curvature_is_ok) {
      /* find points +/- h from xyz0 along r */
      h = delta_radian / ABS(kr0);
      drs[0] = h;
      drs[1] = 0.0;
      duv[0] = duvdrs[0] * drs[0] + duvdrs[2] * drs[1];
      duv[1] = duvdrs[1] * drs[0] + duvdrs[3] * drs[1];
      uv[0] = uv0[0] + duv[0];
      uv[1] = uv0[1] + duv[1];
      if (REF_SUCCESS !=
          ref_egads_face_curvature_at(ref_geom, ref_geom_id(ref_geom, geom),
                                      ref_geom_degen(ref_geom, geom), uv, &kr,
                                      r, &ks, s))
        return REF_SUCCESS;
      if (ABS(kr) < kmin * ABS(kr0) || kmax * ABS(kr0) < ABS(kr)) {
        if (verbose) printf("kr0 %f kr+ %f\n", kr0, kr);
        *usable = REF_FALSE;
        return REF_SUCCESS;
      }
      h = delta_radian / ABS(kr0);
      drs[0] = -h;
      drs[1] = 0.0;
      duv[0] = duvdrs[0] * drs[0] + duvdrs[2] * drs[1];
      duv[1] = duvdrs[1] * drs[0] + duvdrs[3] * drs[1];
      uv[0] = uv0[0] + duv[0];
      uv[1] = uv0[1] + duv[1];
      if (REF_SUCCESS !=
          ref_egads_face_curvature_at(ref_geom, ref_geom_id(ref_geom, geom),
                                      ref_geom_degen(ref_geom, geom), uv, &kr,
                                      r, &ks, s))
        return REF_SUCCESS;
      if (ABS(kr) < kmin * ABS(kr0) || kmax * ABS(kr0) < ABS(kr)) {
        if (verbose) printf("kr0 %f kr- %f\n", kr0, kr);
        *usable = REF_FALSE;
        return REF_SUCCESS;
      }
    }
    if (ABS(ks0) > curvature_is_ok) {
      /* find points +/- h from xyz0 along s */
      h = delta_radian / ks0;
      drs[0] = 0.0;
      drs[1] = h;
      duv[0] = duvdrs[0] * drs[0] + duvdrs[2] * drs[1];
      duv[1] = duvdrs[1] * drs[0] + duvdrs[3] * drs[1];
      uv[0] = uv0[0] + duv[0];
      uv[1] = uv0[1] + duv[1];
      if (REF_SUCCESS !=
          ref_egads_face_curvature_at(ref_geom, ref_geom_id(ref_geom, geom),
                                      ref_geom_degen(ref_geom, geom), uv, &kr,
                                      r, &ks, s))
        return REF_SUCCESS;
      if (ABS(ks) < kmin * ABS(ks0) || kmax * ABS(ks0) < ABS(ks)) {
        if (verbose)
          printf("ks0 %f ks+ %f duv %f %f\n", ks0, ks, duv[0], duv[1]);
        *usable = REF_FALSE;
        return REF_SUCCESS;
      }
      if (verbose)
        printf("keep ks0 %f ks+ %f duv %f %f\n", ks0, ks, duv[0], duv[1]);
      h = delta_radian / ks0;
      drs[0] = 0.0;
      drs[1] = -h;
      duv[0] = duvdrs[0] * drs[0] + duvdrs[2] * drs[1];
      duv[1] = duvdrs[1] * drs[0] + duvdrs[3] * drs[1];
      uv[0] = uv0[0] + duv[0];
      uv[1] = uv0[1] + duv[1];
      if (REF_SUCCESS !=
          ref_egads_face_curvature_at(ref_geom, ref_geom_id(ref_geom, geom),
                                      ref_geom_degen(ref_geom, geom), uv, &kr,
                                      r, &ks, s))
        return REF_SUCCESS;
      if (ABS(ks) < kmin * ABS(ks0) || kmax * ABS(ks0) < ABS(ks)) {
        if (verbose)
          printf("ks0 %f ks- %f duv %f %f\n", ks0, ks, duv[0], duv[1]);
        *usable = REF_FALSE;
        return REF_SUCCESS;
      }
      if (verbose)
        printf("keep ks0 %f ks- %f duv %f %f\n", ks0, ks, duv[0], duv[1]);
    }
  }

  *usable = REF_TRUE;
  return REF_SUCCESS;
}

REF_STATUS ref_geom_reliability(REF_GEOM ref_geom, REF_INT geom,
                                REF_DBL *slop) {
  REF_DBL tol, gap;
  *slop = 0.0;
  RSS(ref_egads_tolerance(ref_geom, ref_geom_type(ref_geom, geom),
                          ref_geom_id(ref_geom, geom), &tol),
      "tol");
  *slop = MAX(*slop, ref_geom_tolerance_protection(ref_geom) * tol);
  /*
    each_ref_geom_having_node(ref_geom, node, item, geom) {
      RSS(ref_geom_tolerance(ref_geom, ref_geom_type(ref_geom, geom),
                             ref_geom_id(ref_geom, geom), &tol),
          "tol");
      *slop = MAX(*slop, ref_geom_tolerance_protection(ref_geom) * tol);
    }
  */
  RSS(ref_egads_gap(ref_geom, ref_geom_node(ref_geom, geom), &gap), "gap");
  *slop = MAX(*slop, ref_geom_gap_protection(ref_geom) * gap);
  return REF_SUCCESS;
}

static REF_STATUS ref_geom_node_min_angle(REF_GRID ref_grid, REF_INT node,
                                          REF_DBL *angle) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  REF_INT item1, cell1, node1;
  REF_INT item2, cell2, node2;
  REF_INT i;
  REF_DBL dot, dx1[3], dx2[3];
  *angle = 180.0;
  each_ref_cell_having_node(ref_cell, node, item1, cell1) {
    each_ref_cell_having_node(ref_cell, node, item2, cell2) {
      if (cell1 == cell2) continue; /* skip same edge */
      node1 = ref_cell_c2n(ref_cell, 0, cell1) +
              ref_cell_c2n(ref_cell, 1, cell1) - node;
      node2 = ref_cell_c2n(ref_cell, 0, cell2) +
              ref_cell_c2n(ref_cell, 1, cell2) - node;
      for (i = 0; i < 3; i++) {
        dx1[i] =
            ref_node_xyz(ref_node, i, node1) - ref_node_xyz(ref_node, i, node);
        dx2[i] =
            ref_node_xyz(ref_node, i, node2) - ref_node_xyz(ref_node, i, node);
      }
      RSS(ref_math_normalize(dx1), "dx1");
      RSS(ref_math_normalize(dx2), "dx2");
      dot = ref_math_dot(dx1, dx2);
      *angle = MIN(*angle, ref_math_in_degrees(acos(dot)));
    }
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_geom_node_short_edge(REF_GRID ref_grid, REF_INT node,
                                           REF_DBL *short_edge,
                                           REF_DBL *short_diag,
                                           REF_INT *short_id) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT item, geom;
  REF_DBL diag, min_diag, max_diag;
  *short_edge = 1.0;
  *short_diag = REF_DBL_MAX;
  *short_id = REF_EMPTY;
  min_diag = REF_DBL_MAX;
  max_diag = REF_DBL_MIN;
  each_ref_geom_having_node(ref_geom, node, item, geom) {
    if (REF_GEOM_EDGE != ref_geom_type(ref_geom, geom)) continue;
    RSS(ref_egads_diagonal(ref_geom, geom, &diag), "edge diag");
    if (diag < min_diag) {
      min_diag = diag;
      *short_diag = diag;
      *short_id = ref_geom_id(ref_geom, geom);
    }
    max_diag = MAX(diag, max_diag);
  }
  if (ref_math_divisible(min_diag, max_diag)) {
    *short_edge = min_diag / max_diag;
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_geom_face_curve_tol(REF_GRID ref_grid, REF_INT faceid,
                                          REF_DBL *curve, REF_DBL *location) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);

  REF_INT edge_geom, node;
  REF_INT item, face_geom;
  REF_DBL hmax, delta_radian, rlimit;
  REF_DBL kr, r[3], ks, s[3];
  REF_DBL hr, hs, slop;
  REF_DBL curvature_ratio = 1.0 / 20.0;

  *curve = 2.0;
  location[0] = 0.0;
  location[1] = 0.0;
  location[2] = 0.0;

  RSS(ref_egads_diagonal(ref_geom, REF_EMPTY, &hmax), "bbox diag");

  each_ref_geom_face(ref_geom, face_geom) {
    if (faceid != ref_geom_id(ref_geom, face_geom)) continue;
    node = ref_geom_node(ref_geom, face_geom);
    each_ref_geom_having_node(ref_geom, node, item, edge_geom) {
      if (REF_GEOM_EDGE != ref_geom_type(ref_geom, edge_geom)) continue;
      RSS(ref_geom_radian_request(ref_geom, edge_geom, &delta_radian), "drad");
      rlimit = hmax / delta_radian; /* h = r*drad, r = h/drad */
      RSS(ref_egads_face_curvature(ref_geom, face_geom, &kr, r, &ks, s),
          "curve");
      /* ignore sign, k is 1 / radius */
      kr = ABS(kr);
      ks = ABS(ks);
      /* limit the aspect ratio of the metric by reducing the largest radius */
      kr = MAX(kr, curvature_ratio * ks);
      ks = MAX(ks, curvature_ratio * kr);
      hr = hmax;
      if (1.0 / rlimit < kr) hr = delta_radian / kr;
      hs = hmax;
      if (1.0 / rlimit < ks) hs = delta_radian / ks;

      RSS(ref_geom_reliability(ref_geom, face_geom, &slop), "edge tol");
      if (hr < slop) {
        if (*curve > hr / slop) {
          *curve = hr / slop;
          location[0] = ref_node_xyz(ref_grid_node(ref_grid), 0, node);
          location[1] = ref_node_xyz(ref_grid_node(ref_grid), 1, node);
          location[2] = ref_node_xyz(ref_grid_node(ref_grid), 2, node);
        }
      }
      if (hs < slop) {
        if (*curve > hs / slop) {
          *curve = hs / slop;
          location[0] = ref_node_xyz(ref_grid_node(ref_grid), 0, node);
          location[1] = ref_node_xyz(ref_grid_node(ref_grid), 1, node);
          location[2] = ref_node_xyz(ref_grid_node(ref_grid), 2, node);
        }
      }
    }
  }
  return REF_SUCCESS;
}

REF_STATUS ref_geom_feedback(REF_GRID ref_grid, const char *filename) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT geom, node, edgeid, faceid;
  REF_DBL angle_tol = 10.0;
  REF_DBL angle;
  REF_DBL short_edge_tol = 1e-3;
  REF_DBL short_edge, diag;
  REF_DBL curve;
  REF_DBL location[3];
  REF_INT nsliver, nshort;

  if (ref_geom_effective(ref_geom)) {
    /* EFFECTIVE */
    if (ref_mpi_once(ref_grid_mpi(ref_grid)))
      printf("ref_geom_feedback does not support Effective EGADS\n");
    return REF_SUCCESS;
  }

  REIS(REF_MIGRATE_SINGLE, ref_grid_partitioner(ref_grid),
       "parallel implementation is incomplete");

  printf("Triaging geometry for issues impacting effiency and robustness\n");

  printf(
      "scaning for discrete face corners with slivers or cusps of %.1f deg or "
      "less\n",
      angle_tol);
  printf(
      "slivers in flat regions constrain the mesh and are efficency "
      "concerns\n");
  printf(
      "cusps (very small angles) in curved regions may be fatal with loose "
      "tolerences\n");
  printf(
      "  because the edge/face topology becomes ambiguous with adaptive "
      "refinement\n");
  nsliver = 0;
  each_ref_geom_node(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    RSS(ref_geom_node_min_angle(ref_grid, node, &angle), "node angle");
    if (angle <= angle_tol) {
      printf("%f %f %f # sliver deg=%f at geom node %d\n",
             ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
             ref_node_xyz(ref_node, 2, node), angle,
             ref_geom_id(ref_geom, geom));
      nsliver++;
    }
  }
  printf("%d geometry nodes with slivers or cusps of %.1f deg or less\n",
         nsliver, angle_tol);

  nshort = 0;
  if (ref_geom_model_loaded(ref_geom)) {
    printf(
        "scaning for geom nodes with shortest/longest edge ratios of %.1f or "
        "less\n",
        short_edge_tol);
    printf("short edges are efficency concerns that rarely create failures.\n");
    each_ref_geom_node(ref_geom, geom) {
      node = ref_geom_node(ref_geom, geom);
      RSS(ref_geom_node_short_edge(ref_grid, node, &short_edge, &diag, &edgeid),
          "short edge");
      if (short_edge <= short_edge_tol) {
        printf("%f %f %f # short edge diagonal %e ratio %e edge id %d\n",
               ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
               ref_node_xyz(ref_node, 2, node), diag, short_edge, edgeid);
        nshort++;
      }
    }
  }

  if (ref_geom_model_loaded(ref_geom)) {
    for (faceid = 1; faceid <= ref_geom->nedge; faceid++) {
      RSS(ref_geom_face_curve_tol(ref_grid, faceid, &curve, location),
          "curved face");
      if (curve < 1.0) {
        printf("%f %f %f# face id %d curve/tol %e\n", location[0], location[1],
               location[2], faceid, curve);
      }
    }
  }

  if (nsliver > 0 || nshort > 0) {
    FILE *file;
    REF_INT n;
    file = fopen(filename, "w");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");
    fprintf(file, "title=\"tecplot refine geometry triage\"\n");
    fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

    if (nsliver > 0) {
      printf("exporting slivers locations to %s\n", filename);
      fprintf(file, "zone t=\"sliver\", i=%d, datapacking=%s\n", nsliver,
              "point");
      n = 0;
      each_ref_geom_node(ref_geom, geom) {
        node = ref_geom_node(ref_geom, geom);
        RSS(ref_geom_node_min_angle(ref_grid, node, &angle), "node angle");
        if (angle <= angle_tol) {
          fprintf(file, "# sliver deg=%f at geom node %d\n %f %f %f\n", angle,
                  ref_geom_id(ref_geom, geom), ref_node_xyz(ref_node, 0, node),
                  ref_node_xyz(ref_node, 1, node),
                  ref_node_xyz(ref_node, 2, node));
          n++;
        }
      }
      REIS(nsliver, n, "tecplot sliver different recount");
    }

    if (nshort > 0) {
      printf("exporting short edge locations to %s\n", filename);
      fprintf(file, "zone t=\"short edge\", i=%d, datapacking=%s\n", nshort,
              "point");
      n = 0;
      each_ref_geom_node(ref_geom, geom) {
        node = ref_geom_node(ref_geom, geom);
        RSS(ref_geom_node_short_edge(ref_grid, node, &short_edge, &diag,
                                     &edgeid),
            "short edge");
        if (short_edge <= short_edge_tol) {
          fprintf(
              file, "# short edge diagonal %e ratio %e edge id %d\n%f %f %f\n",
              diag, short_edge, edgeid, ref_node_xyz(ref_node, 0, node),
              ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node));
          n++;
        }
      }
      REIS(nshort, n, "tecplot short edge different recount");
    }

    fclose(file);
  }
  return REF_SUCCESS;
}

REF_STATUS ref_geom_has_jump(REF_GEOM ref_geom, REF_INT node,
                             REF_BOOL *has_jump) {
  REF_INT item, geom;
  *has_jump = REF_FALSE;
  each_ref_geom_having_node(ref_geom, node, item, geom) {
    if (0 != ref_geom_jump(ref_geom, geom)) {
      *has_jump = REF_TRUE;
      return REF_SUCCESS;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_edge_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, nedg, sens;
  REF_INT jump_geom = REF_EMPTY;
  REF_DBL *t, tvalue;
  REF_DBL radius, normal[3], xyz[3], gap;

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_geom_edge(ref_geom, geom) {
    if (id == ref_geom_id(ref_geom, geom)) {
      RSS(ref_dict_store(ref_dict, ref_geom_node(ref_geom, geom), geom),
          "mark nodes");
      if (0 != ref_geom_jump(ref_geom, geom)) {
        REIS(REF_EMPTY, jump_geom, "should be only one jump per edge");
        jump_geom = geom;
      }
    }
  }
  nnode = ref_dict_n(ref_dict);
  if (REF_EMPTY != jump_geom) nnode++;

  nedg = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[2]) {
      nedg++;
    }
  }

  /* skip degenerate */
  if (0 == nnode || 0 == nedg) {
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(
      file,
      "zone t=\"edge%d\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      id, nnode, nedg, "point", "felineseg");

  ref_malloc(t, nnode, REF_DBL);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[2]) {
      RSB(ref_dict_location(ref_dict, nodes[0], &local), "localize", {
        printf("edg %d %d id %d no edge geom\n", nodes[0], nodes[1], nodes[2]);
        RSS(ref_node_location(ref_node, nodes[0]), "loc");
        RSS(ref_geom_tattle(ref_geom, nodes[0]), "tatt");
      });
      RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_EDGE, &tvalue,
                            &sens),
          "from");
      if (-1 == sens) local = nnode - 1;
      t[local] = tvalue;
      RSS(ref_dict_location(ref_dict, nodes[1], &local), "localize");
      RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_EDGE, &tvalue,
                            &sens),
          "from");
      if (-1 == sens) local = nnode - 1;
      t[local] = tvalue;
    }
  }

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    radius = 0;
    gap = 0;
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_edge_curvature(ref_geom, geom, &radius, normal), "curve");
      radius = ABS(radius);
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_EDGE, id, &(t[item]), xyz, NULL),
          "eval at");
      RSS(ref_egads_gap(ref_geom, node, &gap), "gap")
    }
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], gap, t[item], 0.0, radius, 0.0);
  }
  if (REF_EMPTY != jump_geom) {
    node = ref_geom_node(ref_geom, jump_geom);
    radius = 0;
    gap = 0;
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_edge_curvature(ref_geom, jump_geom, &radius, normal),
          "curve");
      radius = ABS(radius);
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_EDGE, id, &(t[nnode - 1]), xyz,
                            NULL),
          "eval at");
      RSS(ref_egads_gap(ref_geom, node, &gap), "gap")
    }
    node = ref_geom_node(ref_geom, jump_geom);
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], gap, t[nnode - 1], 0.0, radius, 0.0);
  }
  ref_free(t);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[2]) {
      RSS(ref_dict_location(ref_dict, nodes[0], &local), "localize");
      RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_EDGE, &tvalue,
                            &sens),
          "from");
      if (-1 == sens) local = nnode - 1;
      fprintf(file, " %d", local + 1);
      RSS(ref_dict_location(ref_dict, nodes[1], &local), "localize");
      RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_EDGE, &tvalue,
                            &sens),
          "from");
      if (-1 == sens) local = nnode - 1;
      fprintf(file, " %d", local + 1);
      fprintf(file, "\n");
    }
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_pcrv_tec_zone(REF_GRID ref_grid, REF_INT edgeid,
                                  REF_INT faceid, REF_INT sense, FILE *file);
REF_STATUS ref_geom_pcrv_tec_zone(REF_GRID ref_grid, REF_INT edgeid,
                                  REF_INT faceid, REF_INT sense, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, nedg, sens;
  REF_INT jump_geom = REF_EMPTY;
  REF_DBL *t, *uv, tvalue;
  REF_DBL radius, normal[3], xyz[3], gap;

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_geom_edge(ref_geom, geom) {
    if (edgeid == ref_geom_id(ref_geom, geom)) {
      RSS(ref_dict_store(ref_dict, ref_geom_node(ref_geom, geom), geom),
          "mark nodes");
      if (0 != ref_geom_jump(ref_geom, geom)) {
        REIS(REF_EMPTY, jump_geom, "should be only one jump per edge");
        jump_geom = geom;
      }
    }
  }
  nnode = ref_dict_n(ref_dict);
  if (REF_EMPTY != jump_geom) nnode++;

  nedg = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (edgeid == nodes[2]) {
      nedg++;
    }
  }

  /* skip degenerate */
  if (0 == nnode || 0 == nedg) {
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(file,
          "zone t=\"pcrv%df%d\", nodes=%d, elements=%d, datapacking=%s, "
          "zonetype=%s\n",
          edgeid, faceid, nnode, nedg, "point", "felineseg");

  ref_malloc(t, nnode, REF_DBL);
  ref_malloc(uv, 2 * nnode, REF_DBL);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (edgeid == nodes[2]) {
      RSB(ref_dict_location(ref_dict, nodes[0], &local), "localize", {
        printf("edg %d %d id %d no edge geom\n", nodes[0], nodes[1], nodes[2]);
        RSS(ref_node_location(ref_node, nodes[0]), "loc");
        RSS(ref_geom_tattle(ref_geom, nodes[0]), "tatt");
      });
      RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_EDGE, &tvalue,
                            &sens),
          "from");
      if (-1 == sens) local = nnode - 1;
      t[local] = tvalue;
      RSS(ref_egads_edge_face_uv(ref_geom, edgeid, faceid, sense, t[local],
                                 &(uv[2 * local])),
          "p-curve uv");
      RSS(ref_dict_location(ref_dict, nodes[1], &local), "localize");
      RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_EDGE, &tvalue,
                            &sens),
          "from");
      if (-1 == sens) local = nnode - 1;
      t[local] = tvalue;
      RSS(ref_egads_edge_face_uv(ref_geom, edgeid, faceid, sense, t[local],
                                 &(uv[2 * local])),
          "p-curve uv");
    }
  }

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    radius = 0;
    gap = 0;
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_edge_curvature(ref_geom, geom, &radius, normal), "curve");
      radius = ABS(radius);
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_EDGE, edgeid, &(t[item]), xyz,
                            NULL),
          "eval at");
      RSS(ref_egads_edge_face_uv(ref_geom, edgeid, faceid, sense, t[item],
                                 &(uv[2 * item])),
          "p-curve uv");
      RSS(ref_egads_gap(ref_geom, node, &gap), "gap")
    }
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], gap, uv[0 + 2 * item], uv[1 + 2 * item], radius,
            t[item]);
  }
  if (REF_EMPTY != jump_geom) {
    node = ref_geom_node(ref_geom, jump_geom);
    radius = 0;
    gap = 0;
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_edge_curvature(ref_geom, jump_geom, &radius, normal),
          "curve");
      radius = ABS(radius);
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_EDGE, edgeid, &(t[nnode - 1]),
                            xyz, NULL),
          "eval at");
      RSS(ref_egads_edge_face_uv(ref_geom, edgeid, faceid, sense, t[nnode - 1],
                                 &(uv[2 * (nnode - 1)])),
          "p-curve uv");
      RSS(ref_egads_gap(ref_geom, node, &gap), "gap")
    }
    node = ref_geom_node(ref_geom, jump_geom);
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], gap, uv[0 + 2 * (nnode - 1)],
            uv[1 + 2 * (nnode - 1)], radius, t[nnode - 1]);
  }
  ref_free(uv);
  ref_free(t);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (edgeid == nodes[2]) {
      RSS(ref_dict_location(ref_dict, nodes[0], &local), "localize");
      RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_EDGE, &tvalue,
                            &sens),
          "from");
      if (-1 == sens) local = nnode - 1;
      fprintf(file, " %d", local + 1);
      RSS(ref_dict_location(ref_dict, nodes[1], &local), "localize");
      RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_EDGE, &tvalue,
                            &sens),
          "from");
      if (-1 == sens) local = nnode - 1;
      fprintf(file, " %d", local + 1);
      fprintf(file, "\n");
    }
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_face_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict, ref_dict_jump, ref_dict_degen;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item2, item, local, node;
  REF_INT nnode, nnode_sens0, nnode_degen, ntri;
  REF_INT sens;
  REF_DBL *uv, param[2];
  REF_DBL kr, r[3], ks, s[3], xyz[3], kmin, kmax, gap;

  RSS(ref_dict_create(&ref_dict), "create dict");
  RSS(ref_dict_create(&ref_dict_jump), "create dict");
  RSS(ref_dict_create(&ref_dict_degen), "create dict");

  each_ref_geom_face(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (id == ref_geom_id(ref_geom, geom)) {
      if (0 == ref_geom_degen(ref_geom, geom)) {
        RSS(ref_dict_store(ref_dict, node, geom), "mark nodes");
        if (0 != ref_geom_jump(ref_geom, geom)) {
          RSS(ref_dict_store(ref_dict_jump, node, geom), "mark jump");
        }
      } else {
        each_ref_cell_having_node(ref_cell, node, item, cell) {
          RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
          if (id == nodes[3]) {
            RSS(ref_dict_store(ref_dict_degen, cell, node), "mark degen");
          }
        }
      }
    }
  }

  nnode_sens0 = ref_dict_n(ref_dict);
  nnode_degen = ref_dict_n(ref_dict) + ref_dict_n(ref_dict_jump);
  nnode = ref_dict_n(ref_dict) + ref_dict_n(ref_dict_jump) +
          ref_dict_n(ref_dict_degen);

  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      ntri++;
    }
  }

  /* skip degenerate */
  if (0 == nnode || 0 == ntri) {
    RSS(ref_dict_free(ref_dict_degen), "free degen");
    RSS(ref_dict_free(ref_dict_jump), "free jump");
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(
      file,
      "zone t=\"face%d\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      id, nnode, ntri, "point", "fetriangle");

  ref_malloc_init(uv, 2 * nnode, REF_DBL, -1.0);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      each_ref_cell_cell_node(ref_cell, node) {
        RSS(ref_geom_find(ref_geom, nodes[node], REF_GEOM_FACE, id, &geom),
            "find");
        RSS(ref_geom_cell_tuv(ref_geom, nodes[node], nodes, REF_GEOM_FACE,
                              param, &sens),
            "cell tuv");
        if (0 == ref_geom_degen(ref_geom, geom)) {
          if (0 == sens || 1 == sens) {
            RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
          } else {
            RSS(ref_dict_location(ref_dict_jump, nodes[node], &local),
                "localize");
            local += nnode_sens0;
          }
        } else {
          RSS(ref_dict_location(ref_dict_degen, cell, &local), "localize");
          local += nnode_degen;
        }
        uv[0 + 2 * local] = param[0];
        uv[1 + 2 * local] = param[1];
      }
    }
  }

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    kr = 0;
    ks = 0;
    gap = 0;
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RXS(ref_egads_face_curvature(ref_geom, geom, &kr, r, &ks, s), REF_FAILURE,
          "curve");
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, id, &(uv[2 * item]), xyz,
                            NULL),
          "eval at");
      RSS(ref_egads_gap(ref_geom, node, &gap), "gap")
    }
    if (ref_geom_meshlinked(ref_geom)) {
      RSS(ref_meshlink_face_curvature(ref_grid, geom, &kr, r, &ks, s), "curve");
    }
    kmax = MAX(ABS(kr), ABS(ks));
    kmin = MIN(ABS(kr), ABS(ks));
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %16e\n", xyz[0],
            xyz[1], xyz[2], gap, uv[0 + 2 * item], uv[1 + 2 * item], kmax,
            kmin);
  }
  each_ref_dict_key_value(ref_dict_jump, item, node, geom) {
    kr = 0;
    ks = 0;
    gap = 0;
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RXS(ref_egads_face_curvature(ref_geom, geom, &kr, r, &ks, s), REF_FAILURE,
          "curve");
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, id,
                            &(uv[2 * (nnode_sens0 + item)]), xyz, NULL),
          "eval at");
      RSS(ref_egads_gap(ref_geom, node, &gap), "gap")
    }
    kmax = MAX(ABS(kr), ABS(ks));
    kmin = MAX(ABS(kr), ABS(ks));
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], gap, uv[0 + 2 * (nnode_sens0 + item)],
            uv[1 + 2 * (nnode_sens0 + item)], kmax, kmin);
  }
  each_ref_dict_key_value(ref_dict_degen, item, cell, node) {
    kr = 0;
    ks = 0;
    gap = 0;
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      each_ref_geom_having_node(ref_geom, node, item2, geom) {
        if (ref_geom_type(ref_geom, geom) == REF_GEOM_FACE &&
            ref_geom_id(ref_geom, geom) == id) {
          RXS(ref_egads_face_curvature(ref_geom, geom, &kr, r, &ks, s),
              REF_FAILURE, "curve");
        }
      }
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, id,
                            &(uv[2 * (nnode_degen + item)]), xyz, NULL),
          "eval at");
      RSS(ref_egads_gap(ref_geom, node, &gap), "gap")
    }
    kmax = MAX(ABS(kr), ABS(ks));
    kmin = MAX(ABS(kr), ABS(ks));
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], gap, uv[0 + 2 * (nnode_degen + item)],
            uv[1 + 2 * (nnode_degen + item)], kmax, kmin);
  }
  ref_free(uv);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      each_ref_cell_cell_node(ref_cell, node) {
        RSS(ref_geom_find(ref_geom, nodes[node], REF_GEOM_FACE, id, &geom),
            "find");
        RSS(ref_geom_cell_tuv(ref_geom, nodes[node], nodes, REF_GEOM_FACE,
                              param, &sens),
            "cell tuv");
        if (0 == ref_geom_degen(ref_geom, geom)) {
          if (0 == sens || 1 == sens) {
            RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
          } else {
            RSS(ref_dict_location(ref_dict_jump, nodes[node], &local),
                "localize");
            local += nnode_sens0;
          }
        } else {
          RSS(ref_dict_location(ref_dict_degen, cell, &local), "localize");
          local += nnode_degen;
        }
        fprintf(file, " %d", local + 1);
      }
      fprintf(file, "\n");
    }
  }

  RSS(ref_dict_free(ref_dict_degen), "free degen");
  RSS(ref_dict_free(ref_dict_jump), "free jump");
  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_norm_tec_zone(REF_GRID ref_grid, REF_INT id, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, ntri;
  REF_DBL r[3], s[3], n[3], uv[2];
  REF_DBL area_sign;

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_geom_face(ref_geom, geom) {
    if (id == ref_geom_id(ref_geom, geom)) {
      RSS(ref_dict_store(ref_dict, ref_geom_node(ref_geom, geom), geom),
          "mark nodes");
    }
  }
  nnode = ref_dict_n(ref_dict);

  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      ntri++;
    }
  }

  /* skip degenerate */
  if (0 == nnode || 0 == ntri) {
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(
      file,
      "zone t=\"norm%d\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      id, nnode, ntri, "point", "fetriangle");

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    RSS(ref_geom_find(ref_geom, node, REF_GEOM_FACE, id, &geom), "not found");
    uv[0] = ref_geom_param(ref_geom, 0, geom);
    uv[1] = ref_geom_param(ref_geom, 1, geom);
    RSS(ref_geom_face_rsn(ref_geom, id, uv, r, s, n), "rsn");
    RSS(ref_geom_uv_area_sign(ref_grid, id, &area_sign), "a sign");
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
            ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
            ref_node_xyz(ref_node, 2, node), -area_sign * n[0],
            -area_sign * n[1], -area_sign * n[2]);
  }

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      for (node = 0; node < 3; node++) {
        RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
        fprintf(file, " %d", local + 1);
      }
      fprintf(file, "\n");
    }
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

static REF_STATUS ref_geom_curve_tec_zone(REF_GRID ref_grid, REF_INT id,
                                          FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, ntri;
  REF_DBL kr, r[3], ks, s[3];
  kr = 0;
  r[0] = 0;
  r[1] = 0;
  r[2] = 0;
  ks = 0;
  s[0] = 0;
  s[1] = 0;
  s[2] = 0;

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_geom_face(ref_geom, geom) {
    if (id == ref_geom_id(ref_geom, geom)) {
      RSS(ref_dict_store(ref_dict, ref_geom_node(ref_geom, geom), geom),
          "mark nodes");
    }
  }
  nnode = ref_dict_n(ref_dict);

  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      ntri++;
    }
  }

  /* skip degenerate */
  if (0 == nnode || 0 == ntri) {
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(file,
          "zone t=\"curve%d\", nodes=%d, elements=%d, datapacking=%s, "
          "zonetype=%s\n",
          id, nnode, ntri, "point", "fetriangle");

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_face_curvature(ref_geom, geom, &kr, r, &ks, s), "curve");
    }
    if (ref_geom_meshlinked(ref_geom)) {
      RSS(ref_meshlink_face_curvature(ref_grid, geom, &kr, r, &ks, s), "curve");
    }
    fprintf(file,
            " %.16e %.16e %.16e %.16e %.16e %.16e %.16e "
            "%.16e %.16e %.16e %.16e\n",
            ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
            ref_node_xyz(ref_node, 2, node), kr, r[0], r[1], r[2], ks, s[0],
            s[1], s[2]);
  }

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      for (node = 0; node < 3; node++) {
        RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
        fprintf(file, " %d", local + 1);
      }
      fprintf(file, "\n");
    }
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tec(REF_GRID ref_grid, const char *filename) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  FILE *file;
  REF_INT geom, id, min_id, max_id;
  REF_INT *edge_faces;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"refine cad coupling in tecplot format\"\n");
  fprintf(
      file,
      "variables = \"x\" \"y\" \"z\" \"gap\" \"p0\" \"p1\" \"k0\" \"k1\"\n");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_edge(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_geom_edge_tec_zone(ref_grid, id, file), "tec edge");

  if (ref_geom_model_loaded(ref_geom)) {
    RSS(ref_egads_edge_faces(ref_geom, &edge_faces), "edge faces");
    for (id = min_id; id <= max_id; id++) {
      if (edge_faces[0 + 2 * (id - 1)] ==
          edge_faces[1 + 2 * (id - 1)]) { /* edge used twice by face */
        RSS(ref_geom_pcrv_tec_zone(ref_grid, id, edge_faces[0 + 2 * (id - 1)],
                                   1, file),
            "tec pcrv");
        RSS(ref_geom_pcrv_tec_zone(ref_grid, id, edge_faces[1 + 2 * (id - 1)],
                                   -1, file),
            "tec pcrv");
      } else { /* edge used by two faces */
        if (REF_EMPTY != edge_faces[0 + 2 * (id - 1)])
          RSS(ref_geom_pcrv_tec_zone(ref_grid, id, edge_faces[0 + 2 * (id - 1)],
                                     0, file),
              "tec pcrv");
        if (REF_EMPTY != edge_faces[1 + 2 * (id - 1)])
          RSS(ref_geom_pcrv_tec_zone(ref_grid, id, edge_faces[1 + 2 * (id - 1)],
                                     0, file),
              "tec pcrv");
      }
    }
    ref_free(edge_faces);
  }

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_geom_face_tec_zone(ref_grid, id, file), "tec face");

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_curve_tec(REF_GRID ref_grid, const char *filename) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  FILE *file;
  REF_INT geom, id, min_id, max_id;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"refine cad coupling in tecplot format\"\n");
  fprintf(file,
          "variables = \"x\" \"y\" \"z\" \"kr\" \"rx\" \"ry\" \"rz\" "
          "\"ks\" \"sx\" \"sy\" \"sz\"\n");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_geom_curve_tec_zone(ref_grid, id, file), "tec face");

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_tec_para_shard(REF_GRID ref_grid,
                                   const char *root_filename) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  char filename[1024];
  if (ref_mpi_para(ref_mpi)) {
    sprintf(filename, "%s_%04d_%04d.tec", root_filename, ref_mpi_n(ref_mpi),
            ref_mpi_rank(ref_mpi));
  } else {
    sprintf(filename, "%s.tec", root_filename);
  }
  RSS(ref_geom_tec(ref_grid, filename), "tec");
  return REF_SUCCESS;
}

REF_STATUS ref_geom_ghost(REF_GEOM ref_geom, REF_NODE ref_node) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT *a_nnode, *b_nnode;
  REF_INT a_nnode_total, b_nnode_total;
  REF_GLOB *a_global, *b_global;
  REF_INT *a_part, *b_part;
  REF_INT *a_ngeom, *b_ngeom;
  REF_INT a_ngeom_total, b_ngeom_total;
  REF_GLOB *a_descr, *b_descr;
  REF_GLOB global;
  REF_DBL *a_param, *b_param;
  REF_INT part, node, degree;
  REF_INT *a_next, *b_next;
  REF_INT local, item, geom, i;
  REF_INT descr[REF_GEOM_DESCR_SIZE];

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_next, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_next, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(a_nnode, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_nnode, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(a_ngeom, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_ngeom, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      a_nnode[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_nnode, b_nnode, REF_INT_TYPE),
      "alltoall nnodes");

  a_nnode_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_nnode_total += a_nnode[part];
  ref_malloc(a_global, a_nnode_total, REF_GLOB);
  ref_malloc(a_part, a_nnode_total, REF_INT);

  b_nnode_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_nnode_total += b_nnode[part];
  ref_malloc(b_global, b_nnode_total, REF_GLOB);
  ref_malloc(b_part, b_nnode_total, REF_INT);

  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    a_next[part] = a_next[part - 1] + a_nnode[part - 1];
  }

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      part = ref_node_part(ref_node, node);
      a_global[a_next[part]] = ref_node_global(ref_node, node);
      a_part[a_next[part]] = ref_mpi_rank(ref_mpi);
      a_next[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_global, a_nnode, b_global, b_nnode, 1,
                        REF_GLOB_TYPE),
      "alltoallv global");
  RSS(ref_mpi_alltoallv(ref_mpi, a_part, a_nnode, b_part, b_nnode, 1,
                        REF_INT_TYPE),
      "alltoallv global");

  for (node = 0; node < b_nnode_total; node++) {
    RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
    part = b_part[node];
    RSS(ref_adj_degree(ref_geom_adj(ref_geom), local, &degree), "deg");
    /* printf("%d: node %d global %d local %d part %d degree %d\n",
       ref_mpi_rank(ref_mpi), node,b_global[node], local, part, degree); */
    b_ngeom[part] += degree;
  }

  RSS(ref_mpi_alltoall(ref_mpi, b_ngeom, a_ngeom, REF_INT_TYPE),
      "alltoall ngeoms");

  a_ngeom_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_ngeom_total += a_ngeom[part];
  ref_malloc(a_descr, REF_GEOM_DESCR_SIZE * a_ngeom_total, REF_GLOB);
  ref_malloc(a_param, 2 * a_ngeom_total, REF_DBL);

  b_ngeom_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_ngeom_total += b_ngeom[part];
  ref_malloc(b_descr, REF_GEOM_DESCR_SIZE * b_ngeom_total, REF_GLOB);
  ref_malloc(b_param, 2 * b_ngeom_total, REF_DBL);

  b_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    b_next[part] = b_next[part - 1] + b_ngeom[part - 1];
  }

  for (node = 0; node < b_nnode_total; node++) {
    RSS(ref_node_local(ref_node, b_global[node], &local), "g2l");
    part = b_part[node];
    each_ref_geom_having_node(ref_geom, local, item, geom) {
      each_ref_descr(ref_geom, i) {
        b_descr[i + REF_GEOM_DESCR_SIZE * b_next[part]] =
            (REF_GLOB)ref_geom_descr(ref_geom, i, geom);
      }
      b_descr[REF_GEOM_DESCR_NODE + REF_GEOM_DESCR_SIZE * b_next[part]] =
          ref_node_global(ref_node, ref_geom_node(ref_geom, geom));
      b_param[0 + 2 * b_next[part]] = ref_geom_param(ref_geom, 0, geom);
      b_param[1 + 2 * b_next[part]] = ref_geom_param(ref_geom, 1, geom);
      b_next[part]++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, b_descr, b_ngeom, a_descr, a_ngeom,
                        REF_GEOM_DESCR_SIZE, REF_GLOB_TYPE),
      "alltoallv descr");
  RSS(ref_mpi_alltoallv(ref_mpi, b_param, b_ngeom, a_param, a_ngeom, 2,
                        REF_DBL_TYPE),
      "alltoallv param");

  for (geom = 0; geom < a_ngeom_total; geom++) {
    each_ref_descr(ref_geom, i) {
      descr[i] = (REF_INT)a_descr[i + REF_GEOM_DESCR_SIZE * geom];
    }
    global = a_descr[REF_GEOM_DESCR_NODE + REF_GEOM_DESCR_SIZE * geom];
    RSS(ref_node_local(ref_node, global, &local), "g2l");
    descr[REF_GEOM_DESCR_NODE] = local;
    RSS(ref_geom_add_with_descr(ref_geom, descr, &(a_param[2 * geom])),
        "add ghost");
  }

  free(b_param);
  free(b_descr);
  free(a_param);
  free(a_descr);
  free(b_part);
  free(b_global);
  free(a_part);
  free(a_global);
  free(b_ngeom);
  free(a_ngeom);
  free(b_nnode);
  free(a_nnode);
  free(b_next);
  free(a_next);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_report_tri_area_normdev(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER], id;
  REF_DBL min_normdev, min_area, max_area, min_uv_area, max_uv_area;
  REF_DBL normdev, area, uv_area, area_sign;
  REF_DBL local;

  min_normdev = 2.0;
  min_area = REF_DBL_MAX;
  max_area = REF_DBL_MIN;
  min_uv_area = REF_DBL_MAX;
  max_uv_area = REF_DBL_MIN;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "norm dev");
    min_normdev = MIN(min_normdev, normdev);
    RSS(ref_node_tri_area(ref_grid_node(ref_grid), nodes, &area), "vol");
    min_area = MIN(min_area, area);
    max_area = MAX(max_area, area);
    id = nodes[ref_cell_node_per(ref_cell)];
    RSS(ref_geom_uv_area_sign(ref_grid, id, &area_sign), "a sign");
    RSS(ref_geom_uv_area(ref_grid_geom(ref_grid), nodes, &uv_area), "uv area");
    uv_area *= area_sign;
    min_uv_area = MIN(min_uv_area, uv_area);
    max_uv_area = MAX(max_uv_area, uv_area);
  }
  local = min_normdev;
  RSS(ref_mpi_min(ref_mpi, &local, &min_normdev, REF_DBL_TYPE), "mpi min");
  local = min_area;
  RSS(ref_mpi_min(ref_mpi, &local, &min_area, REF_DBL_TYPE), "mpi min");
  local = min_uv_area;
  RSS(ref_mpi_min(ref_mpi, &local, &min_uv_area, REF_DBL_TYPE), "mpi min");
  local = max_area;
  RSS(ref_mpi_max(ref_mpi, &local, &max_area, REF_DBL_TYPE), "mpi max");
  local = max_uv_area;
  RSS(ref_mpi_max(ref_mpi, &local, &max_uv_area, REF_DBL_TYPE), "mpi max");

  if (ref_mpi_once(ref_mpi)) {
    printf("normdev %f area %.5e  %.5e uv area  %.5e  %.5e\n", min_normdev,
           min_area, max_area, min_uv_area, max_uv_area);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_bary3(REF_GEOM ref_geom, REF_INT *nodes, REF_DBL *uv,
                          REF_DBL *bary) {
  REF_DBL uv0[2], uv1[2], uv2[2];
  REF_INT sens;
  REF_DBL total;
  REF_INT geom0, geom1, geom2;
  RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_FACE, uv0, &sens),
      "uv0");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_FACE, uv1, &sens),
      "uv1");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[2], nodes, REF_GEOM_FACE, uv2, &sens),
      "uv2");

  RSS(ref_geom_find(ref_geom, nodes[0], REF_GEOM_FACE, nodes[3], &geom0), "g0");
  if (0 != ref_geom_degen(ref_geom, geom0)) {
    if (0 < ref_geom_degen(ref_geom, geom0)) {
      uv0[1] = uv[1];
    } else {
      uv0[0] = uv[0];
    }
  }

  RSS(ref_geom_find(ref_geom, nodes[1], REF_GEOM_FACE, nodes[3], &geom1), "g1");
  if (0 != ref_geom_degen(ref_geom, geom1)) {
    if (0 < ref_geom_degen(ref_geom, geom1)) {
      uv1[1] = uv[1];
    } else {
      uv1[0] = uv[0];
    }
  }

  RSS(ref_geom_find(ref_geom, nodes[2], REF_GEOM_FACE, nodes[3], &geom2), "g2");
  if (0 != ref_geom_degen(ref_geom, geom2)) {
    if (0 < ref_geom_degen(ref_geom, geom2)) {
      uv2[1] = uv[1];
    } else {
      uv2[0] = uv[0];
    }
  }

  bary[0] = -uv1[0] * uv[1] + uv2[0] * uv[1] + uv[0] * uv1[1] -
            uv2[0] * uv1[1] - uv[0] * uv2[1] + uv1[0] * uv2[1];
  bary[1] = -uv[0] * uv0[1] + uv2[0] * uv0[1] + uv0[0] * uv[1] -
            uv2[0] * uv[1] - uv0[0] * uv2[1] + uv[0] * uv2[1];
  bary[2] = -uv1[0] * uv0[1] + uv[0] * uv0[1] + uv0[0] * uv1[1] -
            uv[0] * uv1[1] - uv0[0] * uv[1] + uv1[0] * uv[1];

  total = bary[0] + bary[1] + bary[2];

  if (ref_math_divisible(bary[0], total) &&
      ref_math_divisible(bary[1], total) &&
      ref_math_divisible(bary[2], total)) {
    bary[0] /= total;
    bary[1] /= total;
    bary[2] /= total;
  } else {
    REF_INT i;
    printf("%s: %d: %s: div zero total %.18e norms %.18e %.18e %.18e\n",
           __FILE__, __LINE__, __func__, total, bary[0], bary[1], bary[2]);
    for (i = 0; i < 3; i++) bary[i] = 0.0;
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tri_uv_bounding_sphere3(REF_GEOM ref_geom, REF_INT *nodes,
                                            REF_DBL *center, REF_DBL *radius) {
  REF_DBL uv0[2], uv1[2], uv2[2];
  REF_INT sens;

  RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_FACE, uv0, &sens),
      "uv0");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_FACE, uv1, &sens),
      "uv1");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[2], nodes, REF_GEOM_FACE, uv2, &sens),
      "uv2");

  center[0] = (uv0[0] + uv1[0] + uv2[0]) / 3.0;
  center[1] = (uv0[1] + uv1[1] + uv2[1]) / 3.0;

  *radius = 0.0;
  *radius = MAX(*radius,
                sqrt(pow(uv0[0] - center[0], 2) + pow(uv0[1] - center[1], 2)));
  *radius = MAX(*radius,
                sqrt(pow(uv1[0] - center[0], 2) + pow(uv1[1] - center[1], 2)));
  *radius = MAX(*radius,
                sqrt(pow(uv2[0] - center[0], 2) + pow(uv2[1] - center[1], 2)));

  return REF_SUCCESS;
}

REF_STATUS ref_geom_bary2(REF_GEOM ref_geom, REF_INT *nodes, REF_DBL t,
                          REF_DBL *bary) {
  REF_DBL t0, t1;
  REF_INT sens;
  REF_DBL total;
  RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_EDGE, &t0, &sens),
      "uv0");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_EDGE, &t1, &sens),
      "uv1");

  bary[0] = t1 - t;
  bary[1] = t - t0;

  total = bary[0] + bary[1];

  if (ref_math_divisible(bary[0], total) &&
      ref_math_divisible(bary[1], total)) {
    bary[0] /= total;
    bary[1] /= total;
  } else {
    REF_INT i;
    printf("%s: %d: %s: div zero total %.18e norms %.18e %.18e\n", __FILE__,
           __LINE__, __func__, total, bary[0], bary[1]);
    for (i = 0; i < 2; i++) bary[i] = 0.0;
    return REF_DIV_ZERO;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_edg_t_bounding_sphere2(REF_GEOM ref_geom, REF_INT *nodes,
                                           REF_DBL *center, REF_DBL *radius) {
  REF_DBL t0, t1;
  REF_INT sens;

  RSS(ref_geom_cell_tuv(ref_geom, nodes[0], nodes, REF_GEOM_EDGE, &t0, &sens),
      "uv0");
  RSS(ref_geom_cell_tuv(ref_geom, nodes[1], nodes, REF_GEOM_EDGE, &t1, &sens),
      "uv1");

  *center = (t0 + t1) * 0.5;

  *radius = 0.0;
  *radius = MAX(*radius, ABS(t0 - (*center)));
  *radius = MAX(*radius, ABS(t1 - (*center)));

  return REF_SUCCESS;
}

REF_STATUS ref_geom_enrich2(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_INT edge, *edge_node, part, node;
  REF_GLOB global;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell, new_cell;

  RSS(ref_node_synchronize_globals(ref_node), "sync glob");

  RSS(ref_geom_constrain_all(ref_grid), "constrain");

  RSS(ref_edge_create(&ref_edge, ref_grid), "edge");
  ref_malloc_init(edge_node, ref_edge_n(ref_edge), REF_INT, REF_EMPTY);
  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (ref_mpi_rank(ref_mpi) == part) {
      RSS(ref_node_next_global(ref_node, &global), "next global");
      RSS(ref_node_add(ref_node, global, &node), "add node");
      edge_node[edge] = node;
      RSS(ref_node_interpolate_edge(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                                    ref_edge_e2n(ref_edge, 1, edge), 0.5, node),
          "new reals");
      RSS(ref_geom_add_constrain_midnode(
              ref_grid, ref_edge_e2n(ref_edge, 0, edge),
              ref_edge_e2n(ref_edge, 1, edge), 0.5, node),
          "new geom");
    }
  }

  RSS(ref_node_shift_new_globals(ref_node), "shift glob");

  each_ref_cell_valid_cell_with_nodes(ref_grid_edg(ref_grid), cell, nodes) {
    nodes[ref_cell_id_index(ref_grid_ed2(ref_grid))] =
        nodes[ref_cell_id_index(ref_grid_edg(ref_grid))];
    RSS(ref_edge_with(ref_edge, nodes[0], nodes[1], &edge), "find edge01");
    nodes[2] = edge_node[edge];
    RSS(ref_cell_add(ref_grid_ed2(ref_grid), nodes, &new_cell), "add");
  }

  each_ref_cell_valid_cell_with_nodes(ref_grid_tri(ref_grid), cell, nodes) {
    nodes[ref_cell_id_index(ref_grid_tr2(ref_grid))] =
        nodes[ref_cell_id_index(ref_grid_tri(ref_grid))];
    RSS(ref_edge_with(ref_edge, nodes[0], nodes[1], &edge), "find edge01");
    nodes[3] = edge_node[edge];
    RSS(ref_edge_with(ref_edge, nodes[1], nodes[2], &edge), "find edge12");
    nodes[4] = edge_node[edge];
    RSS(ref_edge_with(ref_edge, nodes[2], nodes[0], &edge), "find edge20");
    nodes[5] = edge_node[edge];
    RSS(ref_cell_add(ref_grid_tr2(ref_grid), nodes, &new_cell), "add");
  }
  each_ref_cell_valid_cell_with_nodes(ref_grid_tet(ref_grid), cell, nodes) {
    RSS(ref_edge_with(ref_edge, nodes[0], nodes[1], &edge), "find edge");
    nodes[4] = edge_node[edge];
    RSS(ref_edge_with(ref_edge, nodes[1], nodes[2], &edge), "find edge");
    nodes[5] = edge_node[edge];
    RSS(ref_edge_with(ref_edge, nodes[0], nodes[2], &edge), "find edge");
    nodes[6] = edge_node[edge];
    RSS(ref_edge_with(ref_edge, nodes[0], nodes[3], &edge), "find edge");
    nodes[7] = edge_node[edge];
    RSS(ref_edge_with(ref_edge, nodes[1], nodes[3], &edge), "find edge");
    nodes[8] = edge_node[edge];
    RSS(ref_edge_with(ref_edge, nodes[2], nodes[3], &edge), "find edge");
    nodes[9] = edge_node[edge];
    RSS(ref_cell_add(ref_grid_te2(ref_grid), nodes, &new_cell), "add");
  }

  ref_free(edge_node);
  RSS(ref_edge_free(ref_edge), "free edge");

  each_ref_cell_valid_cell(ref_grid_edg(ref_grid), cell) {
    ref_cell_remove(ref_grid_edg(ref_grid), cell);
  }
  each_ref_cell_valid_cell(ref_grid_tri(ref_grid), cell) {
    ref_cell_remove(ref_grid_tri(ref_grid), cell);
  }
  each_ref_cell_valid_cell(ref_grid_tet(ref_grid), cell) {
    ref_cell_remove(ref_grid_tet(ref_grid), cell);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_enrich3(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_INT edge, *edge_node, part, node;
  REF_GLOB global;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell, new_cell;
  REF_DBL t;

  RSS(ref_geom_constrain_all(ref_grid), "constrain");

  RSS(ref_edge_create(&ref_edge, ref_grid), "edge");
  ref_malloc_init(edge_node, 2 * ref_edge_n(ref_edge), REF_INT, REF_EMPTY);
  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (ref_mpi_rank(ref_mpi) == part) {
      t = 1.0 / 3.0;
      RSS(ref_node_next_global(ref_node, &global), "next global");
      RSS(ref_node_add(ref_node, global, &node), "add node");
      edge_node[0 + 2 * edge] = node;
      RSS(ref_node_interpolate_edge(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                                    ref_edge_e2n(ref_edge, 1, edge), t, node),
          "new node");
      RSS(ref_geom_add_constrain_midnode(
              ref_grid, ref_edge_e2n(ref_edge, 0, edge),
              ref_edge_e2n(ref_edge, 1, edge), t, node),
          "new geom");
      t = 2.0 / 3.0;
      RSS(ref_node_next_global(ref_node, &global), "next global");
      RSS(ref_node_add(ref_node, global, &node), "add node");
      edge_node[1 + 2 * edge] = node;
      RSS(ref_node_interpolate_edge(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                                    ref_edge_e2n(ref_edge, 1, edge), t, node),
          "new node");
      RSS(ref_geom_add_constrain_midnode(
              ref_grid, ref_edge_e2n(ref_edge, 0, edge),
              ref_edge_e2n(ref_edge, 1, edge), t, node),
          "new geom");
    }
  }

  each_ref_cell_valid_cell_with_nodes(ref_grid_edg(ref_grid), cell, nodes) {
    nodes[ref_cell_id_index(ref_grid_ed3(ref_grid))] =
        nodes[ref_cell_id_index(ref_grid_edg(ref_grid))];

    RSS(ref_edge_with(ref_edge, nodes[0], nodes[1], &edge), "find edge01");
    if (nodes[0] == ref_edge_e2n(ref_edge, 0, edge)) {
      nodes[2] = edge_node[0 + 2 * edge]; /* forward edge, tri side direction */
      nodes[3] = edge_node[1 + 2 * edge];
    } else {
      nodes[2] = edge_node[1 + 2 * edge]; /* reverse edge, tri side direction */
      nodes[3] = edge_node[0 + 2 * edge];
    }

    RSS(ref_cell_add(ref_grid_ed3(ref_grid), nodes, &new_cell), "add");
  }

  each_ref_cell_valid_cell_with_nodes(ref_grid_tri(ref_grid), cell, nodes) {
    RSS(ref_node_next_global(ref_node, &global), "next global");
    RSS(ref_node_add(ref_node, global, &node), "add node");
    RSS(ref_node_interpolate_face(ref_node, nodes[0], nodes[1], nodes[2], node),
        "new node");
    RSS(ref_geom_add_constrain_inside_midnode(ref_grid, nodes, node),
        "new node");
    nodes[9] = node;

    nodes[ref_cell_id_index(ref_grid_tr3(ref_grid))] =
        nodes[ref_cell_id_index(ref_grid_tri(ref_grid))];

    RSS(ref_edge_with(ref_edge, nodes[0], nodes[1], &edge), "find edge01");
    if (nodes[0] == ref_edge_e2n(ref_edge, 0, edge)) {
      nodes[3] = edge_node[0 + 2 * edge]; /* forward edge, tri side direction */
      nodes[4] = edge_node[1 + 2 * edge];
    } else {
      nodes[3] = edge_node[1 + 2 * edge]; /* reverse edge, tri side direction */
      nodes[4] = edge_node[0 + 2 * edge];
    }
    RSS(ref_edge_with(ref_edge, nodes[1], nodes[2], &edge), "find edge12");
    if (nodes[1] == ref_edge_e2n(ref_edge, 0, edge)) {
      nodes[5] = edge_node[0 + 2 * edge]; /* forward edge, tri side direction */
      nodes[6] = edge_node[1 + 2 * edge];
    } else {
      nodes[5] = edge_node[1 + 2 * edge]; /* reverse edge, tri side direction */
      nodes[6] = edge_node[0 + 2 * edge];
    }
    RSS(ref_edge_with(ref_edge, nodes[2], nodes[0], &edge), "find edge20");
    if (nodes[2] == ref_edge_e2n(ref_edge, 0, edge)) {
      nodes[7] = edge_node[0 + 2 * edge]; /* forward edge, tri side direction */
      nodes[8] = edge_node[1 + 2 * edge];
    } else {
      nodes[7] = edge_node[1 + 2 * edge]; /* reverse edge, tri side direction */
      nodes[8] = edge_node[0 + 2 * edge];
    }

    RSS(ref_cell_add(ref_grid_tr3(ref_grid), nodes, &new_cell), "add");
  }

  ref_free(edge_node);
  RSS(ref_edge_free(ref_edge), "free edge");

  each_ref_cell_valid_cell(ref_grid_edg(ref_grid), cell) {
    ref_cell_remove(ref_grid_edg(ref_grid), cell);
  }
  each_ref_cell_valid_cell(ref_grid_tri(ref_grid), cell) {
    ref_cell_remove(ref_grid_tri(ref_grid), cell);
  }

  return REF_SUCCESS;
}
