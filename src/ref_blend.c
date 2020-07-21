
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

#include "ref_blend.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_egads.h"
#include "ref_malloc.h"

#define ref_blend_grid(ref_blend) ((ref_blend)->grid)
#define ref_blend_geom(ref_blend) (ref_grid_geom(ref_blend_grid(ref_blend)))
#define ref_blend_displacement(ref_blend, ixyz, geom) \
  ((ref_blend)->displacement[(ixyz) + 3 * (geom)])
#define ref_blend_strong_bc(ref_blend, geom) ((ref_blend)->strong_bc[(geom)])
#define ref_blend_edge_search(ref_blend, iedge) \
  ((ref_blend)->edge_search[(iedge)])
#define ref_blend_face_search(ref_blend, iface) \
  ((ref_blend)->face_search[(iface)])

static REF_STATUS ref_blend_cache_search(REF_BLEND ref_blend) {
  REF_INT nedge, iedge, nface, iface;
  REF_GRID ref_grid = ref_blend_grid(ref_blend);
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL center[3], radius, scale = 2.0;

  ref_blend->edge_search = NULL;

  nedge = ref_geom->nedge;

  if (0 < nedge) {
    ref_cell = ref_grid_edg(ref_grid);

    ref_malloc_init(ref_blend->edge_search, nedge, REF_SEARCH, NULL);
    for (iedge = 0; iedge < nedge; iedge++) {
      RSS(ref_search_create(&(ref_blend_edge_search(ref_blend, iedge)),
                            ref_cell_n(ref_cell)),
          "create edge search");
    }

    /* cache each t edg */
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_geom_edg_t_bounding_sphere2(ref_geom, nodes, center, &radius),
          "bound with circle");
      center[1] = 0.0;
      center[2] = 0.0;
      iedge = nodes[2] - 1;
      RSS(ref_search_insert(ref_blend_edge_search(ref_blend, iedge), cell,
                            center, scale * radius),
          "ins");
    }
  }

  ref_blend->face_search = NULL;

  nface = ref_geom->nface;

  if (0 < nface) {
    ref_cell = ref_grid_tri(ref_grid);

    ref_malloc_init(ref_blend->face_search, nface, REF_SEARCH, NULL);
    for (iface = 0; iface < nface; iface++) {
      RSS(ref_search_create(&(ref_blend_face_search(ref_blend, iface)),
                            ref_cell_n(ref_cell)),
          "create face search");
    }

    /* cache each uv tri */
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_geom_tri_uv_bounding_sphere3(ref_geom, nodes, center, &radius),
          "bound with circle");
      center[2] = 0.0;
      iface = nodes[3] - 1;
      RSS(ref_search_insert(ref_blend_face_search(ref_blend, iface), cell,
                            center, scale * radius),
          "ins");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_blend_create(REF_BLEND *ref_blend_ptr, REF_GRID ref_grid) {
  REF_BLEND ref_blend;
  REF_INT n;

  ref_malloc(*ref_blend_ptr, 1, REF_BLEND_STRUCT);

  ref_blend = *ref_blend_ptr;

  RSS(ref_grid_deep_copy(&ref_blend_grid(ref_blend), ref_grid), "deep copy");
  n = ref_geom_max(ref_blend_geom(ref_blend));
  ref_malloc_init(ref_blend->displacement, 3 * n, REF_DBL, 0.0);
  ref_malloc_init(ref_blend->strong_bc, n, REF_BOOL, REF_FALSE);

  RSS(ref_blend_cache_search(ref_blend), "cache tri uv");

  return REF_SUCCESS;
}

REF_STATUS ref_blend_free(REF_BLEND ref_blend) {
  REF_INT i;
  if (NULL == (void *)ref_blend) return REF_NULL;

  if (NULL != ref_blend->edge_search) {
    for (i = 0; i < ref_blend_geom(ref_blend)->nedge; i++)
      ref_search_free(ref_blend_edge_search(ref_blend, i));
    ref_free(ref_blend->edge_search);
  }
  if (NULL != ref_blend->face_search) {
    for (i = 0; i < ref_blend_geom(ref_blend)->nface; i++)
      ref_search_free(ref_blend_face_search(ref_blend, i));
    ref_free(ref_blend->face_search);
  }
  ref_free(ref_blend->strong_bc);
  ref_free(ref_blend->displacement);
  ref_grid_free(ref_blend_grid(ref_blend));
  /* geom is a pointer to the orig */
  ref_free(ref_blend);

  return REF_SUCCESS;
}

static REF_STATUS ref_blend_solve(REF_BLEND ref_blend) {
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  REF_GRID ref_grid = ref_blend_grid(ref_blend);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT center_geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, cell_node, other_geom, center;
  REF_DBL disp[3], hits;

  /* psudo laplace */
  each_ref_geom_face(ref_geom, center_geom) {
    if (!ref_blend_strong_bc(ref_blend, center_geom)) {
      center = ref_geom_node(ref_geom, center_geom);
      hits = 0.0;
      disp[0] = 0.0;
      disp[1] = 0.0;
      disp[2] = 0.0;
      each_ref_cell_having_node(ref_cell, center, item, cell) {
        RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
        each_ref_cell_cell_node(ref_cell, cell_node) {
          if (nodes[cell_node] != center) {
            RSS(ref_geom_find(ref_geom, nodes[cell_node], REF_GEOM_FACE,
                              ref_geom_id(ref_geom, center_geom), &other_geom),
                "other geom");
            hits += 0.5;
            disp[0] += 0.5 * ref_blend_displacement(ref_blend, 0, other_geom);
            disp[1] += 0.5 * ref_blend_displacement(ref_blend, 1, other_geom);
            disp[2] += 0.5 * ref_blend_displacement(ref_blend, 2, other_geom);
          }
        }
        if (hits > 0.1) {
          ref_blend_displacement(ref_blend, 0, center_geom) = disp[0] / hits;
          ref_blend_displacement(ref_blend, 1, center_geom) = disp[1] / hits;
          ref_blend_displacement(ref_blend, 2, center_geom) = disp[2] / hits;
        }
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_blend_initialize(REF_BLEND ref_blend) {
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  REF_DBL edge_xyz[3], face_xyz[3];
  REF_INT i, edge_geom, face_geom, node, item;

  /* also displace edge and face to geom nodes */

  each_ref_geom_edge(ref_geom, edge_geom) {
    RSS(ref_egads_eval(ref_geom, edge_geom, edge_xyz, NULL), "eval edge");
    node = ref_geom_node(ref_geom, edge_geom);
    each_ref_geom_having_node(ref_geom, node, item, face_geom) {
      if (REF_GEOM_FACE == ref_geom_type(ref_geom, face_geom)) {
        RSS(ref_egads_eval(ref_geom, face_geom, face_xyz, NULL), "eval face");
        ref_blend_strong_bc(ref_blend, face_geom) = REF_TRUE;
        for (i = 0; i < 3; i++) {
          ref_blend_displacement(ref_blend, i, face_geom) =
              edge_xyz[i] - face_xyz[i];
        }
      }
    }
  }

  for (i = 0; i < 10; i++) RSS(ref_blend_solve(ref_blend), "solve");

  return REF_SUCCESS;
}

REF_STATUS ref_blend_attach(REF_GRID ref_grid) {
  REF_BLEND ref_blend;
  RSS(ref_blend_create(&ref_blend, ref_grid), "create");
  ref_geom_blend(ref_grid_geom(ref_grid)) = ref_blend;
  RSS(ref_blend_initialize(ref_blend), "init disp");
  return REF_SUCCESS;
}

REF_STATUS ref_blend_enclosing(REF_BLEND ref_blend, REF_INT type, REF_INT id,
                               REF_DBL *param, REF_INT *cell, REF_DBL *bary) {
  REF_GRID ref_grid = ref_blend_grid(ref_blend);
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  REF_LIST ref_list;
  REF_SEARCH ref_search;
  REF_CELL ref_cell;
  REF_DBL parampad[3];
  REF_DBL fuzz = 1.0e-12;
  REF_INT item, candidate, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL current_bary[3], best_bary, min_bary;
  REF_INT best_candidate;

  *cell = REF_EMPTY;
  bary[0] = 0.0;
  bary[1] = 0.0;
  bary[2] = 0.0;

  REIS(REF_GEOM_FACE, type, "only implemented for face uv");
  RSS(ref_list_create(&ref_list), "create list");
  ref_search = ref_blend_face_search(ref_blend, id - 1);
  ref_cell = ref_grid_tri(ref_grid);
  parampad[0] = param[0];
  parampad[1] = param[1];
  parampad[2] = 0.0;

  RSS(ref_search_touching(ref_search, ref_list, parampad, fuzz), "touching");
  RAS(0 < ref_list_n(ref_list), "list empty");
  best_candidate = REF_EMPTY;
  best_bary = -999.0;
  each_ref_list_item(ref_list, item) {
    candidate = ref_list_value(ref_list, item);
    RSS(ref_cell_nodes(ref_cell, candidate, nodes), "cell");
    RSS(ref_geom_bary3(ref_geom, nodes, param, current_bary), "bary");
    min_bary = MIN(MIN(current_bary[0], current_bary[1]), current_bary[2]);
    if (REF_EMPTY == best_candidate || min_bary > best_bary) {
      best_candidate = candidate;
      best_bary = min_bary;
    }
  }

  RUS(REF_EMPTY, best_candidate, "failed to find cell");

  *cell = best_candidate;
  RSS(ref_cell_nodes(ref_cell, best_candidate, nodes), "cell");
  RSS(ref_geom_bary3(ref_geom, nodes, param, bary), "bary");

  RSS(ref_list_free(ref_list), "free list");

  return REF_SUCCESS;
}

REF_STATUS ref_blend_eval_at(REF_BLEND ref_blend, REF_INT type, REF_INT id,
                             REF_DBL *params, REF_DBL *xyz,
                             REF_DBL *dxyz_dtuv) {
  REF_GRID ref_grid = ref_blend_grid(ref_blend);
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL bary[3];
  RSS(ref_egads_eval_at(ref_geom, type, id, params, xyz, dxyz_dtuv),
      "egads eval");
  if (REF_GEOM_FACE == type) {
    RSS(ref_blend_enclosing(ref_blend, type, id, params, &cell, bary),
        "enclose");
    RSS(ref_cell_nodes(ref_grid_tri(ref_grid), cell, nodes), "nodes");

    RSS(ref_geom_find(ref_geom, nodes[0], REF_GEOM_FACE, id, &geom), "find 0");
    xyz[0] += bary[0] * ref_blend_displacement(ref_blend, 0, geom);
    xyz[1] += bary[0] * ref_blend_displacement(ref_blend, 1, geom);
    xyz[2] += bary[0] * ref_blend_displacement(ref_blend, 2, geom);

    RSS(ref_geom_find(ref_geom, nodes[1], REF_GEOM_FACE, id, &geom), "find 1");
    xyz[0] += bary[1] * ref_blend_displacement(ref_blend, 0, geom);
    xyz[1] += bary[1] * ref_blend_displacement(ref_blend, 1, geom);
    xyz[2] += bary[1] * ref_blend_displacement(ref_blend, 2, geom);

    RSS(ref_geom_find(ref_geom, nodes[2], REF_GEOM_FACE, id, &geom), "find 2");
    xyz[0] += bary[2] * ref_blend_displacement(ref_blend, 0, geom);
    xyz[1] += bary[2] * ref_blend_displacement(ref_blend, 1, geom);
    xyz[2] += bary[2] * ref_blend_displacement(ref_blend, 2, geom);
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_blend_face_tec_zone(REF_BLEND ref_blend, REF_INT id,
                                          FILE *file) {
  REF_GRID ref_grid = ref_blend_grid(ref_blend);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  REF_DICT ref_dict, ref_dict_jump, ref_dict_degen;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, nnode_sens0, nnode_degen, ntri;
  REF_INT sens;
  REF_DBL *uv, param[2];
  REF_DBL xyz[3];

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
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);

    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], uv[0 + 2 * item], uv[1 + 2 * item],
            ref_blend_displacement(ref_blend, 0, geom),
            ref_blend_displacement(ref_blend, 1, geom),
            ref_blend_displacement(ref_blend, 2, geom));
  }
  each_ref_dict_key_value(ref_dict_jump, item, node, geom) {
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], uv[0 + 2 * (nnode_sens0 + item)],
            uv[1 + 2 * (nnode_sens0 + item)],
            ref_blend_displacement(ref_blend, 0, geom),
            ref_blend_displacement(ref_blend, 1, geom),
            ref_blend_displacement(ref_blend, 2, geom));
  }
  each_ref_dict_key_value(ref_dict_degen, item, cell, node) {
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], uv[0 + 2 * (nnode_degen + item)],
            uv[1 + 2 * (nnode_degen + item)],
            ref_blend_displacement(ref_blend, 0, geom),
            ref_blend_displacement(ref_blend, 1, geom),
            ref_blend_displacement(ref_blend, 2, geom));
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

REF_STATUS ref_blend_tec(REF_BLEND ref_blend, const char *filename) {
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  FILE *file;
  REF_INT geom, id, min_id, max_id;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"refine watertight cad displacement\"\n");
  fprintf(file,
          "variables = \"x\" \"y\" \"z\" \"u\" \"v\" \"dx\" \"dy\" \"dz\"\n");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_blend_face_tec_zone(ref_blend, id, file), "tec face");

  fclose(file);
  return REF_SUCCESS;
}
