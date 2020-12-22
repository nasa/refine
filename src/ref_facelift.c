
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

#include "ref_facelift.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_egads.h"
#include "ref_export.h"
#include "ref_import.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_node.h"
#include "ref_recon.h"

#define ref_facelift_direct(ref_facelift) (NULL == (ref_facelift)->displacement)
#define ref_facelift_geom(ref_facelift) \
  (ref_grid_geom(ref_facelift_grid(ref_facelift)))
#define ref_facelift_edg(ref_facelift) ((ref_facelift)->edg_cell)
#define ref_facelift_tri(ref_facelift) ((ref_facelift)->tri_cell)
#define ref_facelift_strong_bc(ref_facelift, geom) \
  ((ref_facelift)->strong_bc[(geom)])
#define ref_facelift_edge_search(ref_facelift, iedge) \
  ((ref_facelift)->edge_search[(iedge)])
#define ref_facelift_face_search(ref_facelift, iface) \
  ((ref_facelift)->face_search[(iface)])

static REF_STATUS ref_facelift_cache_search(REF_FACELIFT ref_facelift) {
  REF_INT nedge, iedge, nface, iface;
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL center[3], radius, scale = 2.0;

  ref_facelift->edge_search = NULL;

  nedge = ref_geom->nedge;

  if (0 < nedge) {
    ref_cell = ref_facelift_edg(ref_facelift);

    ref_malloc_init(ref_facelift->edge_search, nedge, REF_SEARCH, NULL);
    for (iedge = 0; iedge < nedge; iedge++) {
      RSS(ref_search_create(&(ref_facelift_edge_search(ref_facelift, iedge)),
                            ref_cell_n(ref_cell)),
          "create edge search");
    }

    /* cache each t edg */
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      nodes[2] = nodes[ref_cell_id_index(ref_cell)]; /* expects P1 */
      RSS(ref_geom_edg_t_bounding_sphere2(ref_geom, nodes, center, &radius),
          "bound with circle");
      center[1] = 0.0;
      center[2] = 0.0;
      iedge = nodes[2] - 1;
      RSS(ref_search_insert(ref_facelift_edge_search(ref_facelift, iedge), cell,
                            center, scale * radius),
          "ins");
    }
  }

  ref_facelift->face_search = NULL;

  nface = ref_geom->nface;

  if (0 < nface) {
    ref_cell = ref_facelift_tri(ref_facelift);

    ref_malloc_init(ref_facelift->face_search, nface, REF_SEARCH, NULL);
    for (iface = 0; iface < nface; iface++) {
      RSS(ref_search_create(&(ref_facelift_face_search(ref_facelift, iface)),
                            ref_cell_n(ref_cell)),
          "create face search");
    }

    /* cache each uv tri */
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      nodes[3] = nodes[ref_cell_id_index(ref_cell)]; /* expects P1 */
      RSS(ref_geom_tri_uv_bounding_sphere3(ref_geom, nodes, center, &radius),
          "bound with circle");
      center[2] = 0.0;
      iface = nodes[3] - 1;
      RSS(ref_search_insert(ref_facelift_face_search(ref_facelift, iface), cell,
                            center, scale * radius),
          "ins");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_facelift_create(REF_FACELIFT *ref_facelift_ptr,
                               REF_GRID freeable_ref_grid, REF_BOOL direct) {
  REF_FACELIFT ref_facelift;
  REF_INT n;

  ref_malloc(*ref_facelift_ptr, 1, REF_FACELIFT_STRUCT);

  ref_facelift = *ref_facelift_ptr;

  ref_facelift_grid(ref_facelift) = freeable_ref_grid;

  ref_facelift_edg(ref_facelift) =
      ref_grid_edg(ref_facelift_grid(ref_facelift));
  ref_facelift_tri(ref_facelift) =
      ref_grid_tri(ref_facelift_grid(ref_facelift));

  if (ref_cell_n(ref_grid_ed2(ref_facelift_grid(ref_facelift))) > 0) {
    ref_facelift_edg(ref_facelift) =
        ref_grid_ed2(ref_facelift_grid(ref_facelift));
  }
  if (ref_cell_n(ref_grid_ed3(ref_facelift_grid(ref_facelift))) > 0) {
    ref_facelift_edg(ref_facelift) =
        ref_grid_ed3(ref_facelift_grid(ref_facelift));
  }
  if (ref_cell_n(ref_grid_tr2(ref_facelift_grid(ref_facelift))) > 0) {
    ref_facelift_tri(ref_facelift) =
        ref_grid_tr2(ref_facelift_grid(ref_facelift));
  }
  if (ref_cell_n(ref_grid_tr3(ref_facelift_grid(ref_facelift))) > 0) {
    ref_facelift_tri(ref_facelift) =
        ref_grid_tr3(ref_facelift_grid(ref_facelift));
  }

  if (direct) {
    ref_facelift->displacement = NULL;
    ref_facelift->strong_bc = NULL;
  } else {
    n = ref_geom_max(ref_facelift_geom(ref_facelift));
    ref_malloc_init(ref_facelift->displacement, 3 * n, REF_DBL, 0.0);
    ref_malloc_init(ref_facelift->strong_bc, n, REF_BOOL, REF_FALSE);
  }

  RSS(ref_facelift_cache_search(ref_facelift), "cache tri uv");

  return REF_SUCCESS;
}

REF_STATUS ref_facelift_free(REF_FACELIFT ref_facelift) {
  REF_INT i;
  if (NULL == (void *)ref_facelift) return REF_NULL;

  if (NULL != ref_facelift->edge_search) {
    for (i = 0; i < ref_facelift_geom(ref_facelift)->nedge; i++)
      ref_search_free(ref_facelift_edge_search(ref_facelift, i));
    ref_free(ref_facelift->edge_search);
  }
  if (NULL != ref_facelift->face_search) {
    for (i = 0; i < ref_facelift_geom(ref_facelift)->nface; i++)
      ref_search_free(ref_facelift_face_search(ref_facelift, i));
    ref_free(ref_facelift->face_search);
  }
  ref_free(ref_facelift->strong_bc);
  ref_free(ref_facelift->displacement);
  ref_grid_free(ref_facelift_grid(ref_facelift));
  /* geom is a pointer to the orig */
  ref_free(ref_facelift);

  return REF_SUCCESS;
}

REF_STATUS ref_facelift_tattle(REF_GRID ref_grid, REF_INT node) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_FACELIFT ref_facelift = ref_geom_facelift(ref_geom);
  REF_INT item, geom;
  REF_INT type, id;
  REF_DBL params[2];
  REF_DBL bary[3];
  REF_INT cell;
  REF_DBL xyz[3];

  if (NULL == ref_facelift) return REF_SUCCESS;
  if (!ref_facelift_direct(ref_facelift)) return REF_SUCCESS;

  printf(" tattle on node = %d %f %f %f\n", node,
         ref_node_xyz(ref_node, 0, node), ref_node_xyz(ref_node, 1, node),
         ref_node_xyz(ref_node, 2, node));
  each_ref_adj_node_item_with_ref(ref_geom_adj(ref_geom), node, item, geom) {
    type = ref_geom_type(ref_geom, geom);
    id = ref_geom_id(ref_geom, geom);
    switch (type) {
      case REF_GEOM_EDGE:
        params[0] = ref_geom_param(ref_geom, 0, geom);
        RSS(ref_facelift_eval_at(ref_facelift, type, id, params, xyz, NULL),
            "eval");
        RSS(ref_facelift_enclosing(ref_facelift, type, id, params, &cell, bary),
            "enclose");
        printf("edg id %d xyz %f %f %f bary %f %f\n", id, xyz[0], xyz[1],
               xyz[2], bary[0], bary[1]);
        break;
      case REF_GEOM_FACE:
        params[0] = ref_geom_param(ref_geom, 0, geom);
        params[1] = ref_geom_param(ref_geom, 1, geom);
        RSS(ref_facelift_eval_at(ref_facelift, type, id, params, xyz, NULL),
            "eval");
        RSS(ref_facelift_enclosing(ref_facelift, type, id, params, &cell, bary),
            "enclose");
        printf("tri id %d xyz %f %f %f bary %f %f %f\n", id, xyz[0], xyz[1],
               xyz[2], bary[0], bary[1], bary[2]);
        break;
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_solve_face(REF_FACELIFT ref_facelift) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_CELL ref_cell = ref_facelift_tri(ref_facelift);
  REF_INT center_geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, cell_node, other_geom, center;
  REF_DBL disp[3], hits;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  /* psudo laplace */
  each_ref_geom_face(ref_geom, center_geom) {
    if (!ref_facelift_strong_bc(ref_facelift, center_geom)) {
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
            disp[0] +=
                0.5 * ref_facelift_displacement(ref_facelift, 0, other_geom);
            disp[1] +=
                0.5 * ref_facelift_displacement(ref_facelift, 1, other_geom);
            disp[2] +=
                0.5 * ref_facelift_displacement(ref_facelift, 2, other_geom);
          }
        }
        if (hits > 0.1) {
          ref_facelift_displacement(ref_facelift, 0, center_geom) =
              disp[0] / hits;
          ref_facelift_displacement(ref_facelift, 1, center_geom) =
              disp[1] / hits;
          ref_facelift_displacement(ref_facelift, 2, center_geom) =
              disp[2] / hits;
        }
      }
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_initialize_face(REF_FACELIFT ref_facelift,
                                               REF_GEOM ref_geom) {
  REF_DBL edge_xyz[3], face_xyz[3];
  REF_INT i, edge_geom, face_geom, node, item;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  /* also displace edge and face to geom nodes */

  each_ref_geom_edge(ref_geom, edge_geom) {
    RSS(ref_egads_eval(ref_geom, edge_geom, edge_xyz, NULL), "eval edge");
    node = ref_geom_node(ref_geom, edge_geom);
    each_ref_geom_having_node(ref_geom, node, item, face_geom) {
      if (REF_GEOM_FACE == ref_geom_type(ref_geom, face_geom) &&
          !ref_facelift_strong_bc(ref_facelift, face_geom)) {
        RSS(ref_egads_eval(ref_geom, face_geom, face_xyz, NULL), "eval face");
        ref_facelift_strong_bc(ref_facelift, face_geom) = REF_TRUE;
        for (i = 0; i < 3; i++) {
          ref_facelift_displacement(ref_facelift, i, face_geom) =
              edge_xyz[i] - face_xyz[i];
        }
      }
    }
  }

  for (i = 0; i < 1000; i++)
    RSS(ref_facelift_solve_face(ref_facelift), "solve");

  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_solve_edge(REF_FACELIFT ref_facelift) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_CELL ref_cell = ref_facelift_edg(ref_facelift);
  REF_INT center_geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, cell_node, other_geom, center;
  REF_DBL disp[3], hits;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  /* psudo laplace */
  each_ref_geom_edge(ref_geom, center_geom) {
    if (!ref_facelift_strong_bc(ref_facelift, center_geom)) {
      center = ref_geom_node(ref_geom, center_geom);
      hits = 0.0;
      disp[0] = 0.0;
      disp[1] = 0.0;
      disp[2] = 0.0;
      each_ref_cell_having_node(ref_cell, center, item, cell) {
        RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
        each_ref_cell_cell_node(ref_cell, cell_node) {
          if (nodes[cell_node] != center) {
            RSS(ref_geom_find(ref_geom, nodes[cell_node], REF_GEOM_EDGE,
                              ref_geom_id(ref_geom, center_geom), &other_geom),
                "other geom");
            hits += 0.5;
            disp[0] +=
                0.5 * ref_facelift_displacement(ref_facelift, 0, other_geom);
            disp[1] +=
                0.5 * ref_facelift_displacement(ref_facelift, 1, other_geom);
            disp[2] +=
                0.5 * ref_facelift_displacement(ref_facelift, 2, other_geom);
          }
        }
        if (hits > 0.1) {
          ref_facelift_displacement(ref_facelift, 0, center_geom) =
              disp[0] / hits;
          ref_facelift_displacement(ref_facelift, 1, center_geom) =
              disp[1] / hits;
          ref_facelift_displacement(ref_facelift, 2, center_geom) =
              disp[2] / hits;
        }
      }
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_initialize_edge(REF_FACELIFT ref_facelift) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_DBL node_xyz[3], edge_xyz[3];
  REF_INT i, node_geom, edge_geom, node, item;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  /* also displace edge and face to geom nodes */

  each_ref_geom_node(ref_geom, node_geom) {
    RSS(ref_egads_eval(ref_geom, node_geom, node_xyz, NULL), "eval node");
    node = ref_geom_node(ref_geom, node_geom);
    each_ref_geom_having_node(ref_geom, node, item, edge_geom) {
      if (REF_GEOM_EDGE == ref_geom_type(ref_geom, edge_geom) &&
          !ref_facelift_strong_bc(ref_facelift, edge_geom)) {
        RSS(ref_egads_eval(ref_geom, edge_geom, edge_xyz, NULL), "eval face");
        ref_facelift_strong_bc(ref_facelift, edge_geom) = REF_TRUE;
        for (i = 0; i < 3; i++) {
          ref_facelift_displacement(ref_facelift, i, edge_geom) =
              node_xyz[i] - edge_xyz[i];
        }
      }
    }
  }

  for (i = 0; i < 100; i++) RSS(ref_facelift_solve_edge(ref_facelift), "solve");

  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_apply(REF_FACELIFT ref_facelift) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_NODE ref_node = ref_grid_node(ref_facelift_grid(ref_facelift));
  REF_DBL xyz[3];
  REF_INT geom;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  each_ref_geom_face(ref_geom, geom) {
    RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, ref_geom_id(ref_geom, geom),
                          &(ref_geom_param(ref_geom, 0, geom)), xyz, NULL),
        "eval at");
    xyz[0] += ref_facelift_displacement(ref_facelift, 0, geom);
    xyz[1] += ref_facelift_displacement(ref_facelift, 1, geom);
    xyz[2] += ref_facelift_displacement(ref_facelift, 2, geom);
    ref_node_xyz(ref_node, 0, ref_geom_node(ref_geom, geom)) = xyz[0];
    ref_node_xyz(ref_node, 1, ref_geom_node(ref_geom, geom)) = xyz[1];
    ref_node_xyz(ref_node, 2, ref_geom_node(ref_geom, geom)) = xyz[2];
  }
  return REF_SUCCESS;
}

REF_STATUS ref_facelift_attach(REF_GRID ref_grid) {
  REF_GRID freeable_ref_grid;
  REF_FACELIFT ref_facelift;
  RSS(ref_grid_deep_copy(&freeable_ref_grid, ref_grid), "deep copy");
  RSS(ref_facelift_create(&ref_facelift, freeable_ref_grid, REF_FALSE),
      "create");
  ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
  RSS(ref_facelift_initialize_edge(ref_facelift), "init disp");
  RSS(ref_facelift_initialize_face(ref_facelift, ref_grid_geom(ref_grid)),
      "init disp");
  RSS(ref_facelift_apply(ref_facelift), "apply");
  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_infer_displacement(REF_FACELIFT ref_facelift) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_NODE ref_node = ref_grid_node(ref_facelift_grid(ref_facelift));
  REF_DBL geom_xyz[3];
  REF_INT i, geom, node;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  each_ref_geom(ref_geom, geom) {
    RSS(ref_egads_eval(ref_geom, geom, geom_xyz, NULL), "eval geom");
    node = ref_geom_node(ref_geom, geom);
    for (i = 0; i < 3; i++) {
      ref_facelift_displacement(ref_facelift, i, geom) =
          ref_node_xyz(ref_node, i, node) - geom_xyz[i];
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_facelift_import(REF_GRID ref_grid, const char *filename) {
  REF_GRID freeable_ref_grid;
  REF_FACELIFT ref_facelift;
  REF_GEOM ref_geom;
  RSS(ref_import_by_extension(&freeable_ref_grid, ref_grid_mpi(ref_grid),
                              filename),
      "import");
  ref_geom = ref_grid_geom(freeable_ref_grid);
  RSS(ref_geom_share_context(ref_geom, ref_grid_geom(ref_grid)),
      "share context");
  if (ref_geom_model_loaded(ref_geom)) {
    RSS(ref_egads_mark_jump_degen(freeable_ref_grid),
        "T and UV jumps; UV degen");
    RSS(ref_geom_verify_topo(freeable_ref_grid), "geom topo");
    RSS(ref_geom_verify_param(freeable_ref_grid), "geom param");
  }
  RSS(ref_facelift_create(&ref_facelift, freeable_ref_grid, REF_FALSE),
      "create");
  ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
  RSS(ref_facelift_infer_displacement(ref_facelift), "infer displacement");
  return REF_SUCCESS;
}

REF_STATUS ref_facelift_surrogate(REF_GRID ref_grid, const char *filename) {
  REF_GRID freeable_ref_grid;
  REF_FACELIFT ref_facelift;
  REF_GEOM ref_geom;
  RSS(ref_import_by_extension(&freeable_ref_grid, ref_grid_mpi(ref_grid),
                              filename),
      "import");
  ref_geom = ref_grid_geom(freeable_ref_grid);
  RSS(ref_geom_share_context(ref_geom, ref_grid_geom(ref_grid)),
      "share context");
  if (ref_geom_model_loaded(ref_geom)) {
    RSS(ref_egads_mark_jump_degen(freeable_ref_grid),
        "T and UV jumps; UV degen");
    RSS(ref_geom_verify_topo(freeable_ref_grid), "geom topo");
    RSS(ref_geom_verify_param(freeable_ref_grid), "geom param");
  }
  RSS(ref_facelift_create(&ref_facelift, freeable_ref_grid, REF_TRUE),
      "create");
  ref_geom_facelift(ref_grid_geom(ref_grid)) = ref_facelift;
  return REF_SUCCESS;
}

REF_STATUS ref_facelift_enclosing(REF_FACELIFT ref_facelift, REF_INT type,
                                  REF_INT id, REF_DBL *param, REF_INT *cell,
                                  REF_DBL *bary) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
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

  if (REF_GEOM_EDGE == type) {
    RSS(ref_list_create(&ref_list), "create list");
    ref_search = ref_facelift_edge_search(ref_facelift, id - 1);
    ref_cell = ref_facelift_edg(ref_facelift);
    parampad[0] = param[0];
    parampad[1] = 0.0;
    parampad[2] = 0.0;

    RSS(ref_search_touching(ref_search, ref_list, parampad, fuzz), "touching");
    if (0 >= ref_list_n(ref_list)) {
      RSS(ref_list_free(ref_list), "free list");
      return REF_SUCCESS;
    }
    best_candidate = REF_EMPTY;
    best_bary = -999.0;
    each_ref_list_item(ref_list, item) {
      candidate = ref_list_value(ref_list, item);
      RSS(ref_cell_nodes(ref_cell, candidate, nodes), "cell");
      nodes[2] = nodes[ref_cell_id_index(ref_cell)]; /* expects P1 */
      RSS(ref_geom_bary2(ref_geom, nodes, param[0], current_bary), "bary");
      min_bary = MIN(current_bary[0], current_bary[1]);
      if (REF_EMPTY == best_candidate || min_bary > best_bary) {
        best_candidate = candidate;
        best_bary = min_bary;
      }
    }

    RUS(REF_EMPTY, best_candidate, "failed to find cell");

    *cell = best_candidate;
    RSS(ref_cell_nodes(ref_cell, best_candidate, nodes), "cell");
    nodes[2] = nodes[ref_cell_id_index(ref_cell)]; /* expects P1 */
    RSS(ref_geom_bary2(ref_geom, nodes, param[0], bary), "bary");

    RSS(ref_list_free(ref_list), "free list");
  }

  if (REF_GEOM_FACE == type) {
    RSS(ref_list_create(&ref_list), "create list");
    ref_search = ref_facelift_face_search(ref_facelift, id - 1);
    ref_cell = ref_facelift_tri(ref_facelift);
    parampad[0] = param[0];
    parampad[1] = param[1];
    parampad[2] = 0.0;

    RSS(ref_search_touching(ref_search, ref_list, parampad, fuzz), "touching");
    if (0 >= ref_list_n(ref_list)) {
      RSS(ref_list_free(ref_list), "free list");
      return REF_SUCCESS;
    }
    best_candidate = REF_EMPTY;
    best_bary = -999.0;
    each_ref_list_item(ref_list, item) {
      candidate = ref_list_value(ref_list, item);
      RSS(ref_cell_nodes(ref_cell, candidate, nodes), "cell");
      nodes[3] = nodes[ref_cell_id_index(ref_cell)]; /* expects P1 */
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
    nodes[3] = nodes[ref_cell_id_index(ref_cell)]; /* expects P1 */
    RSS(ref_geom_bary3(ref_geom, nodes, param, bary), "bary");

    RSS(ref_list_free(ref_list), "free list");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_displacement_at(REF_FACELIFT ref_facelift,
                                               REF_INT type, REF_INT id,
                                               REF_DBL *params,
                                               REF_DBL *displacement) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL bary[3], clip[3];
  displacement[0] = 0.0;
  displacement[1] = 0.0;
  displacement[2] = 0.0;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  if (REF_GEOM_EDGE == type) {
    RSS(ref_facelift_enclosing(ref_facelift, type, id, params, &cell, bary),
        "enclose");
    if (REF_EMPTY == cell) return REF_SUCCESS;
    RSS(ref_node_clip_bary2(bary, clip), "clip edge bary");
    RSS(ref_cell_nodes(ref_facelift_edg(ref_facelift), cell, nodes), "nodes");

    RSS(ref_geom_find(ref_geom, nodes[0], REF_GEOM_EDGE, id, &geom), "find 0");
    displacement[0] +=
        clip[0] * ref_facelift_displacement(ref_facelift, 0, geom);
    displacement[1] +=
        clip[0] * ref_facelift_displacement(ref_facelift, 1, geom);
    displacement[2] +=
        clip[0] * ref_facelift_displacement(ref_facelift, 2, geom);

    RSS(ref_geom_find(ref_geom, nodes[1], REF_GEOM_EDGE, id, &geom), "find 1");
    displacement[0] +=
        clip[1] * ref_facelift_displacement(ref_facelift, 0, geom);
    displacement[1] +=
        clip[1] * ref_facelift_displacement(ref_facelift, 1, geom);
    displacement[2] +=
        clip[1] * ref_facelift_displacement(ref_facelift, 2, geom);
  }
  if (REF_GEOM_FACE == type) {
    RSS(ref_facelift_enclosing(ref_facelift, type, id, params, &cell, bary),
        "enclose");
    if (REF_EMPTY == cell) return REF_SUCCESS;
    RSS(ref_node_clip_bary3(bary, clip), "clip face bary");
    RSS(ref_cell_nodes(ref_facelift_tri(ref_facelift), cell, nodes), "nodes");

    RSS(ref_geom_find(ref_geom, nodes[0], REF_GEOM_FACE, id, &geom), "find 0");
    displacement[0] +=
        clip[0] * ref_facelift_displacement(ref_facelift, 0, geom);
    displacement[1] +=
        clip[0] * ref_facelift_displacement(ref_facelift, 1, geom);
    displacement[2] +=
        clip[0] * ref_facelift_displacement(ref_facelift, 2, geom);

    RSS(ref_geom_find(ref_geom, nodes[1], REF_GEOM_FACE, id, &geom), "find 1");
    displacement[0] +=
        clip[1] * ref_facelift_displacement(ref_facelift, 0, geom);
    displacement[1] +=
        clip[1] * ref_facelift_displacement(ref_facelift, 1, geom);
    displacement[2] +=
        clip[1] * ref_facelift_displacement(ref_facelift, 2, geom);

    RSS(ref_geom_find(ref_geom, nodes[2], REF_GEOM_FACE, id, &geom), "find 2");
    displacement[0] +=
        clip[2] * ref_facelift_displacement(ref_facelift, 0, geom);
    displacement[1] +=
        clip[2] * ref_facelift_displacement(ref_facelift, 1, geom);
    displacement[2] +=
        clip[2] * ref_facelift_displacement(ref_facelift, 2, geom);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_facelift_eval_at(REF_FACELIFT ref_facelift, REF_INT type,
                                REF_INT id, REF_DBL *params, REF_DBL *xyz,
                                REF_DBL *dxyz_dtuv) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_DBL displacement[3];
  RSS(ref_egads_eval_at(ref_geom, type, id, params, xyz, dxyz_dtuv),
      "egads eval");
  if (ref_facelift_direct(ref_facelift)) {
    REF_NODE ref_node = ref_grid_node(ref_facelift_grid(ref_facelift));
    REF_INT i, cell_node, cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_DBL bary[3], clip[3];
    REF_DBL shape[REF_CELL_MAX_NODE_PER];
    REF_CELL ref_cell = NULL;
    RSS(ref_facelift_enclosing(ref_facelift, type, id, params, &cell, bary),
        "enclose");
    if (REF_EMPTY == cell) return REF_SUCCESS;
    if (REF_GEOM_EDGE == type) {
      ref_cell = ref_facelift_edg(ref_facelift);
      RSS(ref_node_clip_bary2(bary, clip), "clip edge bary");
    }
    if (REF_GEOM_FACE == type) {
      ref_cell = ref_facelift_tri(ref_facelift);
      RSS(ref_node_clip_bary3(bary, clip), "clip face bary");
    }
    RSS(ref_cell_shape(ref_cell, clip, shape), "shape");
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    for (i = 0; i < 3; i++) {
      xyz[i] = 0.0;
      each_ref_cell_cell_node(ref_cell, cell_node) {
        xyz[i] +=
            shape[cell_node] * ref_node_xyz(ref_node, i, nodes[cell_node]);
      }
    }
  } else {
    RSS(ref_facelift_displacement_at(ref_facelift, type, id, params,
                                     displacement),
        "facelift displacement");
    xyz[0] += displacement[0];
    xyz[1] += displacement[1];
    xyz[2] += displacement[2];
  }
  return REF_SUCCESS;
}

REF_STATUS ref_facelift_inverse_eval(REF_FACELIFT ref_facelift, REF_INT type,
                                     REF_INT id, REF_DBL *xyz, REF_DBL *param) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_DBL geom_xyz[3];
  REF_DBL displacement[3];

  geom_xyz[0] = xyz[0];
  geom_xyz[1] = xyz[1];
  geom_xyz[2] = xyz[2];
  RSS(ref_egads_inverse_eval(ref_geom, type, id, geom_xyz, param),
      "inv eval before");

  if (ref_facelift_direct(ref_facelift)) {
  } else {
    RSS(ref_facelift_displacement_at(ref_facelift, type, id, param,
                                     displacement),
        "facelift displacement");

    geom_xyz[0] -= displacement[0];
    geom_xyz[1] -= displacement[1];
    geom_xyz[2] -= displacement[2];
    RSS(ref_egads_inverse_eval(ref_geom, type, id, geom_xyz, param),
        "inv eval before");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_facelift_edge_face_uv(REF_FACELIFT ref_facelift, REF_INT edgeid,
                                     REF_INT faceid, REF_INT sense, REF_DBL t,
                                     REF_DBL *uv) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_CELL ref_cell = ref_facelift_edg(ref_facelift);
  REF_INT i, cell_node, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL bary[2], clip[2], faceuv[2];
  REF_DBL shape[REF_CELL_MAX_NODE_PER];
  RSS(ref_egads_edge_face_uv(ref_geom, edgeid, faceid, sense, t, uv),
      "edge uv");

  if (ref_facelift_direct(ref_facelift)) {
    RAS(0 == sense, "implement sense != 0 for uv jumps");
    RSS(ref_facelift_enclosing(ref_facelift, REF_GEOM_EDGE, edgeid, &t, &cell,
                               bary),
        "enclose");
    if (REF_EMPTY == cell) return REF_SUCCESS;
    RSS(ref_node_clip_bary2(bary, clip), "clip edge bary");
    RSS(ref_cell_shape(ref_cell, clip, shape), "shape");
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    for (i = 0; i < 2; i++) {
      uv[i] = 0.0;
      each_ref_cell_cell_node(ref_cell, cell_node) {
        RSS(ref_geom_tuv(ref_geom, nodes[cell_node], REF_GEOM_FACE, faceid,
                         faceuv),
            "face uv");
        uv[i] += shape[cell_node] * faceuv[i];
      }
    }

  } else {
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_edge_tec_zone(REF_FACELIFT ref_facelift,
                                             REF_INT id, FILE *file) {
  REF_GRID ref_grid = ref_facelift_grid(ref_facelift);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_facelift_edg(ref_facelift);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, nedg, sens;
  REF_INT jump_geom = REF_EMPTY;
  REF_DBL *t, tvalue;
  REF_DBL radius, normal[3], xyz[3];

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

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
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_edge_curvature(ref_geom, geom, &radius, normal), "curve");
      radius = ABS(radius);
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_EDGE, id, &(t[item]), xyz, NULL),
          "eval at");
    }
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
            xyz[0], xyz[1], xyz[2], ref_facelift_distance(ref_facelift, geom),
            t[item], 0.0, ref_facelift_displacement(ref_facelift, 0, geom),
            ref_facelift_displacement(ref_facelift, 1, geom),
            ref_facelift_displacement(ref_facelift, 2, geom));
  }
  if (REF_EMPTY != jump_geom) {
    node = ref_geom_node(ref_geom, jump_geom);
    radius = 0;
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
    }
    node = ref_geom_node(ref_geom, jump_geom);
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
            xyz[0], xyz[1], xyz[2],
            ref_facelift_distance(ref_facelift, jump_geom), t[nnode - 1], 0.0,
            ref_facelift_displacement(ref_facelift, 0, jump_geom),
            ref_facelift_displacement(ref_facelift, 1, jump_geom),
            ref_facelift_displacement(ref_facelift, 2, jump_geom));
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

static REF_STATUS ref_facelift_face_tec_zone(REF_FACELIFT ref_facelift,
                                             REF_INT id, FILE *file) {
  REF_GRID ref_grid = ref_facelift_grid(ref_facelift);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_facelift_tri(ref_facelift);
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_DICT ref_dict, ref_dict_jump, ref_dict_degen;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, nnode_sens0, nnode_degen, ntri;
  REF_INT sens;
  REF_DBL *uv, param[2];
  REF_DBL xyz[3];

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

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
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, id, &(uv[2 * item]), xyz,
                            NULL),
          "eval at");
    }
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
            xyz[0], xyz[1], xyz[2], ref_facelift_distance(ref_facelift, geom),
            uv[0 + 2 * item], uv[1 + 2 * item],
            ref_facelift_displacement(ref_facelift, 0, geom),
            ref_facelift_displacement(ref_facelift, 1, geom),
            ref_facelift_displacement(ref_facelift, 2, geom));
  }
  each_ref_dict_key_value(ref_dict_jump, item, node, geom) {
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, id,
                            &(uv[2 * (nnode_sens0 + item)]), xyz, NULL),
          "eval at");
    }
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
            xyz[0], xyz[1], xyz[2], ref_facelift_distance(ref_facelift, geom),
            uv[0 + 2 * (nnode_sens0 + item)], uv[1 + 2 * (nnode_sens0 + item)],
            ref_facelift_displacement(ref_facelift, 0, geom),
            ref_facelift_displacement(ref_facelift, 1, geom),
            ref_facelift_displacement(ref_facelift, 2, geom));
  }
  each_ref_dict_key_value(ref_dict_degen, item, cell, node) {
    RSS(ref_geom_find(ref_geom, node, REF_GEOM_FACE, id, &geom), "find");
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    if (ref_geom_model_loaded(ref_geom)) {
      RSS(ref_egads_eval_at(ref_geom, REF_GEOM_FACE, id,
                            &(uv[2 * (nnode_degen + item)]), xyz, NULL),
          "eval at");
    }
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
            xyz[0], xyz[1], xyz[2], ref_facelift_distance(ref_facelift, geom),
            uv[0 + 2 * (nnode_degen + item)], uv[1 + 2 * (nnode_degen + item)],
            ref_facelift_displacement(ref_facelift, 0, geom),
            ref_facelift_displacement(ref_facelift, 1, geom),
            ref_facelift_displacement(ref_facelift, 2, geom));
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

REF_STATUS ref_facelift_tec(REF_FACELIFT ref_facelift, const char *filename) {
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  FILE *file;
  REF_INT geom, id, min_id, max_id;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"refine watertight cad displacement\"\n");
  fprintf(
      file,
      "variables = \"x\" \"y\" \"z\" \"d\" \"u\" \"v\" \"dx\" \"dy\" \"dz\"\n");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_edge(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_facelift_edge_tec_zone(ref_facelift, id, file), "tec edge");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_facelift_face_tec_zone(ref_facelift, id, file), "tec face");

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_facelift_max_distance(REF_FACELIFT ref_facelift,
                                     REF_DBL *distance) {
  REF_GRID ref_grid = ref_facelift_grid(ref_facelift);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_facelift_geom(ref_facelift);
  REF_INT geom, node;

  RAS(!ref_facelift_direct(ref_facelift), "requires indirect facelift");

  each_ref_node_valid_node(ref_node, node) { distance[node] = 0.0; }

  each_ref_geom(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    distance[node] =
        MAX(distance[node], ref_facelift_distance(ref_facelift, geom));
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_complexity(REF_DBL *metric, REF_GRID ref_grid,
                                          REF_DBL *complexity) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell_node, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL area, det, m2[3], r[3], s[3], n[3];
  ref_cell = ref_grid_tri(ref_grid);
  *complexity = 0.0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      if (ref_node_owned(ref_node, nodes[cell_node])) {
        RSS(ref_recon_rsn(ref_grid, nodes[cell_node], r, s, n), "rsn");
        RSS(ref_matrix_extract2(&(metric[6 * nodes[cell_node]]), r, s, m2),
            "extract");
        RSS(ref_matrix_det_m2(m2, &det), "2x2 det");
        if (det > 0.0) {
          (*complexity) +=
              sqrt(det) * area / ((REF_DBL)ref_cell_node_per(ref_cell));
        }
      }
    }
  }
  RSS(ref_mpi_allsum(ref_grid_mpi(ref_grid), complexity, 1, REF_DBL_TYPE),
      "dbl sum");

  return REF_SUCCESS;
}

static REF_STATUS ref_facelift_gradation_at_complexity(REF_DBL *metric,
                                                       REF_GRID ref_grid,
                                                       REF_DBL gradation,
                                                       REF_DBL complexity) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT relaxations;
  REF_DBL current_complexity;
  REF_DBL complexity_scale;
  REF_INT node, i;

  complexity_scale = 1.0; /* "2D" along surf */

  for (relaxations = 0; relaxations < 20; relaxations++) {
    RSS(ref_facelift_complexity(metric, ref_grid, &current_complexity), "cmp");
    if (!ref_math_divisible(complexity, current_complexity)) {
      return REF_DIV_ZERO;
    }
    each_ref_node_valid_node(ref_node, node) {
      for (i = 0; i < 6; i++) {
        metric[i + 6 * node] *=
            pow(complexity / current_complexity, complexity_scale);
      }
    }
    if (gradation < 1.0) {
      RSS(ref_metric_mixed_space_gradation(metric, ref_grid, -1.0, -1.0),
          "gradation");
    } else {
      RSS(ref_metric_metric_space_gradation(metric, ref_grid, gradation),
          "gradation");
    }
  }
  RSS(ref_facelift_complexity(metric, ref_grid, &current_complexity), "cmp");
  if (!ref_math_divisible(complexity, current_complexity)) {
    return REF_DIV_ZERO;
  }
  each_ref_node_valid_node(ref_node, node) {
    for (i = 0; i < 6; i++) {
      metric[i + 6 * node] *=
          pow(complexity / current_complexity, complexity_scale);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_facelift_multiscale(REF_GRID ref_grid,
                                   REF_DBL target_complexity) {
  REF_FACELIFT ref_facelift;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_DBL *hess, *metric;
  REF_INT node, i;
  REF_INT dimension = 2, p_norm = 2;
  REF_INT gradation = -1.0;
  REF_DBL det, exponent, complexity, area, complexity_scale;
  REF_DBL diag_system2[6];
  REF_DBL diag_system[12];
  REF_DBL m[6], combined[6];
  REF_DBL r[3], s[3], n[3];
  REF_BOOL verbose = REF_FALSE;
  REF_BOOL steps = REF_FALSE;
  REF_INT cell_node, cell, nodes[REF_CELL_MAX_SIZE_PER];

  /* reset facelift to match grid */
  ref_facelift = ref_geom_facelift(ref_grid_geom(ref_grid));
  if (NULL != ref_facelift) ref_facelift_free(ref_facelift);
  ref_facelift = NULL;
  RSS(ref_facelift_attach(ref_grid), "attach");
  ref_facelift = ref_geom_facelift(ref_grid_geom(ref_grid));
  ref_cell = ref_facelift_tri(ref_facelift);

  exponent = -1.0 / ((REF_DBL)(2 * p_norm + dimension));

  ref_malloc(metric, 6 * ref_node_max(ref_node), REF_DBL);
  ref_malloc_init(hess, 3 * ref_node_max(ref_grid_node(ref_grid)), REF_DBL,
                  0.0);
  RSS(ref_recon_rsn_hess_face(ref_grid, hess), "rsn");
  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_matrix_det_m2(&(hess[3 * node]), &det), "2x2 det");
    if (det > 0.0) { /* local scaling */
      for (i = 0; i < 3; i++) hess[i + 3 * node] *= pow(det, exponent);
    }
  }

  complexity = 0.0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      if (ref_node_owned(ref_node, nodes[cell_node])) {
        RSS(ref_matrix_det_m2(&(hess[3 * nodes[cell_node]]), &det), "2x2 det");
        if (det > 0.0) {
          complexity +=
              sqrt(det) * area / ((REF_DBL)ref_cell_node_per(ref_cell));
        }
      }
    }
  }

  if (!ref_math_divisible(target_complexity, complexity)) {
    RSS(REF_DIV_ZERO, "zero compelxity");
  }

  complexity_scale = 1.0;
  each_ref_node_valid_node(ref_node, node) { /* global scaling */
    for (i = 0; i < 3; i++)
      hess[i + 3 * node] *=
          pow(target_complexity / complexity, complexity_scale);
  }

  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_recon_rsn(ref_grid, node, r, s, n), "rsn");
    RSS(ref_matrix_diag_m2(&(hess[3 * node]), diag_system2), "decomp");
    ref_matrix_eig(diag_system, 0) = ref_matrix_eig2(diag_system2, 0);
    for (i = 0; i < 3; i++) {
      ref_matrix_vec(diag_system, i, 0) =
          ref_matrix_vec2(diag_system2, 0, 0) * r[i] +
          ref_matrix_vec2(diag_system2, 1, 0) * s[i];
    }
    ref_matrix_eig(diag_system, 1) = ref_matrix_eig2(diag_system2, 1);
    for (i = 0; i < 3; i++) {
      ref_matrix_vec(diag_system, i, 1) =
          ref_matrix_vec2(diag_system2, 0, 1) * r[i] +
          ref_matrix_vec2(diag_system2, 1, 1) * s[i];
    }
    ref_matrix_eig(diag_system, 2) =
        MAX(ref_matrix_eig2(diag_system2, 0), ref_matrix_eig2(diag_system2, 1));
    for (i = 0; i < 3; i++) ref_matrix_vec(diag_system, i, 2) = n[i];

    RSS(ref_matrix_form_m(diag_system, &(metric[6 * node])), "form m");
    if (verbose) {
      ref_matrix_show_diag_sys(diag_system);
      ref_metric_show(&(metric[6 * node]));
      printf("n %f %f %f\n", n[0], n[1], n[2]);
      printf("v %f %f %f\n", diag_system[0 + 9], diag_system[1 + 9],
             diag_system[2 + 9]);
    }
  }
  ref_free(hess);

  RSS(ref_recon_roundoff_limit(metric, ref_grid), "floor eigs above zero");

  if (steps) {
    RSS(ref_metric_to_node(metric,
                           ref_grid_node(ref_facelift_grid(ref_facelift))),
        "to");
    RSS(ref_export_tec_metric_ellipse(ref_facelift_grid(ref_facelift),
                                      "ref_facelift_raw"),
        "al");
  }

  RSS(ref_facelift_gradation_at_complexity(metric, ref_grid, gradation,
                                           target_complexity),
      "gradation at complexity");

  each_ref_node_valid_node(ref_node, node) {
    RSS(ref_node_metric_get(ref_node, node, m), "curve metric");
    RSS(ref_matrix_intersect(&(metric[6 * node]), m, combined), "intersect");
    if (verbose) {
      ref_metric_show(&(metric[6 * node]));
      ref_metric_show(m);
      ref_metric_show(combined);
      printf("\n");
    }
    for (i = 0; i < 6; i++) metric[i + 6 * node] = combined[i];
  }
  RSS(ref_metric_to_node(metric, ref_grid_node(ref_grid)), "to");

  if (steps) {
    RSS(ref_metric_to_node(metric,
                           ref_grid_node(ref_facelift_grid(ref_facelift))),
        "to");
    RSS(ref_export_tec_metric_ellipse(ref_facelift_grid(ref_facelift),
                                      "ref_facelift_grad"),
        "al");
  }

  ref_free(metric);

  return REF_SUCCESS;
}
