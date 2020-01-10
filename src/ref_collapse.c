
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

#include "ref_collapse.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_adapt.h"
#include "ref_cavity.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_gather.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_mpi.h"
#include "ref_sort.h"

#define MAX_CELL_COLLAPSE (100)
#define MAX_NODE_LIST (1000)

REF_STATUS ref_collapse_diagnostics(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_EDGE ref_edge;
  REF_INT edge, node0, node1;
  REF_DBL edge_ratio;
  REF_BOOL has_side;
  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");
  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);
    RSS(ref_node_ratio(ref_node, node0, node1, &edge_ratio), "ratio");
    if (edge_ratio < 0.1 + ref_grid_adapt(ref_grid, post_min_ratio)) {
      RSS(ref_cell_has_side(ref_cell, node0, node1, &has_side), "side");
      printf("ratio %f bound nodes %d %d tri edge %d\n", edge_ratio,
             !ref_cell_node_empty(ref_cell, node0),
             !ref_cell_node_empty(ref_cell, node1), has_side);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_pass(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *order;
  REF_INT ntarget, *target, *node2target;
  REF_INT node, node0, node1;
  REF_INT i, edge;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL edge_ratio;

  if (ref_grid_surf(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  ref_malloc_init(ratio, ref_node_max(ref_node), REF_DBL,
                  2.0 * ref_grid_adapt(ref_grid, collapse_ratio));

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);
    RSS(ref_node_ratio(ref_node, node0, node1, &edge_ratio), "ratio");
    ratio[node0] = MIN(ratio[node0], edge_ratio);
    ratio[node1] = MIN(ratio[node1], edge_ratio);
  }

  ref_malloc(target, ref_node_n(ref_node), REF_INT);
  ref_malloc_init(node2target, ref_node_max(ref_node), REF_INT, REF_EMPTY);

  ntarget = 0;
  for (node = 0; node < ref_node_max(ref_node); node++)
    if (ratio[node] < ref_grid_adapt(ref_grid, collapse_ratio)) {
      node2target[node] = ntarget;
      target[ntarget] = node;
      ratio[ntarget] = ratio[node];
      ntarget++;
    }

  ref_malloc(order, ntarget, REF_INT);

  RSS(ref_sort_heap_dbl(ntarget, ratio, order), "sort lengths");

  for (i = 0; i < ntarget; i++) {
    if (ratio[order[i]] > ref_grid_adapt(ref_grid, collapse_ratio)) continue;
    node1 = target[order[i]];
    if (!ref_node_valid(ref_node, node1)) continue;
    RSS(ref_collapse_to_remove_node1(ref_grid, &node0, node1), "collapse rm");
    if (!ref_node_valid(ref_node, node1)) {
      ref_node_age(ref_node, node0) = 0;
      each_ref_cell_having_node(ref_cell, node0, item, cell) {
        RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
        for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
          if (REF_EMPTY != node2target[nodes[node]]) {
            ratio[node2target[nodes[node]]] =
                2.0 * ref_grid_adapt(ref_grid, collapse_ratio);
          }
        }
      }
    }
  }

  ref_free(order);
  ref_free(node2target);
  ref_free(target);
  ref_free(ratio);

  ref_edge_free(ref_edge);

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_to_remove_node1(REF_GRID ref_grid,
                                        REF_INT *actual_node0, REF_INT node1) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT nnode, node;
  REF_INT node_to_collapse[MAX_NODE_LIST];
  REF_INT order[MAX_NODE_LIST];
  REF_DBL ratio_to_collapse[MAX_NODE_LIST];
  REF_INT node0;
  REF_BOOL allowed, local, have_geometry_support;
  REF_CAVITY ref_cavity = (REF_CAVITY)NULL;
  REF_BOOL valid_cavity;
  REF_BOOL allowed_cavity_ratio;
  REF_DBL min_del, min_add;
  REF_BOOL audit = REF_FALSE;

  *actual_node0 = REF_EMPTY;
  RAS(ref_node_valid(ref_node, node1), "node1 is invalid");

  if (ref_grid_surf(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }

  RSS(ref_cell_node_list_around(ref_cell, node1, MAX_NODE_LIST, &nnode,
                                node_to_collapse),
      "da hood");
  for (node = 0; node < nnode; node++) {
    RSS(ref_node_ratio(ref_node, node_to_collapse[node], node1,
                       &(ratio_to_collapse[node])),
        "ratio");
  }

  RSS(ref_sort_heap_dbl(nnode, ratio_to_collapse, order), "sort lengths");

  /* audit = (nnode > 0 && ratio_to_collapse[order[0]] < 0.2); */
  if (audit) {
    printf("node1 %d %f %f %f\n", node1, ref_node_xyz(ref_node, 0, node1),
           ref_node_xyz(ref_node, 1, node1), ref_node_xyz(ref_node, 2, node1));
  }

  for (node = 0; node < nnode; node++) {
    node0 = node_to_collapse[order[node]];
    if (audit)
      printf(" %d node0 %d ratio %f\n", nnode - node, node0,
             ratio_to_collapse[order[node]]);

    RSS(ref_collapse_edge_mixed(ref_grid, node0, node1, &allowed), "col mixed");
    if (!allowed && audit) printf("   mixed\n");
    if (!allowed) continue;

    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");
    if (!allowed && audit) printf("   geom\n");
    if (!allowed) continue;

    RSS(ref_collapse_edge_manifold(ref_grid, node0, node1, &allowed),
        "col manifold");
    if (!allowed && audit) printf("   manifold\n");
    if (!allowed) continue;

    RSS(ref_collapse_edge_chord_height(ref_grid, node0, node1, &allowed),
        "col edge chord height");
    if (!allowed && audit) printf("   chord\n");
    if (!allowed) continue;

    RSS(ref_collapse_edge_ratio(ref_grid, node0, node1, &allowed), "ratio");
    if (!allowed && audit) printf("   ratio\n");
    if (!allowed) continue;

    RSS(ref_collapse_surf_ratio(ref_grid, node0, node1, &allowed),
        "surf ratio");
    if (!allowed && audit) printf("   ratio (surf)\n");
    if (!allowed) continue;

    RSS(ref_geom_supported(ref_grid_geom(ref_grid), node0,
                           &have_geometry_support),
        "geom");
    if (have_geometry_support) {
      RSS(ref_collapse_edge_normdev(ref_grid, node0, node1, &allowed),
          "normdev");
      if (!allowed && audit) printf("   normdev\n");
      if (!allowed) continue;
    } else {
      RSS(ref_collapse_edge_same_normal(ref_grid, node0, node1, &allowed),
          "normal deviation");
      if (!allowed && audit) printf("   same normal\n");
      if (!allowed) continue;
    }

    RSS(ref_collapse_edge_tri_quality(ref_grid, node0, node1, &allowed),
        "tri qual");
    if (!allowed && audit) printf("   tri qual\n");
    if (!allowed) continue;

    RSS(ref_collapse_edge_tet_quality(ref_grid, node0, node1, &allowed),
        "tet qual");
    if (!allowed && audit) printf("   tet qual\n");

    RSS(ref_collapse_edge_local_cell(ref_grid, node0, node1, &local), "colloc");
    if (!local) {
      if (allowed) {
        ref_node_age(ref_node, node0)++;
        ref_node_age(ref_node, node1)++;
      }
      continue;
    }

    if (!allowed) {
      RSS(ref_cavity_create(&ref_cavity), "cav create");
      if ((REF_SUCCESS ==
           ref_cavity_form_edge_collapse(ref_cavity, ref_grid, node0, node1)) &&
          (REF_CAVITY_INCONSISTENT != ref_cavity_state(ref_cavity))) {
        RSS(ref_cavity_enlarge_visible(ref_cavity), "enlarge");
        if (REF_CAVITY_VISIBLE == ref_cavity_state(ref_cavity)) {
          RSS(ref_cavity_ratio(ref_cavity, &allowed_cavity_ratio),
              "cavity ratio");
          RSS(ref_cavity_change(ref_cavity, &min_del, &min_add),
              "cavity change");
          valid_cavity =
              allowed_cavity_ratio &&
              (min_add > ref_grid_adapt(ref_grid, collapse_quality_absolute));
          if (REF_FALSE && valid_cavity)
            printf("new %f old %f\n", min_add, min_del);
          if (valid_cavity) {
            *actual_node0 = node0;
            RSS(ref_cavity_replace(ref_cavity), "cav replace");
            RSS(ref_cavity_free(ref_cavity), "cav free");
            ref_cavity = (REF_CAVITY)NULL;
            return REF_SUCCESS;
          }
        }
        if (REF_CAVITY_PARTITION_CONSTRAINED == ref_cavity_state(ref_cavity)) {
          ref_node_age(ref_node, node0)++;
          ref_node_age(ref_node, node1)++;
        }
      }
      RSS(ref_cavity_free(ref_cavity), "cav free");
      ref_cavity = (REF_CAVITY)NULL;
      if (!allowed && audit) printf("   cav unsuccessful\n");
      continue;
    }

    *actual_node0 = node0;
    RSS(ref_collapse_edge(ref_grid, node0, node1), "col!");
    return REF_SUCCESS;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge(REF_GRID ref_grid, REF_INT node0, REF_INT node1)
/*                               keep node0,  remove node1 */
{
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT ncell, cell_in_list;
  REF_INT cell_to_collapse[MAX_CELL_COLLAPSE];

  ref_cell = ref_grid_tet(ref_grid);
  RSS(ref_cell_list_with2(ref_cell, node0, node1, MAX_CELL_COLLAPSE, &ncell,
                          cell_to_collapse),
      "lst");

  for (cell_in_list = 0; cell_in_list < ncell; cell_in_list++) {
    cell = cell_to_collapse[cell_in_list];
    RSS(ref_cell_remove(ref_cell, cell), "remove");
  }
  RSS(ref_cell_replace_node(ref_cell, node1, node0), "replace node");

  ref_cell = ref_grid_tri(ref_grid);
  RSS(ref_cell_list_with2(ref_cell, node0, node1, MAX_CELL_COLLAPSE, &ncell,
                          cell_to_collapse),
      "lst");

  for (cell_in_list = 0; cell_in_list < ncell; cell_in_list++) {
    cell = cell_to_collapse[cell_in_list];
    RSS(ref_cell_remove(ref_cell, cell), "remove");
  }
  RSS(ref_cell_replace_node(ref_cell, node1, node0), "replace node");

  ref_cell = ref_grid_edg(ref_grid);
  RSS(ref_cell_list_with2(ref_cell, node0, node1, MAX_CELL_COLLAPSE, &ncell,
                          cell_to_collapse),
      "lst");

  for (cell_in_list = 0; cell_in_list < ncell; cell_in_list++) {
    cell = cell_to_collapse[cell_in_list];
    RSS(ref_cell_remove(ref_cell, cell), "remove");
  }
  RSS(ref_cell_replace_node(ref_cell, node1, node0), "replace node");

  RSS(ref_node_remove(ref_grid_node(ref_grid), node1), "rm");
  RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), node1), "rm");

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_geometry(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_BOOL *allowed) {
  REF_CELL ref_tri = ref_grid_tri(ref_grid);
  REF_CELL ref_edg = ref_grid_edg(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT degree1;
  REF_INT ids1[3];

  REF_INT ncell;
  REF_INT cell_to_collapse[MAX_CELL_COLLAPSE];
  REF_INT id0, id1;

  REF_BOOL geom_node1;
  REF_BOOL geom_edge1;
  REF_BOOL edge_side;

  *allowed = REF_FALSE;

  /* don't remove a geometry CAD node */
  RSS(ref_geom_is_a(ref_geom, node1, REF_GEOM_NODE, &geom_node1), "node check");
  if (geom_node1) {
    *allowed = REF_FALSE;
    return REF_SUCCESS;
  }

  /* geometery edge allowed if collapse is on edge */
  RSS(ref_geom_is_a(ref_geom, node1, REF_GEOM_EDGE, &geom_edge1),
      "edge check 1");
  if (geom_edge1) {
    RSS(ref_cell_has_side(ref_edg, node0, node1, &edge_side),
        "allowed if a side of a triangle");
    if (!edge_side) return REF_SUCCESS;
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }

  /* ids1 is a list of degree1 face ids for node1 */
  RXS(ref_cell_id_list_around(ref_tri, node1, 3, &degree1, ids1),
      REF_INCREASE_LIMIT, "count faceids");

  switch (degree1) {
    case 3: /* geometry node never allowed to move */
      *allowed = REF_FALSE;
      break;
    case 2: /* geometery edge allowed if collapse is on edge */
      RSS(ref_cell_list_with2(ref_tri, node0, node1, MAX_CELL_COLLAPSE, &ncell,
                              cell_to_collapse),
          "list");
      if (2 != ncell) {
        *allowed = REF_FALSE;
        break;
      }
      RSS(ref_cell_nodes(ref_tri, cell_to_collapse[0], nodes), "nodes");
      id0 = nodes[3];
      RSS(ref_cell_nodes(ref_tri, cell_to_collapse[1], nodes), "nodes");
      id1 = nodes[3];
      if ((id0 == ids1[0] && id1 == ids1[1]) ||
          (id1 == ids1[0] && id0 == ids1[1]))
        *allowed = REF_TRUE;
      break;
    case 1: /* geometry face allowed if on that face */
      RSS(ref_cell_has_side(ref_tri, node0, node1, allowed),
          "allowed if a side of a triangle");
      break;
    case 0: /* volume node always allowed */
      *allowed = REF_TRUE;
      break;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_manifold(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_BOOL *allowed) {
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, new_cell;
  REF_BOOL will_be_collapsed;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);

  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) will_be_collapsed = REF_TRUE;
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = node0;

    RXS(ref_cell_with(ref_cell, nodes, &new_cell), REF_NOT_FOUND,
        "with node0 failed");
    if (REF_EMPTY != new_cell) {
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  ref_cell = ref_grid_edg(ref_grid);

  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) will_be_collapsed = REF_TRUE;
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = node0;

    RXS(ref_cell_with(ref_cell, nodes, &new_cell), REF_NOT_FOUND,
        "with node0 failed");
    if (REF_EMPTY != new_cell) {
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_chord_height(REF_GRID ref_grid, REF_INT node0,
                                          REF_INT node1, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_BOOL will_be_collapsed;
  REF_INT i;
  REF_DBL n0[3], n1[3], cross_prod[3];
  REF_DBL new_length, cross_length, chord, chord_ratio;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_edg(ref_grid);

  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) will_be_collapsed = REF_TRUE;
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node1 == nodes[node]) {
        /* triangle node0, node1, nodes[node] */
        /* nodes[node] = node0; */
        for (i = 0; i < 3; i++) {
          n0[i] = ref_node_xyz(ref_node, i, node1) -
                  ref_node_xyz(ref_node, i, node0);
          n1[i] = ref_node_xyz(ref_node, i, nodes[1 - node]) -
                  ref_node_xyz(ref_node, i, node0);
        }
        new_length = sqrt(ref_math_dot(n1, n1));
        ref_math_cross_product(n0, n1, cross_prod);
        cross_length = sqrt(ref_math_dot(cross_prod, cross_prod));
        /* |n0 x n1| = |n0||n1|sin(t) */
        if (ref_math_divisible(cross_length, new_length)) {
          chord = cross_length / new_length;
          if (ref_math_divisible(chord, new_length)) {
            chord_ratio = chord / new_length;
            if (chord_ratio > 0.1) return REF_SUCCESS;
          } else {
            return REF_SUCCESS;
          }
        } else {
          return REF_SUCCESS;
        }
      }
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_same_normal(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL n0[3], n1[3];
  REF_DBL dot;
  REF_STATUS status;

  *allowed = REF_TRUE;

  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    /* a triangle with node0 and node1 will be removed */
    if (node0 == ref_cell_c2n(ref_cell, 0, cell) ||
        node0 == ref_cell_c2n(ref_cell, 1, cell) ||
        node0 == ref_cell_c2n(ref_cell, 2, cell))
      continue;
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    RSS(ref_node_tri_normal(ref_node, nodes, n0), "orig normal");
    RSS(ref_math_normalize(n0), "original triangle has zero area");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node1 == nodes[node]) {
        nodes[node] = node0;
      }
    }
    RSS(ref_node_tri_normal(ref_node, nodes, n1), "new normal");
    status = ref_math_normalize(n1);
    if (REF_DIV_ZERO == status) { /* new triangle face has zero area */
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
    RSS(status, "new normal length")
    dot = ref_math_dot(n0, n1);
    if (dot < ref_node_same_normal_tol(ref_node)) {
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_mixed(REF_GRID ref_grid, REF_INT node0,
                                   REF_INT node1, REF_BOOL *allowed) {
  SUPRESS_UNUSED_COMPILER_WARNING(node0);

  *allowed = (ref_adj_empty(ref_cell_adj(ref_grid_pyr(ref_grid)), node1) &&
              ref_adj_empty(ref_cell_adj(ref_grid_pri(ref_grid)), node1) &&
              ref_adj_empty(ref_cell_adj(ref_grid_hex(ref_grid)), node1));

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_local_cell(REF_GRID ref_grid, REF_INT node0,
                                        REF_INT node1, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, node;

  *allowed = REF_FALSE;

  /* may be able to relax node0 local if geom constraint is o.k. */

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        return REF_SUCCESS;
      }
    }
  }
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        return REF_SUCCESS;
      }
    }
  }

  /* for parallel surf if ever needed */
  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        return REF_SUCCESS;
      }
    }
  }
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        return REF_SUCCESS;
      }
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

/* triangle can not be fully constrained by CAD edges after a collapse */
REF_STATUS ref_collapse_edge_cad_constrained(REF_GRID ref_grid, REF_INT node0,
                                             REF_INT node1, REF_BOOL *allowed) {
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT item, cell, node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL will_be_collapsed;
  REF_BOOL edge0, edge1, edge2;

  *allowed = REF_TRUE;

  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) will_be_collapsed = REF_TRUE;
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = node0;
    RSS(ref_geom_is_a(ref_geom, nodes[0], REF_GEOM_EDGE, &edge0), "e0");
    RSS(ref_geom_is_a(ref_geom, nodes[1], REF_GEOM_EDGE, &edge1), "e1");
    RSS(ref_geom_is_a(ref_geom, nodes[2], REF_GEOM_EDGE, &edge2), "e2");
    if (edge0 && edge1 && edge2) {
      *allowed = REF_FALSE;
      break;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_tri_quality(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_DBL quality;
  REF_BOOL will_be_collapsed;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node0 == nodes[node]) {
        will_be_collapsed = REF_TRUE;
      }
    }
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = node0;
    RSS(ref_node_tri_quality(ref_node, nodes, &quality), "qual");
    if (quality < ref_grid_adapt(ref_grid, collapse_quality_absolute))
      return REF_SUCCESS;
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_tet_quality(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, ntri;
  REF_DBL quality;
  REF_BOOL will_be_collapsed;

  *allowed = REF_FALSE;

  if (ref_grid_surf(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }

  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) will_be_collapsed = REF_TRUE;
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = node0;
    if (ref_grid_surf(ref_grid)) {
      RSS(ref_node_tri_quality(ref_node, nodes, &quality), "qual");
    } else {
      RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    }
    if (quality < ref_grid_adapt(ref_grid, collapse_quality_absolute))
      return REF_SUCCESS;
    if (!ref_grid_surf(ref_grid)) {
      RSS(ref_cell_ntri_with_tet_nodes(ref_grid_tri(ref_grid), nodes, &ntri),
          "count boundary triangles");
      if (ntri > 1) return REF_SUCCESS;
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_ratio(REF_GRID ref_grid, REF_INT node0,
                                   REF_INT node1, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_BOOL will_be_collapsed;
  REF_DBL edge_ratio;
  REF_DBL old_max, old_min, new_max, new_min;

  *allowed = REF_FALSE;

  if (ref_grid_surf(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }

  old_max = -1.0;
  old_min = REF_DBL_MAX;
  new_max = -1.0;
  new_min = REF_DBL_MAX;
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node1 != nodes[node]) {
        RSS(ref_node_ratio(ref_node, node1, nodes[node], &edge_ratio), "ratio");
        old_max = MAX(old_max, edge_ratio);
        old_min = MIN(old_min, edge_ratio);
      }
    }

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) will_be_collapsed = REF_TRUE;
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 != nodes[node]) {
        RSS(ref_node_ratio(ref_node, node0, nodes[node], &edge_ratio), "ratio");
        new_max = MAX(new_max, edge_ratio);
        new_min = MIN(new_min, edge_ratio);
      }
  }

  if ((new_min >= ref_grid_adapt(ref_grid, post_min_ratio)) &&
      (new_max <= ref_grid_adapt(ref_grid, post_max_ratio))) {
    *allowed = REF_TRUE;
  }

  if (new_min >= old_min && new_max <= old_max) {
    *allowed = REF_TRUE;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_surf_ratio(REF_GRID ref_grid, REF_INT node0,
                                   REF_INT node1, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_BOOL has_side;
  REF_INT node, nnode, node_list[MAX_NODE_LIST];
  REF_DBL ratio;

  *allowed = REF_TRUE;
  RSS(ref_cell_has_side(ref_cell, node0, node1, &has_side), "side");
  if (has_side) {
    RSS(ref_cell_node_list_around(ref_cell, node1, MAX_NODE_LIST, &nnode,
                                  node_list),
        "da hood");
    for (node = 0; node < nnode; node++) {
      RSS(ref_node_ratio(ref_node, node_list[node], node1, &ratio), "ratio");
      if (ratio < ref_grid_adapt(ref_grid, collapse_ratio)) return REF_SUCCESS;
    }
    *allowed = REF_TRUE;
  }
  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_normdev(REF_GRID ref_grid, REF_INT node0,
                                     REF_INT node1, REF_BOOL *allowed) {
  REF_CELL ref_cell;
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER], new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_BOOL will_be_collapsed;

  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_BOOL node0_support, node1_support, supported;
  REF_DBL orig_dev, new_dev;
  REF_DBL sign_uv_area, orig_uv_area, new_uv_area;

  RSS(ref_geom_supported(ref_geom, node0, &node0_support), "support0");
  RSS(ref_geom_supported(ref_geom, node1, &node1_support), "support1");
  if (!ref_geom_model_loaded(ref_geom) || !node0_support || !node1_support) {
    *allowed = REF_TRUE;
    return REF_SUCCESS;
  }

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node0 == nodes[node]) {
        will_be_collapsed = REF_TRUE;
      }
    }
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      new_nodes[node] = nodes[node];
      if (node1 == new_nodes[node]) new_nodes[node] = node0;
    }
    new_nodes[ref_cell_node_per(ref_cell)] = nodes[ref_cell_node_per(ref_cell)];

    /* see if new config is below limit */
    RSS(ref_geom_cell_tuv_supported(ref_geom, nodes, REF_GEOM_FACE, &supported),
        "tuv support fororiginal configuration");
    RAS(supported,
        "original configuration before collapse does not support cell tuv");
    RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &orig_dev), "orig");
    RSS(ref_geom_cell_tuv_supported(ref_geom, new_nodes, REF_GEOM_FACE,
                                    &supported),
        "tuv support for swapped configuration");
    if (!supported) { /* abort collapse, cell tuv not supported */
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
    RSS(ref_geom_tri_norm_deviation(ref_grid, new_nodes, &new_dev), "new");
    RSS(ref_geom_uv_area_sign(ref_grid, nodes[ref_cell_node_per(ref_cell)],
                              &sign_uv_area),
        "sign");
    RSS(ref_geom_uv_area(ref_geom, nodes, &orig_uv_area), "uv area");
    RSS(ref_geom_uv_area(ref_geom, new_nodes, &new_uv_area), "uv area");
    /* allow if improvement */
    if (((new_dev < ref_grid_adapt(ref_grid, post_min_normdev)) &&
         (new_dev < orig_dev)) ||
        ((sign_uv_area * new_uv_area < ref_node_min_uv_area(ref_node)) &&
         (sign_uv_area * new_uv_area < sign_uv_area * orig_uv_area))) {
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_face(REF_GRID ref_grid, REF_INT keep, REF_INT remove) {
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT ncell, cell_in_list;
  REF_INT cell_to_collapse[MAX_CELL_COLLAPSE];

  ref_cell = ref_grid_edg(ref_grid);

  RSS(ref_cell_list_with2(ref_cell, keep, remove, MAX_CELL_COLLAPSE, &ncell,
                          cell_to_collapse),
      "lst");

  for (cell_in_list = 0; cell_in_list < ncell; cell_in_list++) {
    cell = cell_to_collapse[cell_in_list];
    RSS(ref_cell_remove(ref_cell, cell), "remove");
  }
  RSS(ref_cell_replace_node(ref_cell, remove, keep), "replace node");

  ref_cell = ref_grid_tri(ref_grid);

  RSS(ref_cell_list_with2(ref_cell, keep, remove, MAX_CELL_COLLAPSE, &ncell,
                          cell_to_collapse),
      "lst");
  for (cell_in_list = 0; cell_in_list < ncell; cell_in_list++) {
    cell = cell_to_collapse[cell_in_list];
    RSS(ref_cell_remove(ref_cell, cell), "remove");
  }
  RSS(ref_cell_replace_node(ref_cell, remove, keep), "replace node");

  RSS(ref_node_remove(ref_grid_node(ref_grid), remove), "rm");
  RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), remove), "rm");

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_face_local_pris(REF_GRID ref_grid, REF_INT keep,
                                        REF_INT remove, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, node;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_pri(ref_grid);

  each_ref_cell_having_node(ref_cell, remove, item, cell) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        return REF_SUCCESS;
      }
    }
  }

  /* may be able to relax keep local if geom constraint is o.k. */
  each_ref_cell_having_node(ref_cell, keep, item, cell) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (!ref_node_owned(ref_node, ref_cell_c2n(ref_cell, node, cell))) {
        return REF_SUCCESS;
      }
    }
  }
  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_face_quality(REF_GRID ref_grid, REF_INT keep,
                                     REF_INT remove, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_DBL quality;
  REF_BOOL will_be_collapsed;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, remove, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (keep == nodes[node]) {
        will_be_collapsed = REF_TRUE;
      }
    }
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (remove == nodes[node]) nodes[node] = keep;
    RSS(ref_node_tri_quality(ref_node, nodes, &quality), "qual");
    if (quality < ref_grid_adapt(ref_grid, collapse_quality_absolute))
      return REF_SUCCESS;
  }

  /* FIXME check quads too */

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_face_ratio(REF_GRID ref_grid, REF_INT keep,
                                   REF_INT remove, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_BOOL will_be_collapsed;
  REF_DBL ratio;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, remove, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (keep == nodes[node]) {
        will_be_collapsed = REF_TRUE;
      }
    }
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (remove != nodes[node]) {
        RSS(ref_node_ratio(ref_node, keep, nodes[node], &ratio), "ratio");
        if (ratio > ref_grid_adapt(ref_grid, post_max_ratio))
          return REF_SUCCESS;
      }
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_face_outward_norm(REF_GRID ref_grid, REF_INT keep,
                                          REF_INT remove, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_BOOL will_be_collapsed, valid;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, remove, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");

    will_be_collapsed = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (keep == nodes[node]) {
        will_be_collapsed = REF_TRUE;
      }
    }
    if (will_be_collapsed) continue;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (remove == nodes[node]) {
        nodes[node] = keep;
      }
    }

    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    if (!valid) return REF_SUCCESS;
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_face_geometry(REF_GRID ref_grid, REF_INT keep,
                                      REF_INT remove, REF_BOOL *allowed) {
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT deg, degree1;
  REF_INT id, ids1[2];
  REF_BOOL already_have_it, mixed;

  degree1 = 0;
  ids1[0] = REF_EMPTY;
  ids1[1] = REF_EMPTY;
  each_ref_cell_having_node(ref_cell, remove, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    id = nodes[ref_cell_node_per(ref_cell)];
    already_have_it = REF_FALSE;
    for (deg = 0; deg < degree1; deg++) {
      if (id == ids1[deg]) {
        already_have_it = REF_TRUE;
      }
    }
    if (!already_have_it) {
      ids1[degree1] = id;
      degree1++;
      if (2 == degree1) break;
    }
  }

  *allowed = REF_FALSE;

  switch (degree1) {
    case 2: /* geometry node never allowed to move */
      *allowed = REF_FALSE;
      break;
    case 1:
      RSS(ref_cell_has_side(ref_grid_qua(ref_grid), keep, remove, &mixed),
          "not allowed if a side of a hex, mixed");
      if (mixed) {
        *allowed = REF_FALSE;
      } else {
        /* geometry face allowed if on that face */
        RSS(ref_cell_has_side(ref_cell, keep, remove, allowed),
            "allowed if a side of a quad");
      }
      break;
    case 0: /* volume node always allowed */
      *allowed = REF_TRUE;
      break;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_face_same_tangent(REF_GRID ref_grid, REF_INT keep,
                                          REF_INT remove, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, cell;
  REF_INT other, node, ixyz;
  REF_DBL tangent0[3], tangent1[3];
  REF_DBL dot;
  REF_STATUS status;

  *allowed = REF_TRUE;

  each_ref_cell_having_node(ref_cell, remove, item, cell) {
    /* a quad with keep and remove will be removed */
    if (keep == ref_cell_c2n(ref_cell, 0, cell) ||
        keep == ref_cell_c2n(ref_cell, 1, cell))
      continue;
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    other = REF_EMPTY;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (nodes[node] != remove) {
        other = nodes[node];
      }
    }
    RAS(REF_EMPTY != other, "other not found");
    for (ixyz = 0; ixyz < 3; ixyz++) {
      tangent0[ixyz] = ref_node_xyz(ref_node, ixyz, remove) -
                       ref_node_xyz(ref_node, ixyz, other);
      tangent1[ixyz] = ref_node_xyz(ref_node, ixyz, keep) -
                       ref_node_xyz(ref_node, ixyz, other);
    }
    RSS(ref_math_normalize(tangent0), "zero length orig quad");
    status = ref_math_normalize(tangent1);
    if (REF_DIV_ZERO == status) { /* new quad face has zero area */
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
    dot = ref_math_dot(tangent0, tangent1);
    if (dot < ref_node_same_normal_tol(ref_node)) {
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_twod_pass(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *order;
  REF_INT ntarget, *target, *node2target;
  REF_INT node, node0, node1;
  REF_INT i, edge;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL edge_ratio;

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  ref_malloc_init(ratio, ref_node_max(ref_node), REF_DBL,
                  2.0 * ref_grid_adapt(ref_grid, collapse_ratio));

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);

    RSS(ref_node_ratio(ref_node, node0, node1, &edge_ratio), "ratio");
    ratio[node0] = MIN(ratio[node0], edge_ratio);
    ratio[node1] = MIN(ratio[node1], edge_ratio);
  }

  ref_malloc(target, ref_node_n(ref_node), REF_INT);
  ref_malloc_init(node2target, ref_node_max(ref_node), REF_INT, REF_EMPTY);

  ntarget = 0;
  for (node = 0; node < ref_node_max(ref_node); node++)
    if (ratio[node] < ref_grid_adapt(ref_grid, collapse_ratio)) {
      node2target[node] = ntarget;
      target[ntarget] = node;
      ratio[ntarget] = ratio[node];
      ntarget++;
    }

  ref_malloc(order, ntarget, REF_INT);

  RSS(ref_sort_heap_dbl(ntarget, ratio, order), "sort lengths");

  for (i = 0; i < ntarget; i++) {
    if (ratio[order[i]] > ref_grid_adapt(ref_grid, collapse_ratio)) continue;
    node1 = target[order[i]];
    RSS(ref_collapse_face_remove_node1(ref_grid, &node0, node1), "collapse rm");
    if (!ref_node_valid(ref_node, node1)) {
      ref_node_age(ref_node, node0) = 0;
      each_ref_cell_having_node(ref_cell, node0, item, cell) {
        RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
        for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
          if (REF_EMPTY != node2target[nodes[node]]) {
            ratio[node2target[nodes[node]]] =
                2.0 * ref_grid_adapt(ref_grid, collapse_ratio);
          }
        }
      }
    }
  }

  ref_free(order);
  ref_free(node2target);
  ref_free(target);
  ref_free(ratio);

  ref_edge_free(ref_edge);

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_face_remove_node1(REF_GRID ref_grid,
                                          REF_INT *actual_node0,
                                          REF_INT node1) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT nnode, node;
  REF_INT node_to_collapse[MAX_NODE_LIST];
  REF_INT order[MAX_NODE_LIST];
  REF_DBL ratio_to_collapse[MAX_NODE_LIST];
  REF_INT node0;
  REF_BOOL allowed;
  REF_BOOL verbose = REF_FALSE;

  *actual_node0 = REF_EMPTY;

  RSS(ref_cell_node_list_around(ref_cell, node1, MAX_NODE_LIST, &nnode,
                                node_to_collapse),
      "da hood");
  for (node = 0; node < nnode; node++) {
    RSS(ref_node_ratio(ref_node, node_to_collapse[node], node1,
                       &(ratio_to_collapse[node])),
        "ratio");
  }

  RSS(ref_sort_heap_dbl(nnode, ratio_to_collapse, order), "sort lengths");

  for (node = 0; node < nnode; node++) {
    node0 = node_to_collapse[order[node]];

    RSS(ref_collapse_face_geometry(ref_grid, node0, node1, &allowed),
        "col geom");
    if (!allowed && verbose) printf("%d geom\n", node);
    if (!allowed) continue;

    RSS(ref_collapse_face_same_tangent(ref_grid, node0, node1, &allowed),
        "tan");
    if (!allowed && verbose) printf("%d tang\n", node);
    if (!allowed) continue;

    RSS(ref_collapse_face_outward_norm(ref_grid, node0, node1, &allowed),
        "norm");
    if (!allowed && verbose) printf("%d outw\n", node);
    if (!allowed) continue;

    RSS(ref_collapse_face_ratio(ref_grid, node0, node1, &allowed), "qual");
    if (!allowed && verbose) printf("%d ratio\n", node);
    if (!allowed) continue;

    RSS(ref_collapse_face_quality(ref_grid, node0, node1, &allowed), "qual");
    if (!allowed && verbose) printf("%d qual\n", node);
    if (!allowed) continue;

    RSS(ref_collapse_face_local_pris(ref_grid, node0, node1, &allowed),
        "colloc");
    if (!allowed && verbose) printf("%d loca\n", node);
    if (!allowed) {
      ref_node_age(ref_node, node0)++;
      ref_node_age(ref_node, node1)++;
      continue;
    }

    if (verbose) printf("%d split!\n", node);

    *actual_node0 = node0;

    RSS(ref_collapse_face(ref_grid, node0, node1), "col!");

    break;
  }

  return REF_SUCCESS;
}
