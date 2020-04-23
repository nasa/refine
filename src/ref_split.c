
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

#include "ref_split.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_adapt.h"
#include "ref_cavity.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_subdiv.h"

#define MAX_CELL_SPLIT (100)

REF_STATUS ref_split_pass(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *edges, *order;
  REF_INT i, n, edge;
  REF_BOOL allowed_ratio, allowed_tri_conformity, allowed_tet_quality;
  REF_BOOL allowed, allowed_local, geom_support, valid_cavity, try_cavity;
  REF_BOOL allowed_cavity_ratio, has_edge;
  REF_BOOL allowed_tri_quality;
  REF_DBL min_del, min_add;
  REF_GLOB global;
  REF_INT new_node;
  REF_CAVITY ref_cavity = (REF_CAVITY)NULL;
  REF_BOOL span_parts;
  REF_LIST para_no_geom = NULL;
  REF_LIST para_cavity = NULL;
  REF_SUBDIV ref_subdiv = NULL;
  REF_STATUS status;
  REF_BOOL transcript = REF_FALSE;
  REF_DBL ratio01, ratio0, ratio1, weight_node1;

  ref_cell = ref_grid_tet(ref_grid);
  if (ref_grid_twod(ref_grid) || ref_grid_surf(ref_grid))
    ref_cell = ref_grid_tri(ref_grid);

  span_parts = ref_mpi_para(ref_grid_mpi(ref_grid)) &&
               !ref_grid_twod(ref_grid) && !ref_grid_surf(ref_grid);

  if (span_parts) {
    RSS(ref_list_create(&para_no_geom), "list for stuck edges");
    RSS(ref_list_create(&para_cavity), "list for stuck cavity");
  }

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  ref_malloc(ratio, ref_edge_n(ref_edge), REF_DBL);
  ref_malloc(order, ref_edge_n(ref_edge), REF_INT);
  ref_malloc(edges, ref_edge_n(ref_edge), REF_INT);

  n = 0;
  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_node_ratio(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                       ref_edge_e2n(ref_edge, 1, edge), &(ratio[n])),
        "ratio");
    if (ratio[n] > ref_grid_adapt(ref_grid, split_ratio)) {
      edges[n] = edge;
      n++;
    }
  }

  RSS(ref_sort_heap_dbl(n, ratio, order), "sort lengths");

  for (i = n - 1; i >= 0; i--) {
    edge = edges[order[i]];

    /* transcript = (ratio[order[i]] > 3.0); */
    if (transcript) printf("transcript on ratio %f\n", ratio[order[i]]);

    RSS(ref_cell_has_side(ref_cell, ref_edge_e2n(ref_edge, 0, edge),
                          ref_edge_e2n(ref_edge, 1, edge), &allowed),
        "has side");
    if (transcript && !allowed) printf("not a side anymore\n");
    if (!allowed) continue;

    /* skip if neither node is owned */
    if (!ref_node_owned(ref_node, ref_edge_e2n(ref_edge, 0, edge)) &&
        !ref_node_owned(ref_node, ref_edge_e2n(ref_edge, 1, edge))) {
      if (transcript) printf("neither node is local\n");
      continue;
    }

    RSS(ref_split_edge_mixed(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                             ref_edge_e2n(ref_edge, 1, edge), &allowed),
        "mixed");
    if (transcript && !allowed) printf("mixed edge\n");
    if (!allowed) continue;

    weight_node1 = 0.5;
    if (ref_grid_twod(ref_grid) || ref_grid_surf(ref_grid)) {
      RSS(ref_node_ratio(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                         ref_edge_e2n(ref_edge, 1, edge), &ratio01),
          "ratio01");
      RSS(ref_node_ratio_node0(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                               ref_edge_e2n(ref_edge, 1, edge), &ratio0),
          "ratio0");
      RSS(ref_node_ratio_node0(ref_node, ref_edge_e2n(ref_edge, 1, edge),
                               ref_edge_e2n(ref_edge, 0, edge), &ratio1),
          "ratio1");
      if (ref_math_divisible(ratio0, ratio1 + ratio0)) {
        if (0.25 < ratio0 / (ratio0 + ratio1) &&
            ratio0 / (ratio0 + ratio1) < 0.75) {
          weight_node1 = 1.0 - ratio0 / (ratio0 + ratio1);
        } else {
          if (ratio0 < ratio1) {
            if (ref_math_divisible(ratio0, ratio01))
              weight_node1 = 1.0 - ratio0 / ratio01;
          } else {
            if (ref_math_divisible(ratio1, ratio01))
              weight_node1 = ratio1 / ratio01;
          }
        }
      }
    }
    weight_node1 = MIN(1, MAX(0, weight_node1));

    RSS(ref_node_next_global(ref_node, &global), "next global");
    RSS(ref_node_add(ref_node, global, &new_node), "new node");
    RSS(ref_node_interpolate_edge(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                                  ref_edge_e2n(ref_edge, 1, edge), weight_node1,
                                  new_node),
        "interp new node");
    RSS(ref_geom_add_between(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                             ref_edge_e2n(ref_edge, 1, edge), weight_node1,
                             new_node),
        "geom new node");
    RSS(ref_geom_constrain(ref_grid, new_node), "geom constraint");
    RSS(ref_metric_interpolate_between(
            ref_grid, ref_edge_e2n(ref_edge, 0, edge),
            ref_edge_e2n(ref_edge, 1, edge), new_node),
        "interp new node metric");
    RSS(ref_geom_supported(ref_grid_geom(ref_grid), new_node, &geom_support),
        "geom support");
    if (transcript && geom_support) printf("geom support\n");

    if (transcript)
      printf("weight_node1 %f xyz %f %f %f\n", weight_node1,
             ref_node_xyz(ref_node, 0, new_node),
             ref_node_xyz(ref_node, 1, new_node),
             ref_node_xyz(ref_node, 2, new_node));

    RSS(ref_split_edge_tet_quality(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                                   ref_edge_e2n(ref_edge, 1, edge), new_node,
                                   &allowed_tet_quality),
        "edge tet qual");
    if (transcript && !allowed_tet_quality) printf("tet quality poor\n");

    RSS(ref_split_edge_tri_quality(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                                   ref_edge_e2n(ref_edge, 1, edge), new_node,
                                   &allowed_tri_quality),
        "quality of new tri");
    if (transcript && !allowed_tri_quality) printf("tri quality poor\n");

    RSS(ref_split_edge_ratio(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                             ref_edge_e2n(ref_edge, 1, edge), new_node,
                             &allowed_ratio),
        "edge tet ratio");
    if (transcript && !allowed_ratio) printf("ratio poor\n");

    RSS(ref_split_edge_tri_conformity(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                                      ref_edge_e2n(ref_edge, 1, edge), new_node,
                                      &allowed_tri_conformity),
        "edge tri qual");
    if (transcript && !allowed_tri_conformity) printf("tri conformity poor\n");

    RSS(ref_cell_has_side(ref_grid_edg(ref_grid),
                          ref_edge_e2n(ref_edge, 0, edge),
                          ref_edge_e2n(ref_edge, 1, edge), &has_edge),
        "check for an edge");
    if (transcript && has_edge) printf("has geom edge\n");

    try_cavity = REF_FALSE;
    if (!allowed_tet_quality || !allowed_ratio || !allowed_tri_conformity ||
        !allowed_tri_quality) {
      if (geom_support) {
        try_cavity = REF_TRUE;
      } else {
        RSS(ref_node_remove(ref_node, new_node), "remove new node");
        RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), new_node), "rm");
        continue;
      }
    }

    if (try_cavity) {
      RSS(ref_cavity_create(&ref_cavity), "cav create");
      RSS(ref_cavity_form_edge_split(ref_cavity, ref_grid,
                                     ref_edge_e2n(ref_edge, 0, edge),
                                     ref_edge_e2n(ref_edge, 1, edge), new_node),
          "form edge split cav");
      if (transcript)
        printf("try cavity status %d\n", (int)ref_cavity_state(ref_cavity));
      if (REF_SUCCESS != ref_cavity_enlarge_combined(ref_cavity)) {
        RSS(ref_node_location(ref_node, ref_edge_e2n(ref_edge, 0, edge)), "n0");
        RSS(ref_node_location(ref_node, ref_edge_e2n(ref_edge, 1, edge)), "n1");
        REF_WHERE("enlarge"); /* note but skip cavity failures */
      }
      if (REF_CAVITY_VISIBLE == ref_cavity_state(ref_cavity)) {
        if (transcript) printf("cavity visible\n");
        RSS(ref_cavity_ratio(ref_cavity, &allowed_cavity_ratio),
            "cavity ratio");
        RSS(ref_cavity_change(ref_cavity, &min_del, &min_add), "cavity change");
        valid_cavity =
            (allowed_cavity_ratio || has_edge) &&
            (min_add > ref_grid_adapt(ref_grid, split_quality_absolute));
        if (transcript && !valid_cavity)
          printf("valid_cavity edge %d ratio %d add %d %f\n", has_edge,
                 allowed_cavity_ratio,
                 (min_add > ref_grid_adapt(ref_grid, split_quality_absolute)),
                 min_add);

        if (valid_cavity) {
          if (transcript) printf("cavity replace\n");
          RSS(ref_cavity_replace(ref_cavity), "cav replace");
          RSS(ref_cavity_free(ref_cavity), "cav free");
          ref_cavity = (REF_CAVITY)NULL;
          ref_node_age(ref_node, ref_edge_e2n(ref_edge, 0, edge)) = 0;
          ref_node_age(ref_node, ref_edge_e2n(ref_edge, 1, edge)) = 0;
          RSS(ref_smooth_post_edge_split(ref_grid, new_node),
              "smooth after split");
          continue;
        }
      } else {
        if (transcript)
          printf("cavity not visible %d\n", (int)ref_cavity_state(ref_cavity));
      }
      if (REF_CAVITY_PARTITION_CONSTRAINED == ref_cavity_state(ref_cavity)) {
        if (span_parts) RSS(ref_list_push(para_cavity, edge), "push");
      }
      RSS(ref_cavity_free(ref_cavity), "cav free");
      ref_cavity = (REF_CAVITY)NULL;
      RSS(ref_node_remove(ref_node, new_node), "remove new node");
      RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), new_node), "rm");
      continue;
    }

    RSS(ref_cell_local_gem(ref_cell, ref_node, ref_edge_e2n(ref_edge, 0, edge),
                           ref_edge_e2n(ref_edge, 1, edge), &allowed_local),
        "local tet");
    if (!allowed_local) {
      if (span_parts) {
        RSS(ref_list_push(para_no_geom, edge), "push");
      } else {
        ref_node_age(ref_node, ref_edge_e2n(ref_edge, 0, edge))++;
        ref_node_age(ref_node, ref_edge_e2n(ref_edge, 1, edge))++;
      }
      RSS(ref_node_remove(ref_node, new_node), "remove new node");
      RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), new_node), "rm");
      continue;
    }

    if (transcript) printf("split\n");
    status = ref_split_edge(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                            ref_edge_e2n(ref_edge, 1, edge), new_node);
    if (REF_INCREASE_LIMIT == status) {
      RSS(ref_node_remove(ref_node, new_node), "remove new node");
      RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), new_node), "rm");
      continue;
    }
    RSS(status, "tet edge split");

    ref_node_age(ref_node, ref_edge_e2n(ref_edge, 0, edge)) = 0;
    ref_node_age(ref_node, ref_edge_e2n(ref_edge, 1, edge)) = 0;

    RSS(ref_smooth_post_edge_split(ref_grid, new_node), "smooth after split");
  }

  ref_free(edges);
  ref_free(order);
  ref_free(ratio);

  if (span_parts) {
    if (ref_grid_adapt(ref_grid, instrument))
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "spl int");
    RSS(ref_subdiv_create(&ref_subdiv, ref_grid), "create");
    ref_subdiv_new_mark_allowed(ref_subdiv) = REF_FALSE;
    ref_subdiv->instrument = REF_TRUE;
    each_ref_list_item(para_no_geom, i) {
      edge = ref_list_value(para_no_geom, i);
      RSS(ref_cell_has_side(ref_grid_tet(ref_grid),
                            ref_edge_e2n(ref_edge, 0, edge),
                            ref_edge_e2n(ref_edge, 1, edge), &allowed),
          "has side");
      if (!allowed) continue;

      RSS(ref_subdiv_mark_to_split(ref_subdiv, ref_edge_e2n(ref_edge, 0, edge),
                                   ref_edge_e2n(ref_edge, 1, edge)),
          "mark edge to para split");
    }

    RSS(ref_subdiv_split(ref_subdiv), "split");
    RSS(ref_subdiv_free(ref_subdiv), "free");

    {
      REF_INT n_cavity;
      n_cavity = ref_list_n(para_cavity);
      RSS(ref_mpi_allsum(ref_mpi, &n_cavity, 1, REF_INT_TYPE), "cavs");
      if (ref_mpi_once(ref_mpi))
        printf("cavities with part constraints %d\n", n_cavity);
    }

    ref_list_free(para_no_geom);
    ref_list_free(para_cavity);
  }

  ref_edge_free(ref_edge);

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                          REF_INT new_node) {
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ncell, cell_in_list;
  REF_INT cell_to_split[MAX_CELL_SPLIT];
  REF_INT node, new_cell;
  REF_STATUS status;

  ref_cell = ref_grid_tet(ref_grid);
  status = ref_cell_list_with2(ref_cell, node0, node1, MAX_CELL_SPLIT, &ncell,
                               cell_to_split);
  if (REF_INCREASE_LIMIT == status) return status;
  RSS(status, "tet list to split");

  for (cell_in_list = 0; cell_in_list < ncell; cell_in_list++) {
    cell = cell_to_split[cell_in_list];
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    RSS(ref_cell_remove(ref_cell, cell), "remove");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node0 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node0;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node1 version");
  }

  ref_cell = ref_grid_tri(ref_grid);
  RSS(ref_cell_list_with2(ref_cell, node0, node1, MAX_CELL_SPLIT, &ncell,
                          cell_to_split),
      "get tri list, should have been smaller then tets");

  for (cell_in_list = 0; cell_in_list < ncell; cell_in_list++) {
    cell = cell_to_split[cell_in_list];
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    RSS(ref_cell_remove(ref_cell, cell), "remove");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node0 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node0;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node1 version");
  }

  ref_cell = ref_grid_edg(ref_grid);
  RSS(ref_cell_list_with2(ref_cell, node0, node1, MAX_CELL_SPLIT, &ncell,
                          cell_to_split),
      "get edg list, should have been smaller then tets");

  for (cell_in_list = 0; cell_in_list < ncell; cell_in_list++) {
    cell = cell_to_split[cell_in_list];
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    RSS(ref_cell_remove(ref_cell, cell), "remove");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node0 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node0;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node1 version");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_split_face(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                          REF_INT node2, REF_INT new_node) {
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER], face_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell0, cell1;
  REF_INT node, new_cell;

  face_nodes[0] = node0;
  face_nodes[1] = node1;
  face_nodes[2] = node2;
  face_nodes[3] = node0;

  ref_cell = ref_grid_tet(ref_grid);
  RSS(ref_cell_with_face(ref_cell, face_nodes, &cell0, &cell1), "get tet(2)");
  if (REF_EMPTY != cell0) {
    cell = cell0;
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    RSS(ref_cell_remove(ref_cell, cell), "remove");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node0 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node0;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node1 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node1;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node2 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node2 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node2;
  }
  if (REF_EMPTY != cell1) {
    cell = cell1;
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    RSS(ref_cell_remove(ref_cell, cell), "remove");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node0 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node0;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node1 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node1;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node2 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node2 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node2;
  }

  ref_cell = ref_grid_tri(ref_grid);
  RXS(ref_cell_with(ref_cell, face_nodes, &cell), REF_NOT_FOUND, "find tri");
  if (REF_EMPTY != cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    RSS(ref_cell_remove(ref_cell, cell), "remove");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node0 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node0;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node1 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node1;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node2 == nodes[node]) nodes[node] = new_node;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add node2 version");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node2;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_mixed(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_BOOL *allowed) {
  REF_BOOL pyr_side, pri_side, hex_side, qua_side;

  RSS(ref_cell_has_side(ref_grid_pyr(ref_grid), node0, node1, &pyr_side),
      "pyr");
  RSS(ref_cell_has_side(ref_grid_pri(ref_grid), node0, node1, &pri_side),
      "pri");
  RSS(ref_cell_has_side(ref_grid_hex(ref_grid), node0, node1, &hex_side),
      "hex");
  RSS(ref_cell_has_side(ref_grid_qua(ref_grid), node0, node1, &qua_side),
      "qua");

  *allowed = (!pyr_side && !pri_side && !hex_side && !qua_side);

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_tet_quality(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_INT new_node,
                                      REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, cell_node;
  REF_INT node;
  REF_DBL quality, quality0, quality1;
  REF_DBL min_existing_quality;

  *allowed = REF_FALSE;

  min_existing_quality = 1.0;
  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");

    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "q");
    min_existing_quality = MIN(min_existing_quality, quality);
  }

  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node0 == nodes[node]) nodes[node] = new_node;
    }
    RSS(ref_node_tet_quality(ref_node, nodes, &quality0), "q0");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (new_node == nodes[node]) nodes[node] = node0;
    }
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node1 == nodes[node]) nodes[node] = new_node;
    }
    RSS(ref_node_tet_quality(ref_node, nodes, &quality1), "q1");

    if (quality0 < ref_grid_adapt(ref_grid, split_quality_absolute) ||
        quality1 < ref_grid_adapt(ref_grid, split_quality_absolute) ||
        quality0 < ref_grid_adapt(ref_grid, split_quality_relative) *
                       min_existing_quality ||
        quality1 < ref_grid_adapt(ref_grid, split_quality_relative) *
                       min_existing_quality) {
      *allowed = REF_FALSE;
      return REF_SUCCESS;
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_ratio(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                                REF_INT new_node, REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, cell_node;
  REF_INT node;
  REF_INT cell_edge, e0, e1;
  REF_DBL ratio;

  *allowed = REF_FALSE;

  if (ref_grid_surf(ref_grid) || ref_grid_twod(ref_grid)) {
    ref_cell = ref_grid_tri(ref_grid);
  } else {
    ref_cell = ref_grid_tet(ref_grid);
  }

  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node0 == nodes[node]) nodes[node] = new_node;
    }

    each_ref_cell_cell_edge(ref_cell, cell_edge) { /* limit ratio */
      e0 = nodes[ref_cell_e2n_gen(ref_cell, 0, cell_edge)];
      e1 = nodes[ref_cell_e2n_gen(ref_cell, 1, cell_edge)];
      if (e0 == new_node || e1 == new_node) {
        RSS(ref_node_ratio(ref_node, e0, e1, &ratio), "ratio node0");
        if (ratio < ref_grid_adapt(ref_grid, post_min_ratio) ||
            ratio > ref_grid_adapt(ref_grid, post_max_ratio)) {
          *allowed = REF_FALSE;
          return REF_SUCCESS;
        }
      }
    }

    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (new_node == nodes[node]) nodes[node] = node0;
    }
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (node1 == nodes[node]) nodes[node] = new_node;
    }

    each_ref_cell_cell_edge(ref_cell, cell_edge) { /* limit ratio */
      e0 = nodes[ref_cell_e2n_gen(ref_cell, 0, cell_edge)];
      e1 = nodes[ref_cell_e2n_gen(ref_cell, 1, cell_edge)];
      if (e0 == new_node || e1 == new_node) {
        RSS(ref_node_ratio(ref_node, e0, e1, &ratio), "ratio node0");
        if (ratio < ref_grid_adapt(ref_grid, post_min_ratio) ||
            ratio > ref_grid_adapt(ref_grid, post_max_ratio)) {
          *allowed = REF_FALSE;
          return REF_SUCCESS;
        }
      }
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_tri_conformity(REF_GRID ref_grid, REF_INT node0,
                                         REF_INT node1, REF_INT new_node,
                                         REF_BOOL *allowed) {
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, cell_node;
  REF_INT node;
  REF_DBL sign_uv_area, uv_area0, uv_area1;
  REF_DBL normdev, normdev0, normdev1;

  *allowed = REF_FALSE;

  if (0 < ref_geom_n(ref_geom)) {
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
      RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "nd");

      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        if (node0 == nodes[node]) nodes[node] = new_node;
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev0), "nd0");

      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        if (new_node == nodes[node]) nodes[node] = node0;

      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        if (node1 == nodes[node]) nodes[node] = new_node;
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev1), "nd1");

      if ((normdev0 <= normdev &&
           normdev0 < ref_grid_adapt(ref_grid, post_min_normdev)) ||
          (normdev1 <= normdev &&
           normdev1 < ref_grid_adapt(ref_grid, post_min_normdev))) {
        *allowed = REF_FALSE;
        return REF_SUCCESS;
      }
    }
  }

  if (0 < ref_geom_n(ref_geom) && !ref_geom_meshlinked(ref_geom)) {
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
      RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &normdev), "nd");

      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        if (node0 == nodes[node]) nodes[node] = new_node;
      RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area0), "uv area");
      RSS(ref_geom_uv_area_sign(ref_grid, nodes[3], &sign_uv_area), "sign");
      uv_area0 *= sign_uv_area;

      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        if (new_node == nodes[node]) nodes[node] = node0;

      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        if (node1 == nodes[node]) nodes[node] = new_node;
      RSS(ref_geom_uv_area(ref_geom, nodes, &uv_area1), "uv area");
      RSS(ref_geom_uv_area_sign(ref_grid, nodes[3], &sign_uv_area), "sign");
      uv_area1 *= sign_uv_area;

      if (ref_node_min_uv_area(ref_node) > uv_area0 ||
          ref_node_min_uv_area(ref_node) > uv_area1) {
        *allowed = REF_FALSE;
        return REF_SUCCESS;
      }
    }
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_tri_quality(REF_GRID ref_grid, REF_INT node0,
                                      REF_INT node1, REF_INT new_node,
                                      REF_BOOL *allowed) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, cell_node;
  REF_INT node;
  REF_DBL quality, quality0, quality1;
  REF_DBL min_existing_quality;

  *allowed = REF_FALSE;

  min_existing_quality = 1.0;
  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");

    RSS(ref_node_tri_quality(ref_node, nodes, &quality), "q");
    min_existing_quality = MIN(min_existing_quality, quality);
  }

  each_ref_cell_having_node2(ref_cell, node0, node1, item, cell_node, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "cell nodes");

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node0 == nodes[node]) nodes[node] = new_node;
    RSS(ref_node_tri_quality(ref_node, nodes, &quality0), "q0");
    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (new_node == nodes[node]) nodes[node] = node0;

    for (node = 0; node < ref_cell_node_per(ref_cell); node++)
      if (node1 == nodes[node]) nodes[node] = new_node;
    RSS(ref_node_tri_quality(ref_node, nodes, &quality1), "q1");

    if (quality0 < ref_grid_adapt(ref_grid, split_quality_absolute) ||
        quality1 < ref_grid_adapt(ref_grid, split_quality_absolute) ||
        quality0 < ref_grid_adapt(ref_grid, split_quality_relative) *
                       min_existing_quality ||
        quality1 < ref_grid_adapt(ref_grid, split_quality_relative) *
                       min_existing_quality)
      return REF_SUCCESS;
  }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_pattern(REF_GRID ref_grid, REF_INT first,
                                  REF_INT skip) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_INT edge;
  REF_BOOL allowed;
  REF_GLOB global;
  REF_INT new_node;

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  for (edge = first; edge < ref_edge_n(ref_edge); edge += skip) {
    RSS(ref_cell_has_side(ref_grid_tet(ref_grid),
                          ref_edge_e2n(ref_edge, 0, edge),
                          ref_edge_e2n(ref_edge, 1, edge), &allowed),
        "has side");
    if (!allowed) continue;

    RSS(ref_split_edge_mixed(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                             ref_edge_e2n(ref_edge, 1, edge), &allowed),
        "mixed");
    if (!allowed) continue;

    RSS(ref_cell_local_gem(ref_grid_tet(ref_grid), ref_node,
                           ref_edge_e2n(ref_edge, 0, edge),
                           ref_edge_e2n(ref_edge, 1, edge), &allowed),
        "local tet");
    if (!allowed) continue;

    RSS(ref_node_next_global(ref_node, &global), "next global");
    RSS(ref_node_add(ref_node, global, &new_node), "new node");
    RSS(ref_node_interpolate_edge(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                                  ref_edge_e2n(ref_edge, 1, edge), 0.5,
                                  new_node),
        "interp new node");
    RSS(ref_geom_add_between(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                             ref_edge_e2n(ref_edge, 1, edge), 0.5, new_node),
        "geom new node");
    RSS(ref_geom_constrain(ref_grid, new_node), "geom constraint");
    RSS(ref_metric_interpolate_between(
            ref_grid, ref_edge_e2n(ref_edge, 0, edge),
            ref_edge_e2n(ref_edge, 1, edge), new_node),
        "interp new node metric");

    RSS(ref_split_edge(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                       ref_edge_e2n(ref_edge, 1, edge), new_node),
        "split");
  }

  ref_edge_free(ref_edge);

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_geometry(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_INT edge;
  REF_BOOL allowed, tri_side;
  REF_GLOB global;
  REF_INT new_node;

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_cell_has_side(ref_grid_tet(ref_grid),
                          ref_edge_e2n(ref_edge, 0, edge),
                          ref_edge_e2n(ref_edge, 1, edge), &allowed),
        "has side");
    if (!allowed) continue;

    RSS(ref_split_edge_mixed(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                             ref_edge_e2n(ref_edge, 1, edge), &allowed),
        "mixed");
    if (!allowed) continue;

    RSS(ref_cell_local_gem(ref_grid_tet(ref_grid), ref_node,
                           ref_edge_e2n(ref_edge, 0, edge),
                           ref_edge_e2n(ref_edge, 1, edge), &allowed),
        "local tet");
    if (!allowed) continue;

    RSS(ref_cell_has_side(ref_grid_tri(ref_grid),
                          ref_edge_e2n(ref_edge, 0, edge),
                          ref_edge_e2n(ref_edge, 1, edge), &tri_side),
        "has side");
    allowed = (!ref_cell_node_empty(ref_grid_tri(ref_grid),
                                    ref_edge_e2n(ref_edge, 0, edge)) &&
               !ref_cell_node_empty(ref_grid_tri(ref_grid),
                                    ref_edge_e2n(ref_edge, 1, edge)) &&
               !tri_side);
    if (!allowed) continue;

    RSS(ref_node_next_global(ref_node, &global), "next global");
    RSS(ref_node_add(ref_node, global, &new_node), "new node");
    RSS(ref_node_interpolate_edge(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                                  ref_edge_e2n(ref_edge, 1, edge), 0.5,
                                  new_node),
        "interp new node");
    RSS(ref_geom_add_between(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                             ref_edge_e2n(ref_edge, 1, edge), 0.5, new_node),
        "geom new node");
    RSS(ref_geom_constrain(ref_grid, new_node), "geom constraint");
    RSS(ref_metric_interpolate_between(
            ref_grid, ref_edge_e2n(ref_edge, 0, edge),
            ref_edge_e2n(ref_edge, 1, edge), new_node),
        "interp new node metric");

    RSS(ref_split_edge(ref_grid, ref_edge_e2n(ref_edge, 0, edge),
                       ref_edge_e2n(ref_edge, 1, edge), new_node),
        "split");
  }

  ref_edge_free(ref_edge);

  return REF_SUCCESS;
}
