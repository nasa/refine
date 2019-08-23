
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_cavity.h"

#include "ref_adapt.h"
#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_clump.h"
#include "ref_collapse.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_face.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_twod.h"
#include "ref_validation.h"

#include "ref_histogram.h"
#include "ref_part.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  if (argc > 2) {
    REF_GRID ref_grid;

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[1]), "examine header");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "read grid");
    RSS(ref_geom_egads_load(ref_grid_geom(ref_grid), argv[2]),
        "load egads geom");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "load geom");
    /* RSS(ref_part_metric(ref_grid_node(ref_grid), argv[2]), "get metric"); */
    RSS(ref_metric_interpolated_curvature(ref_grid), "interp curve");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "metric");

    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "stats");

    RSS(ref_cavity_pass(ref_grid), "smooth pass");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt cav");

    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "stats");

    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* add same face 3, raise error */
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];

    RSS(ref_cavity_create(&ref_cavity), "create");

    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = 3;
    RSS(ref_cavity_insert_face(ref_cavity, nodes), "insert first");
    REIS(REF_INVALID, ref_cavity_insert_face(ref_cavity, nodes),
         "insert second");

    nodes[0] = 2;
    nodes[1] = 3;
    nodes[2] = 1;
    REIS(REF_INVALID, ref_cavity_insert_face(ref_cavity, nodes),
         "insert second");

    nodes[0] = 3;
    nodes[1] = 1;
    nodes[2] = 2;
    REIS(REF_INVALID, ref_cavity_insert_face(ref_cavity, nodes),
         "insert second");

    RSS(ref_cavity_free(ref_cavity), "free");
  }

  { /* add opposite face 3, mutual destruction */
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];

    RSS(ref_cavity_create(&ref_cavity), "create");
    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = 3;
    RSS(ref_cavity_insert_face(ref_cavity, nodes), "insert first");
    nodes[0] = 1;
    nodes[1] = 3;
    nodes[2] = 2;
    RSS(ref_cavity_insert_face(ref_cavity, nodes), "insert opposite");

    REIS(0, ref_cavity_nface(ref_cavity), "cancel");

    RSS(ref_cavity_free(ref_cavity), "free");
  }

  { /* find face 3 */
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];
    REF_INT face;
    REF_BOOL reversed;

    RSS(ref_cavity_create(&ref_cavity), "create");
    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = 3;
    RSS(ref_cavity_insert_face(ref_cavity, nodes), "insert first");

    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = 3;
    RSS(ref_cavity_find_face(ref_cavity, nodes, &face, &reversed), "find same");
    REIS(0, face, "found");
    REIS(REF_FALSE, reversed, "not rev");
    nodes[0] = 3;
    nodes[1] = 1;
    nodes[2] = 2;
    RSS(ref_cavity_find_face(ref_cavity, nodes, &face, &reversed), "find same");
    REIS(0, face, "found");
    REIS(REF_FALSE, reversed, "not rev");
    nodes[0] = 2;
    nodes[1] = 3;
    nodes[2] = 1;
    RSS(ref_cavity_find_face(ref_cavity, nodes, &face, &reversed), "find same");
    REIS(0, face, "found");
    REIS(REF_FALSE, reversed, "not rev");

    nodes[0] = 2;
    nodes[1] = 1;
    nodes[2] = 3;
    RSS(ref_cavity_find_face(ref_cavity, nodes, &face, &reversed),
        "find reversed");
    REIS(0, face, "found");
    REIS(REF_TRUE, reversed, "not rev");
    nodes[0] = 3;
    nodes[1] = 2;
    nodes[2] = 1;
    RSS(ref_cavity_find_face(ref_cavity, nodes, &face, &reversed),
        "find reversed");
    REIS(0, face, "found");
    REIS(REF_TRUE, reversed, "not rev");
    nodes[0] = 1;
    nodes[1] = 3;
    nodes[2] = 2;
    RSS(ref_cavity_find_face(ref_cavity, nodes, &face, &reversed),
        "find reversed");
    REIS(0, face, "found");
    REIS(REF_TRUE, reversed, "not rev");

    nodes[0] = 3;
    nodes[1] = 4;
    nodes[2] = 5;
    REIS(REF_NOT_FOUND,
         ref_cavity_find_face(ref_cavity, nodes, &face, &reversed), "missing");
    REIS(REF_EMPTY, face, "found");

    RSS(ref_cavity_free(ref_cavity), "free");
  }

  { /* find seg */
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];
    REF_INT faceid = 10;
    REF_INT seg;
    REF_BOOL reversed;

    RSS(ref_cavity_create(&ref_cavity), "create");
    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = faceid;
    RSS(ref_cavity_insert_seg(ref_cavity, nodes), "insert first");

    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = faceid;
    RSS(ref_cavity_find_seg(ref_cavity, nodes, &seg, &reversed), "find same");
    REIS(0, seg, "found");
    REIS(faceid, ref_cavity_s2n(ref_cavity, 2, seg), "faceid");
    REIS(REF_FALSE, reversed, "not rev");

    nodes[0] = 2;
    nodes[1] = 1;
    nodes[2] = faceid;
    RSS(ref_cavity_find_seg(ref_cavity, nodes, &seg, &reversed),
        "find reversed");
    REIS(0, seg, "found");
    REIS(REF_TRUE, reversed, "not rev");

    nodes[0] = 3;
    nodes[1] = 4;
    nodes[2] = faceid;
    REIS(REF_NOT_FOUND, ref_cavity_find_seg(ref_cavity, nodes, &seg, &reversed),
         "missing");
    REIS(REF_EMPTY, seg, "found");

    RSS(ref_cavity_free(ref_cavity), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* add tet */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];
    REF_INT face;
    REF_BOOL reversed;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "pri");

    /* for parallel, skip test for part with no tets */
    if (ref_cell_valid(ref_grid_tet(ref_grid), 0)) {
      RSS(ref_cavity_create(&ref_cavity), "create");
      RSS(ref_cavity_form_empty(ref_cavity, ref_grid, REF_EMPTY), "form empty");

      RSS(ref_cavity_add_tet(ref_cavity, 0), "insert first");

      if (ref_node_owned(ref_grid_node(ref_grid), 0) &&
          ref_node_owned(ref_grid_node(ref_grid), 1) &&
          ref_node_owned(ref_grid_node(ref_grid), 2) &&
          ref_node_owned(ref_grid_node(ref_grid), 3)) {
        nodes[0] = 0;
        nodes[1] = 1;
        nodes[2] = 2;
        RSS(ref_cavity_find_face(ref_cavity, nodes, &face, &reversed),
            "find 0");
        REIS(REF_FALSE, reversed, "not rev");
      } else {
        REIS(REF_CAVITY_PARTITION_CONSTRAINED, ref_cavity_state(ref_cavity),
             "expected part constraint");
      }

      RSS(ref_cavity_free(ref_cavity), "free");
    }

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* add tri */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];
    REF_INT faceid = 10;
    REF_INT seg;
    REF_BOOL reversed;

    RSS(ref_fixture_tri_surf_grid(&ref_grid, ref_mpi), "pri");

    /* for parallel, skip test for part with no tets */
    if (ref_cell_valid(ref_grid_tri(ref_grid), 0)) {
      RSS(ref_cavity_create(&ref_cavity), "create");
      RSS(ref_cavity_form_empty(ref_cavity, ref_grid, REF_EMPTY), "form empty");

      RSS(ref_cavity_add_tri(ref_cavity, 0), "insert first");

      nodes[0] = 1;
      nodes[1] = 2;
      nodes[2] = faceid;
      RSS(ref_cavity_find_seg(ref_cavity, nodes, &seg, &reversed), "opp 0");
      REIS(REF_FALSE, reversed, "not rev");
      nodes[0] = 2;
      nodes[1] = 0;
      nodes[2] = faceid;
      RSS(ref_cavity_find_seg(ref_cavity, nodes, &seg, &reversed), "op 1");
      REIS(REF_FALSE, reversed, "not rev");
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = faceid;
      RSS(ref_cavity_find_seg(ref_cavity, nodes, &seg, &reversed), "opp 2");
      REIS(REF_FALSE, reversed, "not rev");

      RSS(ref_cavity_free(ref_cavity), "free");
    }

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* insert tet node */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_GLOB global;
    REF_INT node;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "pri");
    ref_node = ref_grid_node(ref_grid);

    if (ref_cell_valid(ref_grid_tet(ref_grid), 0)) {
      RSS(ref_cavity_create(&ref_cavity), "create");
      RSS(ref_node_next_global(ref_node, &global), "next global");
      RSS(ref_node_add(ref_node, global, &node), "new node");
      RSS(ref_cavity_form_empty(ref_cavity, ref_grid, node), "form empty");

      RSS(ref_cavity_add_tet(ref_cavity, 0), "insert first");

      ref_node_xyz(ref_node, 0, node) = 0.1;
      ref_node_xyz(ref_node, 1, node) = 0.2;
      ref_node_xyz(ref_node, 2, node) = 0.3;

      RSS(ref_cavity_replace(ref_cavity), "replace");

      REIS(5, ref_node_n(ref_grid_node(ref_grid)), "nodes");
      REIS(4, ref_cell_n(ref_grid_tet(ref_grid)), "cells");

      RSS(ref_cavity_free(ref_cavity), "free");
    }

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* visible three node face */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_GLOB global;
    REF_INT node, face;
    REF_BOOL visible;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "pri");
    ref_node = ref_grid_node(ref_grid);

    if (ref_cell_valid(ref_grid_tet(ref_grid), 0)) {
      RSS(ref_cavity_create(&ref_cavity), "create");
      RSS(ref_node_next_global(ref_node, &global), "next global");
      RSS(ref_node_add(ref_node, global, &node), "new node");
      RSS(ref_cavity_form_empty(ref_cavity, ref_grid, node), "form empty");

      RSS(ref_cavity_add_tet(ref_cavity, 0), "insert first");

      ref_node_xyz(ref_node, 0, node) = 0.1;
      ref_node_xyz(ref_node, 1, node) = 0.2;
      ref_node_xyz(ref_node, 2, node) = 0.3;
      face = 0;
      RSS(ref_cavity_visible(ref_cavity, face, &visible), "viz");
      REIS(REF_TRUE, visible, "vis");

      ref_node_xyz(ref_node, 0, node) = 1.0;
      ref_node_xyz(ref_node, 1, node) = 1.0;
      ref_node_xyz(ref_node, 2, node) = 1.0;
      face = 0;
      RSS(ref_cavity_visible(ref_cavity, face, &visible), "viz");
      REIS(REF_FALSE, visible, "vis");

      RSS(ref_cavity_free(ref_cavity), "free");
    }

    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* enlarge threed face */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;

    RSS(ref_fixture_tet2_grid(&ref_grid, ref_mpi), "brick");

    RSS(ref_cavity_create(&ref_cavity), "create");
    RSS(ref_cavity_form_empty(ref_cavity, ref_grid, REF_EMPTY), "form empty");
    RSS(ref_cavity_add_tet(ref_cavity, 0), "insert first tri");
    REIS(4, ref_cavity_nface(ref_cavity), "n");
    REIS(1, ref_list_n(ref_cavity_tet_list(ref_cavity)), "l");
    RSS(ref_cavity_enlarge_face(ref_cavity, 0), "enl face 1");
    REIS(6, ref_cavity_nface(ref_cavity), "n");
    REIS(2, ref_list_n(ref_cavity_tet_list(ref_cavity)), "l");
    RSS(ref_cavity_free(ref_cavity), "free");
    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split edge of tet */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;
    REF_INT new_node;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "pri");
    RSS(ref_cavity_create(&ref_cavity), "create");

    RSS(ref_node_add(ref_grid_node(ref_grid), 4, &new_node), "new");
    RSS(ref_cavity_form_edge_split(ref_cavity, ref_grid, 0, 3, new_node),
        "insert edge");
    REIS(2, ref_cavity_nface(ref_cavity), "n");
    REIS(1, ref_list_n(ref_cavity_tet_list(ref_cavity)), "l");

    RSS(ref_cavity_free(ref_cavity), "free");
    RSS(ref_grid_free(ref_grid), "free");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split edge of tet2 */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;
    REF_INT new_node;

    RSS(ref_fixture_tet2_grid(&ref_grid, ref_mpi), "pri");
    RSS(ref_cavity_create(&ref_cavity), "create");

    RSS(ref_node_add(ref_grid_node(ref_grid), 5, &new_node), "new");
    RSS(ref_cavity_form_edge_split(ref_cavity, ref_grid, 1, 2, new_node),
        "insert edge");
    REIS(4, ref_cavity_nface(ref_cavity), "n");
    REIS(2, ref_list_n(ref_cavity_tet_list(ref_cavity)), "l");

    RSS(ref_cavity_free(ref_cavity), "free");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* tet brick insert (on every part) */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT node, nnode;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_node = ref_grid_node(ref_grid);
    nnode = ref_node_n(ref_node);

    node = 39;
    ref_node_xyz(ref_node, 0, node) = 0.5;
    RSS(ref_cavity_create(&ref_cavity), "create");
    RSS(ref_cavity_form_ball(ref_cavity, ref_grid, node), "insert first");
    RSS(ref_cavity_enlarge_visible(ref_cavity), "insert first");
    REIS(REF_CAVITY_VISIBLE, ref_cavity_state(ref_cavity),
         "enlarge not successful");
    RSS(ref_cavity_replace(ref_cavity), "free");
    RSS(ref_cavity_free(ref_cavity), "free");

    RAS(nnode > ref_node_n(ref_node), "node count did not decrease");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* tet brick can not insert outside of boundary (on every part) */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT node;

    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
    ref_node = ref_grid_node(ref_grid);

    node = 39;
    ref_node_xyz(ref_node, 0, node) = -0.5;
    RSS(ref_cavity_create(&ref_cavity), "create");
    RSS(ref_cavity_form_ball(ref_cavity, ref_grid, node), "insert first");
    RSS(ref_cavity_enlarge_visible(ref_cavity), "enlarge");
    REIS(REF_CAVITY_BOUNDARY_CONSTRAINED, ref_cavity_state(ref_cavity),
         "enlarge wrong state");
    RSS(ref_cavity_free(ref_cavity), "free");

    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* replace tet */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "pri");
    RSS(ref_cavity_create(&ref_cavity), "create");

    if (!ref_mpi_para(ref_mpi)) {
      RSS(ref_cavity_form_ball(ref_cavity, ref_grid, 0), "insert ball");

      REIS(2, ref_cavity_nface(ref_cavity), "n");
      REIS(1, ref_list_n(ref_cavity_tet_list(ref_cavity)), "l");

      RSS(ref_cavity_replace(ref_cavity), "replace");
    }

    RSS(ref_cavity_free(ref_cavity), "free");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* replace surf tri */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "pri");
    RSS(ref_cavity_create(&ref_cavity), "create");

    if (!ref_mpi_para(ref_mpi)) {
      RSS(ref_cavity_form_empty(ref_cavity, ref_grid, 0), "insert ball");
      RSS(ref_cavity_add_tri(ref_cavity, 0), "insert tri");

      REIS(0, ref_cavity_nface(ref_cavity), "n");
      REIS(0, ref_list_n(ref_cavity_tet_list(ref_cavity)), "l");
      REIS(3, ref_cavity_nseg(ref_cavity), "n");
      REIS(1, ref_list_n(ref_cavity_tri_list(ref_cavity)), "l");

      RSS(ref_cavity_replace(ref_cavity), "replace");
    }

    RSS(ref_cavity_free(ref_cavity), "free");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
