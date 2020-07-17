
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

#include "ref_subdiv.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_export.h"
#include "ref_face.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_sort.h"
#include "ref_validation.h"

static REF_STATUS set_up_tet_for_subdiv(REF_SUBDIV *ref_subdiv_ptr,
                                        REF_MPI ref_mpi) {
  REF_GRID ref_grid;

  RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "tet");
  RSS(ref_subdiv_create(ref_subdiv_ptr, ref_grid), "create");

  return REF_SUCCESS;
}

static REF_STATUS set_up_pyramid_for_subdiv(REF_SUBDIV *ref_subdiv_ptr,
                                            REF_MPI ref_mpi) {
  REF_GRID ref_grid;

  RSS(ref_fixture_pyr_grid(&ref_grid, ref_mpi), "pri");
  RSS(ref_subdiv_create(ref_subdiv_ptr, ref_grid), "create");

  return REF_SUCCESS;
}

static REF_STATUS set_up_prism_for_subdiv(REF_SUBDIV *ref_subdiv_ptr,
                                          REF_MPI ref_mpi) {
  REF_GRID ref_grid;

  RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "pri");
  RSS(ref_subdiv_create(ref_subdiv_ptr, ref_grid), "create");

  return REF_SUCCESS;
}

static REF_STATUS tear_down(REF_SUBDIV ref_subdiv) {
  REF_GRID ref_grid;

  ref_grid = ref_subdiv_grid(ref_subdiv);

  RSS(ref_subdiv_free(ref_subdiv), "free");

  RSS(ref_grid_free(ref_grid), "free");

  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  if (2 == argc) {
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_EDGE ref_edge;
    REF_INT edge;

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[1]), "import");
    RSS(ref_migrate_to_balance(ref_grid), "new part");
    RSS(ref_subdiv_create(&ref_subdiv, ref_grid), "create");

    ref_node = ref_grid_node(ref_grid);
    ref_edge = ref_subdiv_edge(ref_subdiv);

    if (ref_mpi_once(ref_mpi))
      printf("orig " REF_GLOB_FMT "\n", ref_node_n_global(ref_node));

    for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
      REF_DBL zmin = 0.00005;
      if (ref_node_xyz(ref_node, 2, ref_edge_e2n(ref_edge, 0, edge)) < zmin &&
          ref_node_xyz(ref_node, 2, ref_edge_e2n(ref_edge, 1, edge)) < zmin)
        RSS(ref_subdiv_mark_to_split(ref_subdiv,
                                     ref_edge_e2n(ref_edge, 0, edge),
                                     ref_edge_e2n(ref_edge, 1, edge)),
            "mark edge");
    }

    if (!ref_mpi_para(ref_mpi)) {
      RSS(ref_subdiv_mark_relax(ref_subdiv), "relax");
      ref_edge_tec_int(ref_subdiv_edge(ref_subdiv), "edge.tec",
                       ref_subdiv->mark);
      RSS(ref_subdiv_mark_verify(ref_subdiv), "vrfy");
    }

    RSS(ref_subdiv_split(ref_subdiv), "split");

    if (ref_mpi_once(ref_mpi))
      printf("split " REF_GLOB_FMT "\n", ref_node_n_global(ref_node));

    RSS(ref_gather_by_extension(ref_grid, "ref_subdiv_test.b8.ugrid"),
        "gather");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (3 == argc) {
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;

    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "import");
    RSS(ref_subdiv_create(&ref_subdiv, ref_grid), "create");

    RSS(ref_subdiv_mark_prism_sides(ref_subdiv), "mark sides");

    if (ref_mpi_once(ref_mpi)) ref_grid_inspect(ref_grid);

    RSS(ref_subdiv_split(ref_subdiv), "split");

    if (ref_mpi_once(ref_mpi)) ref_grid_inspect(ref_grid);

    RSS(ref_export_by_extension(ref_grid, argv[2]), "export");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 4) { /* split stack in two */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node0, node1;

    RSS(ref_fixture_pri_stack_grid(&ref_grid, ref_mpi), "stack");
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_subdiv_create(&ref_subdiv, ref_grid), "create");

    REIS(12, ref_node_n_global(ref_node), "start with 12");

    if (REF_SUCCESS == ref_node_local(ref_node, 0, &node0) &&
        REF_SUCCESS == ref_node_local(ref_node, 1, &node1))
      if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, node0))
        RSS(ref_subdiv_mark_to_split(ref_subdiv, node0, node1), "mark edge");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(16, ref_node_n_global(ref_node), "where my nodes?");

    RSS(ref_validation_cell_node(ref_grid), "validate");
    RSS(ref_validation_unused_node(ref_grid), "validate");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 6) { /* mark and relax prism */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RES(0, ref_subdiv_mark(ref_subdiv, 0), "init mark");
    RES(0, ref_subdiv_mark(ref_subdiv, 6), "init mark");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0-1");

    RES(1, ref_subdiv_mark(ref_subdiv, 0), "been marked");
    RES(0, ref_subdiv_mark(ref_subdiv, 6), "not been marked");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    RES(1, ref_subdiv_mark(ref_subdiv, 0), "been marked");
    RES(1, ref_subdiv_mark(ref_subdiv, 6), "been relaxed");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 6) { /* mark and relax prism */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge 1-2");

    RES(0, ref_subdiv_mark(ref_subdiv, 1), "no yet");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    RES(1, ref_subdiv_mark(ref_subdiv, 1), "yet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 4) { /* relax tet */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge 1-2");

    REIS(1, ref_subdiv_mark(ref_subdiv, 0), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 3), "yet");

    REIS(0, ref_subdiv_mark(ref_subdiv, 1), "no yet");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_subdiv_mark(ref_subdiv, 1), "yet");

    REIS(0, ref_subdiv_mark(ref_subdiv, 2), "no yet");
    REIS(0, ref_subdiv_mark(ref_subdiv, 4), "no yet");
    REIS(0, ref_subdiv_mark(ref_subdiv, 5), "no yet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 4) { /* relax tet opposite edges */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 3), "mark edge 5");

    REIS(1, ref_subdiv_mark(ref_subdiv, 0), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 5), "yet");

    REIS(0, ref_subdiv_mark(ref_subdiv, 1), "no yet");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_subdiv_mark(ref_subdiv, 1), "yet");

    REIS(1, ref_subdiv_mark(ref_subdiv, 2), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 4), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 5), "yet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 4) { /* unrelax 2 to 1, first two edge */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 3), "mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 1), "mark edge");

    REIS(1, ref_subdiv_mark(ref_subdiv, 5), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 4), "yet");

    RSS(ref_subdiv_unmark_relax(ref_subdiv), "split");

    if (!ref_mpi_para(ref_mpi)) {
      REIS(0, ref_subdiv_mark(ref_subdiv, 5), "yet");
      REIS(1, ref_subdiv_mark(ref_subdiv, 4), "yet");
    }

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 4) { /* unrelax 2 to 1, first and last edge */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 3), "mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge");

    REIS(1, ref_subdiv_mark(ref_subdiv, 5), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 3), "yet");

    RSS(ref_subdiv_unmark_relax(ref_subdiv), "split");

    REIS(0, ref_subdiv_mark(ref_subdiv, 5), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 3), "yet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 4) { /* unrelax 2 to 1, last two edges */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 1), "mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge");

    REIS(1, ref_subdiv_mark(ref_subdiv, 4), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 3), "yet");

    RSS(ref_subdiv_unmark_relax(ref_subdiv), "split");

    REIS(0, ref_subdiv_mark(ref_subdiv, 4), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 3), "yet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 4) { /* unrelax para 0 2 */
    REF_SUBDIV ref_subdiv;
    REF_BOOL again = REF_FALSE;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 3), "mark edge");

    REIS(1, ref_subdiv_mark(ref_subdiv, 0), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 2), "yet");

    RSS(ref_subdiv_unmark_tet(ref_subdiv, 0, &again), "unmark");

    REIS(1, ref_subdiv_mark(ref_subdiv, 0), "yet");
    REIS(0, ref_subdiv_mark(ref_subdiv, 2), "yet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 4) { /* unrelax para 0 5 */
    REF_SUBDIV ref_subdiv;
    REF_BOOL again = REF_FALSE;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 3), "mark edge");

    REIS(1, ref_subdiv_mark(ref_subdiv, 0), "yet");
    REIS(1, ref_subdiv_mark(ref_subdiv, 5), "yet");

    RSS(ref_subdiv_unmark_tet(ref_subdiv, 0, &again), "unmark");

    REIS(1, ref_subdiv_mark(ref_subdiv, 0), "yet");
    REIS(0, ref_subdiv_mark(ref_subdiv, 5), "yet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 6) { /* new prism nodes */
    REF_SUBDIV ref_subdiv;
    REF_NODE ref_node;
    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 4), "mark edge 6 3-4");

    REIS(6, ref_node_n(ref_node), "6 node prism");
    RSS(ref_subdiv_split(ref_subdiv), "split");

    if (!ref_mpi_para(ref_mpi)) {
      REIS(8, ref_node_n(ref_node), "two new nodes");
      REIS(6, ref_subdiv_node(ref_subdiv, 0), "new 6");
      REIS(7, ref_subdiv_node(ref_subdiv, 6), "new 6");
    }

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 6) { /* split prism in two */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0-1");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    if (!ref_mpi_para(ref_mpi)) {
      REIS(2, ref_cell_n(ref_grid_pri(ref_grid)), "two");
    }

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 6) { /* split prism in two */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0-1");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    if (!ref_mpi_para(ref_mpi)) {
      REIS(2, ref_cell_n(ref_grid_pri(ref_grid)), "two pri");
      REIS(0, ref_cell_n(ref_grid_qua(ref_grid)), "no tri");
      REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "no qua");
    }

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (ref_mpi_n(ref_mpi) <= 6) { /* split prism in four */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge 1-2");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 0), "mark edge 2-0");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    if (!ref_mpi_para(ref_mpi)) {
      REIS(4, ref_cell_n(ref_grid_pri(ref_grid)), "two pri");
      REIS(0, ref_cell_n(ref_grid_qua(ref_grid)), "no tri");
      REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "no qua");
    }

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* cleave prism across quads */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_prism_sides(ref_subdiv), "sides");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    if (!ref_mpi_para(ref_mpi)) {
      REIS(2, ref_cell_n(ref_grid_pri(ref_grid)), "two");
      REIS(9, ref_node_n(ref_grid_node(ref_grid)), "9");
    }

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split tet in two, map 1 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0-1");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(2, ref_cell_n(ref_grid_tet(ref_grid)), "two tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)), "two tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)), "still one edg");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split tet in two, map 4 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 3), "mark edge 0-3");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(2, ref_cell_n(ref_grid_tet(ref_grid)), "two tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)), "still one tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)), "still one edg");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split tet in 4 around node 3 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge 1-2");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)), "four tet");
    REIS(4, ref_cell_n(ref_grid_tri(ref_grid)), "four tri");
    REIS(2, ref_cell_n(ref_grid_edg(ref_grid)), "two edg");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split tet in 4 around node 0 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 2), "mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)), "four tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(2, ref_cell_n(ref_grid_edg(ref_grid)), "two edg");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split tet in 4 around node 1 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 2), "mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 2), "mark edge");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)), "four tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)), "still one edg");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split tet in 4 around node 2 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 1), "mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)), "four tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)), "still one edg");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  if (!ref_mpi_para(ref_mpi)) { /* split tet in 8 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 3), "mark edge 5");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(8, ref_cell_n(ref_grid_tet(ref_grid)), "eight tet");
    REIS(4, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(2, ref_cell_n(ref_grid_edg(ref_grid)), "two edg");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* unsplit pyramid */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "tri");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* unsplit pyramid */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* relax and split pyramid in two, mark e0:n0n1, promoting opp */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_subdiv_mark(ref_subdiv, 7), "promoted");

    REIS(2, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* relax and split pyramid in two, mark e2:n0n3, promoting opp */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 3), "mark edge 2");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_subdiv_mark(ref_subdiv, 4), "promoted");

    REIS(2, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* relax and split pyramid in two, mark e6:n2n4, no promotion */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 4), "mark edge 6");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(2, ref_cell_n(ref_grid_tet(ref_grid)), "tet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* relax and split pyramid in to pyr and pri, mark e6, e3 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge 3");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 4, 2), "mark edge 6");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(1, ref_cell_n(ref_grid_pri(ref_grid)), "pri");
    REIS(6, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(2, ref_cell_n(ref_grid_qua(ref_grid)), "qua");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* relax and split pyramid in to pyr and pri e1, e5 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 2), "mark edge 1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 2), "mark edge 5");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(1, ref_cell_n(ref_grid_pri(ref_grid)), "pri");
    REIS(6, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(2, ref_cell_n(ref_grid_qua(ref_grid)), "qua");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* relax and split pyramid in to pyr and 3 pri, all but e2,e4 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 2), "mark edge 1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge 3");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 2), "mark edge 5");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 4, 2), "mark edge 6");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 4), "mark edge 7");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(3, ref_cell_n(ref_grid_pri(ref_grid)), "pri");
    REIS(10, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(4, ref_cell_n(ref_grid_qua(ref_grid)), "qua");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* split pyramid in to pyrs and tets, e0e1e3 base and e7 top */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 2), "mark edge 1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge 3");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 4), "mark edge 7");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(2, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)), "tet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* split pyramid in to pyrs and tets, e5e6e7 base and e0 top */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 3), "mark edge 5");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 4), "mark edge 6");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 3, 4), "mark edge 7");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 1), "mark edge 0");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(2, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)), "tet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* split pyramid in to pyrs and tets, e3e4e6 base and e2 top */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 2), "mark edge 3");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 1, 4), "mark edge 4");
    RSS(ref_subdiv_mark_to_split(ref_subdiv, 2, 4), "mark edge 6");

    RSS(ref_subdiv_mark_to_split(ref_subdiv, 0, 3), "mark edge 2");

    RSS(ref_subdiv_split(ref_subdiv), "split");

    REIS(2, ref_cell_n(ref_grid_pyr(ref_grid)), "pyr");
    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)), "tet");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* split prism by metric, o.k. metric */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_INT n;

    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);
    ref_node = ref_grid_node(ref_grid);

    each_ref_node_valid_node(ref_node, node) {
      RSS(ref_node_metric_form(ref_node, node, 1.0, 0, 0, 1.0 / (0.8 * 0.8), 0,
                               1),
          "add");
    }

    RSS(ref_subdiv_mark_prism_by_metric(ref_subdiv), "mark metric");

    RSS(ref_subdiv_mark_relax(ref_subdiv), "relax");

    RSS(ref_subdiv_mark_n(ref_subdiv, &n), "relax");

    REIS(0, n, "marks");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  { /* split prism by metric, x-split metric */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;
    REF_INT n;

    RSS(set_up_prism_for_subdiv(&ref_subdiv, ref_mpi), "set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);
    ref_node = ref_grid_node(ref_grid);

    each_ref_node_valid_node(ref_node, node) {
      RSS(ref_node_metric_form(ref_node, node, 1.0 / (0.8 * 0.8), 0, 0,
                               1.0 / (0.8 * 0.8), 0, 1),
          "add");
    }

    RSS(ref_subdiv_mark_prism_by_metric(ref_subdiv), "mark metric");

    RSS(ref_subdiv_mark_relax(ref_subdiv), "relax");

    RSS(ref_subdiv_mark_n(ref_subdiv, &n), "relax");

    REIS(2, n, "marks");

    RSS(tear_down(ref_subdiv), "tear down");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
