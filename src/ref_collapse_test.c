
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
#include <string.h>

#include "ref_adapt.h"
#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_histogram.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_validation.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  if (argc > 2) {
    REF_GRID ref_grid;
    REF_INT pass;
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[1]), "examine header");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "read grid");
    RSS(ref_part_metric(ref_grid_node(ref_grid), argv[2]), "get metric");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "read metric");
    if (argc > 3) {
      RSS(ref_geom_egads_load(ref_grid_geom(ref_grid), argv[3]),
          "load egads geom");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "load geom");
    }

    RSS(ref_validation_cell_volume(ref_grid), "vol");
    RSS(ref_histogram_quality(ref_grid), "gram");
    RSS(ref_histogram_ratio(ref_grid), "gram");
    ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "stats");

    for (pass = 0; pass < 5; pass++) {
      RSS(ref_collapse_pass(ref_grid), "col pass");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "adapt col");

      RSS(ref_validation_cell_volume(ref_grid), "vol");
      RSS(ref_histogram_quality(ref_grid), "gram");
      RSS(ref_histogram_ratio(ref_grid), "gram");
      ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "stats");
    }

    RSS(ref_grid_free(ref_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* collapse tet into triangle, keep renumbered edge */
    REF_GRID ref_grid;
    REF_INT node0, node1;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    node0 = 0;
    node1 = 3;

    RSS(ref_collapse_edge(ref_grid, node0, node1), "collapse");

    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)), "edg");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* collapse tet into nothing, keep edge */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_INT nodes[REF_CELL_MAX_SIZE_PER];

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    node0 = 0;
    node1 = 1;

    RSS(ref_collapse_edge(ref_grid, node0, node1), "collapse");

    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet");
    REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)), "edg");

    RSS(ref_cell_nodes(ref_grid_edg(ref_grid), 0, nodes), "nodes");
    REIS(0, nodes[0], "1");
    REIS(2, nodes[1], "2");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* collapse tet into nothing, no more edge */
    REF_GRID ref_grid;
    REF_INT node0, node1;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    node0 = 1;
    node1 = 2;

    RSS(ref_collapse_edge(ref_grid, node0, node1), "collapse");

    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet");
    REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(0, ref_cell_n(ref_grid_edg(ref_grid)), "edg");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* collapse removes node */
    REF_GRID ref_grid;
    REF_INT node0, node1;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    node0 = 0;
    node1 = 3;

    RSS(ref_collapse_edge(ref_grid, node0, node1), "collapse");

    REIS(3, ref_node_n(ref_grid_node(ref_grid)), "n");
    REIS(REF_FALSE, ref_node_valid(ref_grid_node(ref_grid), 3), "val");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* collapse tet, renumber triangle, renumber edge */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_INT nodes[REF_CELL_MAX_SIZE_PER];

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    node0 = 3;
    node1 = 1;

    RSS(ref_collapse_edge(ref_grid, node0, node1), "collapse");

    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)), "edg");

    RSS(ref_cell_nodes(ref_grid_tri(ref_grid), 0, nodes), "nodes");
    REIS(0, nodes[0], "0");
    REIS(3, nodes[1], "1");
    REIS(2, nodes[2], "2");

    RSS(ref_cell_nodes(ref_grid_edg(ref_grid), 0, nodes), "nodes");
    REIS(3, nodes[0], "1");
    REIS(2, nodes[1], "2");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* geometry: collapse of volume node? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    RSS(ref_cell_remove(ref_grid_tri(ref_grid), 0), "remove tri");

    node0 = 0;
    node1 = 1;
    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");

    REIS(REF_TRUE, allowed, "interior edge allowed?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* geometry: collapse allowed on face? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_INT nodes[4];
    REF_INT tri1, tri2;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    /*
    2---4
    |\ 1|\
    |0\ | \
    |  \|2 \
    0---1---5
    */
    nodes[0] = 1;
    nodes[1] = 4, nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &tri1), "add tri");
    nodes[0] = 1;
    nodes[1] = 5, nodes[2] = 4;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &tri2), "add tri");

    node0 = 0;
    node1 = 3;
    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");
    REIS(REF_TRUE, allowed, "interior node to face?");

    node0 = 0;
    node1 = 1;
    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");
    REIS(REF_TRUE, allowed, "parallel and interior to face?");

    ref_cell_c2n(ref_grid_tri(ref_grid), 3, tri2) = 20;

    node0 = 4;
    node1 = 1;
    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");
    REIS(REF_TRUE, allowed, "parallel along geom edge?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* geometry: collapse of geom node? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;
    REF_INT node, type, id;
    REF_DBL params[2];

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    RSS(ref_cell_remove(ref_grid_tri(ref_grid), 0), "remove tri");

    node = 1;
    type = REF_GEOM_NODE;
    id = 25;
    RSS(ref_geom_add(ref_grid_geom(ref_grid), node, type, id, params),
        "add node geom");

    node0 = 0;
    node1 = 1;
    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");

    REIS(REF_FALSE, allowed, "interior edge allowed?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* geometry: collapse allowed on edge with jump? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_INT nodes[4];
    REF_INT tri1, tri2;
    REF_BOOL allowed;
    REF_INT node, type, id;
    REF_DBL params[2];

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    /*
    2---4
    |\ 1|\
    |0\ | \   edge 20, 1-2
    |  \|2 \  face all 10
    0---1---5
    */
    nodes[0] = 1;
    nodes[1] = 4, nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &tri1), "add tri");
    nodes[0] = 1;
    nodes[1] = 5, nodes[2] = 4;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &tri2), "add tri");

    node = 1;
    type = REF_GEOM_EDGE;
    id = 20;
    params[0] = 1.0;
    RSS(ref_geom_add(ref_grid_geom(ref_grid), node, type, id, params),
        "add node geom");
    node = 2;
    type = REF_GEOM_EDGE;
    id = 20;
    params[0] = 2.0;
    RSS(ref_geom_add(ref_grid_geom(ref_grid), node, type, id, params),
        "add node geom");
    /* add jump? */

    node0 = 1;
    node1 = 4;
    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");
    REIS(REF_TRUE, allowed, "pull to edge?");

    node0 = 4;
    node1 = 1;
    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");
    REIS(REF_FALSE, allowed, "pull off edge?");

    node0 = 1;
    node1 = 2;
    RSS(ref_collapse_edge_geometry(ref_grid, node0, node1, &allowed),
        "col geom");
    REIS(REF_TRUE, allowed, "parallel along geom edge?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* same normal after collapse? */
    REF_GRID ref_grid;
    REF_INT node;
    REF_INT node0, node1;
    REF_INT nodes[4];
    REF_INT tri1, tri2;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");
    /*
    2---4      y
    |\ 1|\
    |0\ | \    ^
    |  \|2 \   |
    0---1---5  +-> x
    */

    RSS(ref_node_add(ref_grid_node(ref_grid), 4, &node), "add node");
    ref_node_xyz(ref_grid_node(ref_grid), 0, node) = 1.0;
    ref_node_xyz(ref_grid_node(ref_grid), 1, node) = 1.0;
    ref_node_xyz(ref_grid_node(ref_grid), 2, node) = 0.0;

    RSS(ref_node_add(ref_grid_node(ref_grid), 5, &node), "add node");
    ref_node_xyz(ref_grid_node(ref_grid), 0, node) = 2.0;
    ref_node_xyz(ref_grid_node(ref_grid), 1, node) = 0.0;
    ref_node_xyz(ref_grid_node(ref_grid), 2, node) = 0.0;

    nodes[0] = 1;
    nodes[1] = 4, nodes[2] = 2;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &tri1), "add tri");
    nodes[0] = 1;
    nodes[1] = 5, nodes[2] = 4;
    nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &tri2), "add tri");

    node0 = 2;
    node1 = 1;
    RSS(ref_collapse_edge_same_normal(ref_grid, node0, node1, &allowed),
        "norm");
    REIS(REF_TRUE, allowed, "normal will be the same");

    ref_node_xyz(ref_grid_node(ref_grid), 2, node) = 0.5;

    RSS(ref_collapse_edge_same_normal(ref_grid, node0, node1, &allowed),
        "norm");
    REIS(REF_FALSE, allowed, "normal would change");

    ref_node_xyz(ref_grid_node(ref_grid), 0, node) = 2.0;
    ref_node_xyz(ref_grid_node(ref_grid), 1, node) = -1.0;
    ref_node_xyz(ref_grid_node(ref_grid), 2, node) = 0.0;

    node0 = 2;
    node1 = 4;
    RSS(ref_collapse_edge_same_normal(ref_grid, node0, node1, &allowed),
        "norm");
    REIS(REF_FALSE, allowed, "tri will have zero area");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* mixed: collapse tet allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");

    node0 = 0;
    node1 = 1;
    RSS(ref_collapse_edge_mixed(ref_grid, node0, node1, &allowed), "col mixed");

    REIS(REF_TRUE, allowed, "pure tet allowed?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }
  { /* mixed: collapse of/near mixed allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_pri_tet_cap_grid(&ref_grid, ref_mpi), "set up");

    node0 = 5;
    node1 = 6;
    RSS(ref_collapse_edge_mixed(ref_grid, node0, node1, &allowed), "col mixed");
    REIS(REF_TRUE, allowed, "tet near mixed allowed?");

    node0 = 6;
    node1 = 5;
    RSS(ref_collapse_edge_mixed(ref_grid, node0, node1, &allowed), "col mixed");
    REIS(REF_FALSE, allowed, "tet changing mixed allowed?");

    node0 = 3;
    node1 = 4;
    RSS(ref_collapse_edge_mixed(ref_grid, node0, node1, &allowed), "col mixed");
    REIS(REF_FALSE, allowed, "mixed collapse allowed?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* local: collapse allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");

    node0 = 0;
    node1 = 1;
    RSS(ref_collapse_edge_local_cell(ref_grid, node0, node1, &allowed),
        "col loc");
    REIS(REF_TRUE, allowed, "local collapse allowed?");

    ref_node_part(ref_grid_node(ref_grid), 2) =
        ref_node_part(ref_grid_node(ref_grid), 2) + 1;

    node0 = 0;
    node1 = 1;
    RSS(ref_collapse_edge_local_cell(ref_grid, node0, node1, &allowed),
        "col loc");
    REIS(REF_FALSE, allowed, "ghost collapse allowed?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* no collapse, close enough */
    REF_GRID ref_grid;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");

    RSS(ref_collapse_pass(ref_grid), "pass");

    REIS(4, ref_node_n(ref_grid_node(ref_grid)), "nodes");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)), "tets");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* top big */
    REF_GRID ref_grid;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up");

    RSS(ref_node_metric_form(ref_grid_node(ref_grid), 3, 1, 0, 0, 1, 0,
                             1.0 / (10.0 * 10.0)),
        "set top z big");

    RSS(ref_collapse_pass(ref_grid), "pass");

    /* ref_export_by_extension(ref_grid,"ref_collapse_test.tec"); */

    REIS(3, ref_node_n(ref_grid_node(ref_grid)), "nodes");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)), "tets");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* collapse prism, keep qua */
    REF_GRID ref_grid;
    REF_INT keep, remove;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "set up");
    keep = 1;
    remove = 2;

    RSS(ref_collapse_face(ref_grid, keep, remove), "split");

    REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(0, ref_cell_n(ref_grid_edg(ref_grid)), "qua");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* collapse prism and qua */
    REF_GRID ref_grid;
    REF_INT keep, remove;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "set up");
    keep = 0;
    remove = 1;

    RSS(ref_collapse_face(ref_grid, keep, remove), "split");

    REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(0, ref_cell_n(ref_grid_qua(ref_grid)), "qua");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* local prism: collapse allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "set up");

    node0 = 0;
    node1 = 1;
    RSS(ref_collapse_face_local_pris(ref_grid, node0, node1, &allowed),
        "col loc");
    REIS(REF_TRUE, allowed, "local collapse allowed?");

    ref_node_part(ref_grid_node(ref_grid), 3) =
        ref_node_part(ref_grid_node(ref_grid), 3) + 1;

    node0 = 0;
    node1 = 1;
    RSS(ref_collapse_face_local_pris(ref_grid, node0, node1, &allowed),
        "col loc");
    REIS(REF_FALSE, allowed, "ghost collapse allowed?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* geometry: collapse of interior face */
    REF_GRID ref_grid;
    REF_INT keep, remove;
    REF_BOOL allowed;

    RSS(ref_fixture_tri_grid(&ref_grid, ref_mpi), "set up");

    keep = 0;
    remove = 1;
    RSS(ref_collapse_face_geometry(ref_grid, keep, remove, &allowed),
        "col geom");

    REIS(REF_TRUE, allowed, "interior edge allowed?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* edge tangent: collapse allowed? */
    REF_GRID ref_grid;
    REF_INT keep, remove;
    REF_BOOL allowed;

    RSS(ref_fixture_tri2_grid(&ref_grid, ref_mpi), "set up");

    keep = 2;
    remove = 0;
    RSS(ref_collapse_face_same_tangent(ref_grid, keep, remove, &allowed),
        "same");
    REIS(REF_TRUE, allowed, "straight tangent collapse allowed?");

    ref_node_xyz(ref_grid_node(ref_grid), 1, remove) = 0.5;

    RSS(ref_collapse_face_same_tangent(ref_grid, keep, remove, &allowed),
        "same");
    REIS(REF_FALSE, allowed, "curved boundary collapse allowed?");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* no collapse, close enough, twod */
    REF_GRID ref_grid;

    RSS(ref_fixture_tri_grid(&ref_grid, ref_mpi), "set up");

    RSS(ref_collapse_twod_pass(ref_grid), "pass");

    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(0, ref_cell_n(ref_grid_qua(ref_grid)), "qua");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  { /* top big, twod */
    REF_GRID ref_grid;

    RSS(ref_fixture_tri_grid(&ref_grid, ref_mpi), "set up");

    RSS(ref_node_metric_form(ref_grid_node(ref_grid), 1, 1, 0, 0,
                             1.0 / (10.0 * 10.0), 0, 1),
        "set top z big");

    RSS(ref_collapse_twod_pass(ref_grid), "pass");

    REIS(2, ref_node_n(ref_grid_node(ref_grid)), "nodes");
    REIS(0, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "pri");

    RSS(ref_grid_free(ref_grid), "free grid");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
