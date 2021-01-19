
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

#include "ref_validation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_export.h"
#include "ref_face.h"
#include "ref_malloc.h"
#include "ref_mpi.h"
#include "ref_sort.h"

REF_STATUS ref_validation_simplex_node(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_BOOL problem;

  ref_cell = ref_grid_tet(ref_grid);
  if (ref_grid_twod(ref_grid) || ref_grid_surf(ref_grid))
    ref_cell = ref_grid_tri(ref_grid);

  problem = REF_FALSE;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_cell_node_empty(ref_cell, node)) {
      problem = REF_TRUE;
      RSS(ref_node_location(ref_node, node), "location");
    }
  }

  return (problem ? REF_FAILURE : REF_SUCCESS);
}

REF_STATUS ref_validation_unused_node(REF_GRID ref_grid) {
  REF_INT node;
  REF_BOOL problem;
  REF_ADJ ref_adj;
  REF_CELL ref_cell;
  REF_INT group, cell, nodes[REF_CELL_MAX_SIZE_PER];

  problem = REF_FALSE;
  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    if (ref_adj_empty(ref_cell_adj(ref_grid_tet(ref_grid)), node) &&
        ref_adj_empty(ref_cell_adj(ref_grid_pyr(ref_grid)), node) &&
        ref_adj_empty(ref_cell_adj(ref_grid_pri(ref_grid)), node) &&
        ref_adj_empty(ref_cell_adj(ref_grid_hex(ref_grid)), node)) {
      problem = REF_TRUE;
      printf(" unused node %d: %e %e %e\n", node,
             ref_node_xyz(ref_grid_node(ref_grid), 0, node),
             ref_node_xyz(ref_grid_node(ref_grid), 1, node),
             ref_node_xyz(ref_grid_node(ref_grid), 2, node));
    }
  }

  ref_adj_create(&ref_adj);

  each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        RSS(ref_adj_add(ref_adj, nodes[node], group + 4 * cell), "add");
    }
  }

  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    if (ref_adj_empty(ref_adj, node)) {
      problem = REF_TRUE;
      printf(" unused node %d\n", node);
    }
  }

  ref_adj_free(ref_adj);

  return (problem ? REF_FAILURE : REF_SUCCESS);
}

REF_STATUS ref_validation_boundary_at_node(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_edge, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT i, c, ns[REF_CELL_MAX_SIZE_PER];
  REF_INT node0, node1, ncell, cell_list[10];

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "get cell");
    each_ref_cell_cell_edge(ref_cell, cell_edge) {
      node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
      node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
      if (!ref_node_owned(ref_node, node0) && !ref_node_owned(ref_node, node1))
        continue;
      RSB(ref_cell_list_with2(ref_cell, node0, node1, 10, &ncell, cell_list),
          "more than ten triangles found with two nodes", {
            printf("nodes %d %d faceid %d\n", node0, node1, nodes[3]);
            ref_node_location(ref_node, node0);
            ref_node_location(ref_node, node1);
          });
      REIB(2, ncell, "two triangles expected", {
        printf("nodes %d %d faceid %d\n", node0, node1, nodes[3]);
        ref_node_location(ref_node, node0);
        ref_node_location(ref_node, node1);
        for (i = 0; i < ncell; i++) {
          c = cell_list[i];
          RSS(ref_cell_nodes(ref_cell, c, ns), "get cell");
          printf(" item %d cell %d nodes %d %d %d id %d\n", i, c, ns[0], ns[1],
                 ns[2], ns[3]);
        }
      });
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_validation_boundary_manifold(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell, cell_edge, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT i, c, ns[REF_CELL_MAX_SIZE_PER];
  REF_INT node0, node1, ncell, cell_list[10];

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    each_ref_cell_cell_edge(ref_cell, cell_edge) {
      node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
      node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
      if (!ref_node_owned(ref_node, node0) && !ref_node_owned(ref_node, node1))
        continue;
      RSB(ref_cell_list_with2(ref_cell, node0, node1, 10, &ncell, cell_list),
          "more than ten triangles found with two nodes", {
            printf("nodes %d %d faceid %d\n", node0, node1, nodes[3]);
            ref_node_location(ref_node, node0);
            ref_node_location(ref_node, node1);
          });
      REIB(2, ncell, "two triangles expected", {
        printf("nodes %d %d faceid %d\n", node0, node1, nodes[3]);
        ref_node_location(ref_node, node0);
        ref_node_location(ref_node, node1);
        for (i = 0; i < ncell; i++) {
          c = cell_list[i];
          RSS(ref_cell_nodes(ref_cell, c, ns), "get cell");
          printf(" item %d cell %d nodes %d %d %d id %d\n", i, c, ns[0], ns[1],
                 ns[2], ns[3]);
        }
      });
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_validation_boundary_face(REF_GRID ref_grid) {
  REF_CELL ref_cell;
  REF_BOOL has_face;
  REF_INT cell;
  REF_INT node;
  REF_INT nodes[4];
  REF_BOOL problem;

  problem = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    for (node = 0; node < 3; node++)
      nodes[node] = ref_cell_c2n(ref_cell, node, cell);
    nodes[3] = nodes[0];
    RSS(ref_grid_cell_has_face(ref_grid, nodes, &has_face), "has_face");
    if (!has_face) {
      problem = REF_TRUE;
      printf("triangle %d nodes %d %d %d global " REF_GLOB_FMT " " REF_GLOB_FMT
             " " REF_GLOB_FMT "\n",
             cell, nodes[0], nodes[1], nodes[2],
             ref_node_global(ref_grid_node(ref_grid), nodes[0]),
             ref_node_global(ref_grid_node(ref_grid), nodes[1]),
             ref_node_global(ref_grid_node(ref_grid), nodes[2]));
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[0]), "n0");
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[1]), "n1");
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[2]), "n2");
    }
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    for (node = 0; node < 4; node++)
      nodes[node] = ref_cell_c2n(ref_cell, node, cell);
    RSS(ref_grid_cell_has_face(ref_grid, nodes, &has_face), "has_face");
    if (!has_face) problem = REF_TRUE;
  }

  return (problem ? REF_FAILURE : REF_SUCCESS);
}

REF_STATUS ref_validation_boundary_all(REF_GRID ref_grid) {
  RSS(ref_validation_boundary_face(ref_grid), "face");
  RSS(ref_validation_boundary_manifold(ref_grid), "manifold");
  return REF_SUCCESS;
}

static REF_STATUS ref_validation_node_cell(REF_CELL ref_cell, REF_INT node) {
  REF_INT item, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  printf("node %d\n", node);
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    printf("cell %d:", cell);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "get cell");
    for (cell_node = 0; cell_node < ref_cell_node_per(ref_cell); cell_node++) {
      printf(" %d", nodes[cell_node]);
    }
    printf("\n");
  }
  return REF_SUCCESS;
}

REF_STATUS ref_validation_cell_face_node(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tet = ref_grid_tet(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_face;
  REF_INT face_node, face_nodes[4];
  REF_INT cell0, cell1, tri_cell;

  each_ref_cell_having_node(tet, node, item, cell) {
    each_ref_cell_cell_face(tet, cell_face) {
      for (face_node = 0; face_node < 4; face_node++) {
        face_nodes[face_node] = ref_cell_f2n(tet, face_node, cell_face, cell);
      }
      if (ref_node_owned(ref_node, face_nodes[0]) ||
          ref_node_owned(ref_node, face_nodes[1]) ||
          ref_node_owned(ref_node, face_nodes[2]) ||
          ref_node_owned(ref_node, face_nodes[3])) {
        RSB(ref_cell_with_face(tet, face_nodes, &cell0, &cell1), "two expected",
            {
              ref_node_location(ref_node, face_nodes[0]);
              ref_node_location(ref_node, face_nodes[1]);
              ref_node_location(ref_node, face_nodes[2]);
              ref_node_location(ref_node, face_nodes[3]);
              ref_validation_node_cell(tet, face_nodes[0]);
              ref_validation_node_cell(tet, face_nodes[1]);
              ref_validation_node_cell(tet, face_nodes[2]);
            });
        if (REF_EMPTY == cell1) { /* should have tri on other side */
          RSB(ref_cell_with(tri, face_nodes, &tri_cell), "matching tri", {
            ref_node_location(ref_node, face_nodes[0]);
            ref_node_location(ref_node, face_nodes[1]);
            ref_node_location(ref_node, face_nodes[2]);
            ref_node_location(ref_node, face_nodes[3]);
          });
        }
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_validation_cell_face(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_FACE ref_face;
  REF_CELL ref_cell;
  REF_INT *hits;
  REF_INT face;
  REF_INT group, cell, cell_face;
  REF_INT node;
  REF_INT nodes[4];
  REF_BOOL problem, report;
  REF_STATUS code;

  problem = REF_FALSE;

  RSS(ref_face_create(&ref_face, ref_grid), "face");

  ref_malloc(hits, ref_face_n(ref_face), REF_INT);

  for (face = 0; face < ref_face_n(ref_face); face++) hits[face] = 0;

  each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell) {
    each_ref_cell_valid_cell(ref_cell, cell) {
      each_ref_cell_cell_face(ref_cell, cell_face) {
        for (node = 0; node < 4; node++) {
          nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);
        }
        RSS(ref_face_with(ref_face, nodes, &face), "find cell face");
        hits[face]++;
      }
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    for (node = 0; node < 3; node++) {
      nodes[node] = ref_cell_c2n(ref_cell, node, cell);
    }
    nodes[3] = nodes[0];
    code = ref_face_with(ref_face, nodes, &face);
    if (REF_SUCCESS != code) {
      ref_node_location(ref_grid_node(ref_grid), nodes[0]);
      ref_node_location(ref_grid_node(ref_grid), nodes[1]);
      ref_node_location(ref_grid_node(ref_grid), nodes[2]);
      ref_node_location(ref_grid_node(ref_grid), nodes[3]);
    }
    RSS(code, "find tri");
    hits[face]++;
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    for (node = 0; node < 4; node++) {
      nodes[node] = ref_cell_c2n(ref_cell, node, cell);
    }
    RSS(ref_face_with(ref_face, nodes, &face), "find qua");
    hits[face]++;
  }

  for (face = 0; face < ref_face_n(ref_face); face++) {
    report = REF_FALSE;
    if (ref_mpi_para(ref_grid_mpi(ref_grid))) {
      report = report || (2 < hits[face]);
      if (ref_node_owned(ref_node, ref_face_f2n(ref_face, 0, face)) ||
          ref_node_owned(ref_node, ref_face_f2n(ref_face, 1, face)) ||
          ref_node_owned(ref_node, ref_face_f2n(ref_face, 2, face)) ||
          ref_node_owned(ref_node, ref_face_f2n(ref_face, 3, face))) {
        report = report || (2 > hits[face]);
      }
    } else {
      report = report || (2 != hits[face]);
    }
    if (report) {
      problem = REF_TRUE;
      printf(" hits %d\n", hits[face]);
      for (node = 0; node < 4; node++) {
        nodes[node] = ref_face_f2n(ref_face, node, face);
      }
      printf("face %d nodes %d %d %d %d global " REF_GLOB_FMT " " REF_GLOB_FMT
             " " REF_GLOB_FMT " " REF_GLOB_FMT "\n",
             face, nodes[0], nodes[1], nodes[2], nodes[3],
             ref_node_global(ref_grid_node(ref_grid), nodes[0]),
             ref_node_global(ref_grid_node(ref_grid), nodes[1]),
             ref_node_global(ref_grid_node(ref_grid), nodes[2]),
             ref_node_global(ref_grid_node(ref_grid), nodes[3]));
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[0]), "n0");
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[1]), "n1");
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[2]), "n2");
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[3]), "n3");
    }
  }

  free(hits);

  RSS(ref_face_free(ref_face), "face free");

  return (problem ? REF_FAILURE : REF_SUCCESS);
}

REF_STATUS ref_validation_cell_node(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT group;
  REF_INT cell, node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL has_local;

  each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      has_local = REF_FALSE;
      for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
        if (!ref_node_valid(ref_grid_node(ref_grid), nodes[node])) {
          RSB(REF_FAILURE, "cell with invalid node", {
            printf("group %d node_per %d\n", group,
                   ref_cell_node_per(ref_cell));
          });
        }
        has_local = has_local || (ref_mpi_rank(ref_grid_mpi(ref_grid)) ==
                                  ref_node_part(ref_node, nodes[node]));
      }
      if (!has_local) {
        RSS(REF_FAILURE, "cell with all ghost nodes");
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_validation_cell_volume(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_DBL volume;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

  ref_cell = ref_grid_tet(ref_grid);
  if (ref_grid_twod(ref_grid)) ref_cell = ref_grid_tri(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (ref_grid_twod(ref_grid)) {
      RSS(ref_node_tri_area(ref_node, nodes, &volume), "area");
    } else {
      RSS(ref_node_tet_vol(ref_node, nodes, &volume), "vol");
    }
    RAB(volume > 0.0, "negative volume tet", {
      REF_INT cell_node;
      printf("cell %d volume %e\n", cell, volume);
      each_ref_cell_cell_node(ref_cell, cell_node)
          ref_node_location(ref_node, nodes[cell_node]);
    });
  }

  return REF_SUCCESS;
}

REF_STATUS ref_validation_cell_volume_at_node(REF_GRID ref_grid, REF_INT node) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_DBL volume;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];

  ref_cell = ref_grid_tet(ref_grid);
  if (ref_grid_twod(ref_grid)) ref_cell = ref_grid_tri(ref_grid);

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
    if (ref_grid_twod(ref_grid)) {
      RSS(ref_node_tri_area(ref_node, nodes, &volume), "area");
    } else {
      RSS(ref_node_tet_vol(ref_node, nodes, &volume), "vol");
    }
    RAB(volume > 0.0, "negative volume tet", {
      REF_INT cell_node;
      printf("cell %d volume %e\n", cell, volume);
      each_ref_cell_cell_node(ref_cell, cell_node)
          ref_node_location(ref_node, nodes[cell_node]);
    });
  }

  return REF_SUCCESS;
}

REF_STATUS ref_validation_all(REF_GRID ref_grid) {
  RSS(ref_validation_unused_node(ref_grid), "unused node");
  RSS(ref_validation_boundary_face(ref_grid), "boundary face");
  RSS(ref_validation_cell_face(ref_grid), "cell face");
  RSS(ref_validation_cell_node(ref_grid), "cell node");
  RSS(ref_validation_cell_volume(ref_grid), "cell volume");

  return REF_SUCCESS;
}

REF_STATUS ref_validation_volume_status(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL min_volume, max_volume;
  REF_DBL volume;

  min_volume = REF_DBL_MAX;
  max_volume = REF_DBL_MIN;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_tet_vol(ref_grid_node(ref_grid), nodes, &volume), "vol");
    min_volume = MIN(min_volume, volume);
    max_volume = MAX(max_volume, volume);
  }
  volume = min_volume;
  RSS(ref_mpi_min(ref_mpi, &volume, &min_volume, REF_DBL_TYPE), "mpi min");
  RSS(ref_mpi_bcast(ref_mpi, &min_volume, 1, REF_DBL_TYPE), "min");
  volume = max_volume;
  RSS(ref_mpi_max(ref_mpi, &volume, &max_volume, REF_DBL_TYPE), "mpi max");
  RSS(ref_mpi_bcast(ref_mpi, &max_volume, 1, REF_DBL_TYPE), "min");

  if (ref_grid_once(ref_grid))
    printf("volume %.5e %.5e\n", min_volume, max_volume);

  return REF_SUCCESS;
}

REF_STATUS ref_validation_twod_orientation(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL correct_orientation, valid = REF_TRUE;

  if (!ref_grid_twod(ref_grid)) return REF_SUCCESS;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &correct_orientation),
        "valid");
    if (!correct_orientation) {
      valid = REF_FALSE;
      printf("tri %d %d %d %d\n", nodes[0], nodes[1], nodes[2], nodes[3]);
    }
  }

  RAS(valid, "incorrect twod tri orientation");

  return REF_SUCCESS;
}

REF_STATUS ref_validation_finite(REF_GRID ref_grid, REF_INT ldim,
                                 REF_DBL *field) {
  REF_INT i, node, n;
  n = 0;
  each_ref_node_valid_node(ref_grid_node(ref_grid), node) {
    for (i = 0; i < ldim; i++) {
      if (!isfinite(field[i + ldim * node])) n++;
    }
  }
  RSS(ref_mpi_allsum(ref_grid_mpi(ref_grid), &n, 1, REF_INT_TYPE), "global");
  if (n > 0) {
    if (ref_mpi_once(ref_grid_mpi(ref_grid))) printf("%d not finite\n", n);
    return REF_INVALID;
  }
  return REF_SUCCESS;
}

REF_STATUS ref_validation_repair(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_FACE ref_face;
  REF_CELL ref_cell;
  REF_INT *hits;
  REF_INT face;
  REF_INT group, cell, cell_face;
  REF_INT node;
  REF_INT nodes[4];
  REF_BOOL report;
  REF_STATUS code;

  RSS(ref_face_create(&ref_face, ref_grid), "face");

  ref_malloc(hits, ref_face_n(ref_face), REF_INT);

  for (face = 0; face < ref_face_n(ref_face); face++) hits[face] = 0;

  each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell) {
    each_ref_cell_valid_cell(ref_cell, cell) {
      each_ref_cell_cell_face(ref_cell, cell_face) {
        for (node = 0; node < 4; node++) {
          nodes[node] = ref_cell_f2n(ref_cell, node, cell_face, cell);
        }
        RSS(ref_face_with(ref_face, nodes, &face), "find cell face");
        hits[face]++;
      }
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    for (node = 0; node < 3; node++) {
      nodes[node] = ref_cell_c2n(ref_cell, node, cell);
    }
    nodes[3] = nodes[0];
    code = ref_face_with(ref_face, nodes, &face);
    if (REF_SUCCESS != code) {
      ref_node_location(ref_grid_node(ref_grid), nodes[0]);
      ref_node_location(ref_grid_node(ref_grid), nodes[1]);
      ref_node_location(ref_grid_node(ref_grid), nodes[2]);
      ref_node_location(ref_grid_node(ref_grid), nodes[3]);
    }
    RSS(code, "find tri");
    hits[face]++;
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    for (node = 0; node < 4; node++) {
      nodes[node] = ref_cell_c2n(ref_cell, node, cell);
    }
    RSS(ref_face_with(ref_face, nodes, &face), "find qua");
    hits[face]++;
  }

  for (face = 0; face < ref_face_n(ref_face); face++) {
    report = REF_FALSE;
    if (ref_mpi_para(ref_grid_mpi(ref_grid))) {
      report = report || (2 < hits[face]);
      if (ref_node_owned(ref_node, ref_face_f2n(ref_face, 0, face)) ||
          ref_node_owned(ref_node, ref_face_f2n(ref_face, 1, face)) ||
          ref_node_owned(ref_node, ref_face_f2n(ref_face, 2, face)) ||
          ref_node_owned(ref_node, ref_face_f2n(ref_face, 3, face))) {
        report = report || (2 > hits[face]);
      }
    } else {
      report = report || (2 != hits[face]);
    }
    if (report) {
      printf(" hits %d\n", hits[face]);
      for (node = 0; node < 4; node++) {
        nodes[node] = ref_face_f2n(ref_face, node, face);
      }
      printf("face %d nodes %d %d %d %d global " REF_GLOB_FMT " " REF_GLOB_FMT
             " " REF_GLOB_FMT " " REF_GLOB_FMT "\n",
             face, nodes[0], nodes[1], nodes[2], nodes[3],
             ref_node_global(ref_grid_node(ref_grid), nodes[0]),
             ref_node_global(ref_grid_node(ref_grid), nodes[1]),
             ref_node_global(ref_grid_node(ref_grid), nodes[2]),
             ref_node_global(ref_grid_node(ref_grid), nodes[3]));
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[0]), "n0");
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[1]), "n1");
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[2]), "n2");
      RSS(ref_node_location(ref_grid_node(ref_grid), nodes[3]), "n3");
      if (nodes[0] != nodes[3]) {
        REF_INT cell0, cell1;
        REF_INT hex_nodes[REF_CELL_MAX_SIZE_PER];
        REF_GLOB global;
        REF_INT new_node, cell_node;
        ref_cell = ref_grid_hex(ref_grid);
        printf("have quad\n");
        RSS(ref_cell_with_face(ref_cell, nodes, &cell0, &cell1), "hex");
        printf("hexes %d %d\n", cell0, cell1);
        RUS(REF_EMPTY, cell0, "cant find hex");
        REIS(REF_EMPTY, cell1, "found 2 hex");
        RSS(ref_cell_nodes(ref_cell, cell0, hex_nodes), "hex");
        RSS(ref_node_next_global(ref_node, &global), "next global");
        RSS(ref_node_add(ref_node, global, &new_node), "new node");
        ref_node_xyz(ref_node, 0, new_node) = 0.0;
        ref_node_xyz(ref_node, 1, new_node) = 0.0;
        ref_node_xyz(ref_node, 2, new_node) = 0.0;
        each_ref_cell_cell_node(ref_cell, cell_node) {
          ref_node_xyz(ref_node, 0, new_node) +=
              ref_node_xyz(ref_node, 0, hex_nodes[cell_node]);
          ref_node_xyz(ref_node, 1, new_node) +=
              ref_node_xyz(ref_node, 1, hex_nodes[cell_node]);
          ref_node_xyz(ref_node, 2, new_node) +=
              ref_node_xyz(ref_node, 2, hex_nodes[cell_node]);
        }
        ref_node_xyz(ref_node, 0, new_node) /=
            (REF_DBL)ref_cell_node_per(ref_cell);
        ref_node_xyz(ref_node, 1, new_node) /=
            (REF_DBL)ref_cell_node_per(ref_cell);
        ref_node_xyz(ref_node, 2, new_node) /=
            (REF_DBL)ref_cell_node_per(ref_cell);
        each_ref_cell_cell_face(ref_cell, cell_face) {
          REF_INT face_nodes[4];
          REF_BOOL split_face;
          REF_INT new_cell;

          face_nodes[0] = ref_cell_f2n(ref_cell, 0, cell_face, cell0);
          face_nodes[1] = ref_cell_f2n(ref_cell, 1, cell_face, cell0);
          face_nodes[2] = ref_cell_f2n(ref_cell, 2, cell_face, cell0);
          face_nodes[3] = ref_cell_f2n(ref_cell, 3, cell_face, cell0);
          RSS(ref_sort_same(4, nodes, face_nodes, &split_face), "same");
          if (split_face) {
            REF_INT tri_face;
            face_nodes[0] = ref_cell_f2n(ref_cell, 0, cell_face, cell0);
            face_nodes[1] = ref_cell_f2n(ref_cell, 1, cell_face, cell0);
            face_nodes[2] = ref_cell_f2n(ref_cell, 2, cell_face, cell0);
            face_nodes[3] = face_nodes[0];
            if (REF_SUCCESS == ref_face_with(ref_face, face_nodes, &tri_face)) {
              face_nodes[0] = ref_cell_f2n(ref_cell, 0, cell_face, cell0);
              face_nodes[1] = ref_cell_f2n(ref_cell, 1, cell_face, cell0);
              face_nodes[2] = ref_cell_f2n(ref_cell, 2, cell_face, cell0);
              face_nodes[3] = new_node;
              RSS(ref_cell_add(ref_grid_tet(ref_grid), face_nodes, &new_cell),
                  "tet");
              face_nodes[0] = ref_cell_f2n(ref_cell, 0, cell_face, cell0);
              face_nodes[1] = ref_cell_f2n(ref_cell, 2, cell_face, cell0);
              face_nodes[2] = ref_cell_f2n(ref_cell, 3, cell_face, cell0);
              face_nodes[3] = new_node;
              RSS(ref_cell_add(ref_grid_tet(ref_grid), face_nodes, &new_cell),
                  "tet");
            }
            face_nodes[0] = ref_cell_f2n(ref_cell, 0, cell_face, cell0);
            face_nodes[1] = ref_cell_f2n(ref_cell, 1, cell_face, cell0);
            face_nodes[2] = ref_cell_f2n(ref_cell, 3, cell_face, cell0);
            face_nodes[3] = face_nodes[0];
            if (REF_SUCCESS == ref_face_with(ref_face, face_nodes, &tri_face)) {
              face_nodes[0] = ref_cell_f2n(ref_cell, 0, cell_face, cell0);
              face_nodes[1] = ref_cell_f2n(ref_cell, 1, cell_face, cell0);
              face_nodes[2] = ref_cell_f2n(ref_cell, 3, cell_face, cell0);
              face_nodes[3] = new_node;
              RSS(ref_cell_add(ref_grid_tet(ref_grid), face_nodes, &new_cell),
                  "tet");
              face_nodes[0] = ref_cell_f2n(ref_cell, 1, cell_face, cell0);
              face_nodes[1] = ref_cell_f2n(ref_cell, 2, cell_face, cell0);
              face_nodes[2] = ref_cell_f2n(ref_cell, 3, cell_face, cell0);
              face_nodes[3] = new_node;
              RSS(ref_cell_add(ref_grid_tet(ref_grid), face_nodes, &new_cell),
                  "tet");
            }
          } else {
            REF_INT pyr_nodes[REF_CELL_MAX_SIZE_PER];
            pyr_nodes[0] = ref_cell_f2n(ref_cell, 0, cell_face, cell0);
            pyr_nodes[3] = ref_cell_f2n(ref_cell, 1, cell_face, cell0);
            pyr_nodes[4] = ref_cell_f2n(ref_cell, 2, cell_face, cell0);
            pyr_nodes[1] = ref_cell_f2n(ref_cell, 3, cell_face, cell0);
            pyr_nodes[2] = new_node;
            RSS(ref_cell_add(ref_grid_pyr(ref_grid), pyr_nodes, &new_cell),
                "pyr");
          }
        }
        RSS(ref_cell_remove(ref_grid_hex(ref_grid), cell0), "hex");
      }
    }
  }

  free(hits);

  RSS(ref_face_free(ref_face), "face free");

  return REF_SUCCESS;
}
