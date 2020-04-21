
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

#include "ref_gather.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_args.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_matrix.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_sort.h"

static REF_STATUS ref_gather_meshb_fixture(REF_MPI ref_mpi,
                                           const char *filename,
                                           REF_INT version) {
  REF_GRID ref_grid;
  REF_GEOM ref_geom;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT type, id, node;
  REF_DBL param[2];

  if (ref_mpi_once(ref_mpi)) {
    RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(ref_grid) = version;
    ref_geom = ref_grid_geom(ref_grid);
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 15;
    RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &cell), "add edge");
    type = REF_GEOM_NODE;
    id = 1;
    node = 0;
    RSS(ref_geom_add(ref_geom, node, type, id, param), "add geom node");
    id = 2;
    node = 1;
    RSS(ref_geom_add(ref_geom, node, type, id, param), "add geom node");
    type = REF_GEOM_EDGE;
    id = 15;
    node = 0;
    param[0] = 10.0;
    RSS(ref_geom_add(ref_geom, node, type, id, param), "add geom edge");
    id = 15;
    node = 1;
    param[0] = 20.0;
    RSS(ref_geom_add(ref_geom, node, type, id, param), "add geom edge");
    type = REF_GEOM_FACE;
    id = 3;
    node = 0;
    param[0] = 10.0;
    param[1] = 20.0;
    RSS(ref_geom_add(ref_geom, node, type, id, param), "add geom face");
    type = REF_GEOM_FACE;
    id = 3;
    node = 1;
    param[0] = 11.0;
    param[1] = 20.0;
    RSS(ref_geom_add(ref_geom, node, type, id, param), "add geom face");
    type = REF_GEOM_FACE;
    id = 3;
    node = 2;
    param[0] = 10.5;
    param[1] = 21.0;
    RSS(ref_geom_add(ref_geom, node, type, id, param), "add geom face");

    RSS(ref_export_by_extension(ref_grid, filename), "export");
    RSS(ref_grid_free(ref_grid), "free");
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_gather_bbox_intersects(REF_DBL *bbox1, REF_DBL *bbox2,
                                             REF_BOOL *intersects) {
  REF_INT i;

  *intersects = REF_FALSE;

  for (i = 0; i < 3; i++) {
    if (bbox1[i + 3] < bbox2[i]) return REF_SUCCESS;
    if (bbox1[i] > bbox2[i + 3]) return REF_SUCCESS;
  }

  *intersects = REF_TRUE;

  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_INT pos;
  REF_MPI ref_mpi;
  REF_BOOL transmesh = REF_FALSE;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  RXS(ref_args_find(argc, argv, "--transmesh", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) transmesh = REF_TRUE;

  RXS(ref_args_find(argc, argv, "--subset", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos == 1 && argc == 12) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CELL ref_cell;
    char *filename;
    REF_INT i, node, ldim, group, cell, cell_node;
    REF_INT nodes[REF_CELL_MAX_SIZE_PER];
    REF_DBL bbox[6], cell_bbox[6];
    REF_DBL *solution;
    REF_BOOL keep, intersects;
    REF_INT *keep_node;

    filename = argv[2];
    if (ref_mpi_once(ref_mpi)) printf("part whole grid %s\n", filename);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, filename), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "part whole grid");
    ref_node = ref_grid_node(ref_grid);

    filename = argv[3];
    if (ref_mpi_once(ref_mpi)) printf("part whole solution %s\n", filename);
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &solution, filename),
        "part whole solution");
    ref_mpi_stopwatch_stop(ref_mpi, "part whole solution");

    for (i = 0; i < 6; i++) {
      bbox[i] = atof(argv[i + 4]);
    }
    if (ref_mpi_once(ref_mpi)) {
      printf("bounding box min %f %f %f\n", bbox[0], bbox[1], bbox[2]);
      printf("bounding dx dy dz %f %f %f\n", bbox[3], bbox[4], bbox[5]);
    }
    for (i = 0; i < 3; i++) {
      bbox[i + 3] += bbox[i];
    }
    if (ref_mpi_once(ref_mpi)) {
      printf("bounding box min %f %f %f\n", bbox[0], bbox[1], bbox[2]);
      printf("bounding box max %f %f %f\n", bbox[3], bbox[4], bbox[5]);
    }

    ref_malloc_init(keep_node, ref_node_max(ref_node), REF_INT, 0);

    each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        for (i = 0; i < 3; i++) {
          cell_bbox[i] = ref_node_xyz(ref_node, i, nodes[0]);
          cell_bbox[i + 3] = ref_node_xyz(ref_node, i, nodes[0]);
        }
        each_ref_cell_cell_node(ref_cell, cell_node) {
          node = nodes[cell_node];
          for (i = 0; i < 3; i++) {
            cell_bbox[i] =
                MIN(cell_bbox[i], ref_node_xyz(ref_node, i, nodes[cell_node]));
            cell_bbox[i + 3] = MAX(cell_bbox[i + 3],
                                   ref_node_xyz(ref_node, i, nodes[cell_node]));
          }
        }
        RSS(ref_gather_bbox_intersects(bbox, cell_bbox, &intersects),
            "compute intersection");
        if (intersects) {
          each_ref_cell_cell_node(ref_cell, cell_node) {
            node = nodes[cell_node];
            keep_node[node] = REF_TRUE;
          }
        }
      }
    }

    RSS(ref_node_ghost_int(ref_node, keep_node, 1), "ghost");
    ref_mpi_stopwatch_stop(ref_mpi, "mark kept nodes");

    each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        keep = REF_TRUE;
        each_ref_cell_cell_node(ref_cell, cell_node) {
          node = nodes[cell_node];
          keep = keep && keep_node[node];
        }
        if (!keep) RSS(ref_cell_remove(ref_cell, cell), "rm");
      }
    }
    ref_free(keep_node);
    ref_mpi_stopwatch_stop(ref_mpi, "prune cells");

    each_ref_node_valid_node(ref_node, node) {
      keep = REF_FALSE;
      each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
        if (!ref_cell_node_empty(ref_cell, node)) keep = REF_TRUE;
      }
      if (!keep) {
        if (ref_node_owned(ref_node, node)) {
          RSS(ref_node_remove(ref_node, node), "rm");
        } else {
          RSS(ref_node_remove_without_global(ref_node, node), "rm");
        }
        RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), node), "rm");
      }
    }
    ref_mpi_stopwatch_stop(ref_mpi, "prune nodes");

    RSS(ref_node_synchronize_globals(ref_node), "sync");
    ref_mpi_stopwatch_stop(ref_mpi, "sync node globals");
    if (ref_mpi_once(ref_mpi))
      printf("global " REF_GLOB_FMT "\n", ref_node_n_global(ref_node));

    filename = argv[10];
    if (ref_mpi_once(ref_mpi)) printf("gather mesh subset %s\n", filename);
    RSS(ref_gather_by_extension(ref_grid, filename), "gather");
    ref_mpi_stopwatch_stop(ref_mpi, "gather mesh subset");

    filename = argv[11];
    if (ref_mpi_once(ref_mpi)) printf("gather solution subset %s\n", filename);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, solution, NULL,
                                       filename),
        "gather");
    ref_mpi_stopwatch_stop(ref_mpi, "gather solution subset");

    ref_free(solution);
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (2 == argc && !transmesh) {
    REF_GRID import_grid;

    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_part_by_extension(&import_grid, ref_mpi, argv[1]), "import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "read");
    RSS(ref_migrate_to_balance(import_grid), "balance");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "balance");

    ref_mpi_stopwatch_start(ref_grid_mpi(import_grid));
    RSS(ref_gather_by_extension(import_grid, "ref_gather_test.meshb"),
        "gather");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "meshb");

    ref_mpi_stopwatch_start(ref_grid_mpi(import_grid));
    RSS(ref_gather_by_extension(import_grid, "ref_gather_test.b8.ugrid"),
        "gather");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "b8.ugrid");

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");

    return 0;
  }

  if (4 == argc) {
    REF_GRID import_grid;
    REF_INT ldim;
    REF_DBL *field;

    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_part_by_extension(&import_grid, ref_mpi, argv[1]), "import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "read grid");
    RSS(ref_part_scalar(ref_grid_node(import_grid), &ldim, &field, argv[2]),
        "field");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "read field");
    RSS(ref_gather_scalar_by_extension(import_grid, ldim, field, NULL, argv[3]),
        "field");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "write tec");

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");

    return 0;
  }

  if (3 == argc) {
    REF_GRID import_grid;

    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_part_by_extension(&import_grid, ref_mpi, argv[1]), "import");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "read grid");
    RSS(ref_part_metric(ref_grid_node(import_grid), argv[2]), "field");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "read metric");
    RSS(ref_gather_tec_movie_record_button(ref_grid_gather(import_grid),
                                           REF_TRUE),
        "movie on");
    ref_gather_low_quality_zone(ref_grid_gather(import_grid)) = REF_TRUE;
    RSS(ref_gather_tec_movie_frame(import_grid, "triage"), "movie frame");
    ref_mpi_stopwatch_stop(ref_grid_mpi(import_grid), "movie frame");

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");

    return 0;
  }

  { /* gather prism lb8 */
    REF_GRID ref_grid;

    RSS(ref_fixture_pri_stack_grid(&ref_grid, ref_mpi), "set up tet");

    RSS(ref_gather_by_extension(ref_grid, "ref_gather_test.lb8.ugrid"),
        "gather");

    RSS(ref_grid_free(ref_grid), "free");
    if (ref_mpi_once(ref_mpi))
      REIS(0, remove("ref_gather_test.lb8.ugrid"), "test clean up");
  }

  { /* gather prism b8 */
    REF_GRID ref_grid;

    RSS(ref_fixture_pri_stack_grid(&ref_grid, ref_mpi), "set up tet");

    RSS(ref_gather_by_extension(ref_grid, "ref_gather_test.b8.ugrid"),
        "gather");

    RSS(ref_grid_free(ref_grid), "free");
    if (ref_mpi_once(ref_mpi))
      REIS(0, remove("ref_gather_test.b8.ugrid"), "test clean up");
  }

  { /* recycle tet brick lb8.ugrid */
    REF_GRID seq_grid, para_grid;
    char seq_file[] = "ref_gather_test_seq.lb8.ugrid";
    char para_file[] = "ref_gather_test_para.lb8.ugrid";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&seq_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(seq_grid, seq_file), "export");
    }
    RSS(ref_part_by_extension(&para_grid, ref_mpi, seq_file), "part");
    RSS(ref_gather_by_extension(para_grid, para_file), "gather");
    RSS(ref_grid_free(para_grid), "free");
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_import_by_extension(&para_grid, ref_mpi, para_file), "import");
      REIS(ref_node_n(ref_grid_node(seq_grid)),
           ref_node_n(ref_grid_node(para_grid)), "nnode mis-match");
      REIS(ref_cell_n(ref_grid_tri(seq_grid)),
           ref_cell_n(ref_grid_tri(para_grid)), "ntri mis-match");
      REIS(ref_cell_n(ref_grid_tet(seq_grid)),
           ref_cell_n(ref_grid_tet(para_grid)), "ntet mis-match");
      RSS(ref_grid_free(para_grid), "free");
      RSS(ref_grid_free(seq_grid), "free");
      REIS(0, remove(seq_file), "test clean up");
      REIS(0, remove(para_file), "test clean up");
    }
  }

  { /* recycle tet brick b8.ugrid */
    REF_GRID seq_grid, para_grid;
    char seq_file[] = "ref_gather_test_seq.b8.ugrid";
    char para_file[] = "ref_gather_test_para.b8.ugrid";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&seq_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(seq_grid, seq_file), "export");
    }
    RSS(ref_part_by_extension(&para_grid, ref_mpi, seq_file), "part");
    RSS(ref_gather_by_extension(para_grid, para_file), "gather");
    RSS(ref_grid_free(para_grid), "free");
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_import_by_extension(&para_grid, ref_mpi, para_file), "import");
      REIS(ref_node_n(ref_grid_node(seq_grid)),
           ref_node_n(ref_grid_node(para_grid)), "nnode mis-match");
      REIS(ref_cell_n(ref_grid_tri(seq_grid)),
           ref_cell_n(ref_grid_tri(para_grid)), "ntri mis-match");
      REIS(ref_cell_n(ref_grid_tet(seq_grid)),
           ref_cell_n(ref_grid_tet(para_grid)), "ntet mis-match");
      RSS(ref_grid_free(para_grid), "free");
      RSS(ref_grid_free(seq_grid), "free");
      REIS(0, remove(seq_file), "test clean up");
      REIS(0, remove(para_file), "test clean up");
    }
  }

  { /* recycle tet brick lb8l.ugrid */
    REF_GRID seq_grid, para_grid;
    char seq_file[] = "ref_gather_test_seq.lb8l.ugrid";
    char para_file[] = "ref_gather_test_para.lb8l.ugrid";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&seq_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(seq_grid, seq_file), "export");
    }
    RSS(ref_part_by_extension(&para_grid, ref_mpi, seq_file), "part");
    RSS(ref_gather_by_extension(para_grid, para_file), "gather");
    RSS(ref_grid_free(para_grid), "free");
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_import_by_extension(&para_grid, ref_mpi, para_file), "import");
      REIS(ref_node_n(ref_grid_node(seq_grid)),
           ref_node_n(ref_grid_node(para_grid)), "nnode mis-match");
      REIS(ref_cell_n(ref_grid_tri(seq_grid)),
           ref_cell_n(ref_grid_tri(para_grid)), "ntri mis-match");
      REIS(ref_cell_n(ref_grid_tet(seq_grid)),
           ref_cell_n(ref_grid_tet(para_grid)), "ntet mis-match");
      RSS(ref_grid_free(para_grid), "free");
      RSS(ref_grid_free(seq_grid), "free");
      REIS(0, remove(seq_file), "test clean up");
      REIS(0, remove(para_file), "test clean up");
    }
  }

  { /* recycle tet brick b8l.ugrid */
    REF_GRID seq_grid, para_grid;
    char seq_file[] = "ref_gather_test_seq.b8l.ugrid";
    char para_file[] = "ref_gather_test_para.b8l.ugrid";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&seq_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(seq_grid, seq_file), "export");
    }
    RSS(ref_part_by_extension(&para_grid, ref_mpi, seq_file), "part");
    RSS(ref_gather_by_extension(para_grid, para_file), "gather");
    RSS(ref_grid_free(para_grid), "free");
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_import_by_extension(&para_grid, ref_mpi, para_file), "import");
      REIS(ref_node_n(ref_grid_node(seq_grid)),
           ref_node_n(ref_grid_node(para_grid)), "nnode mis-match");
      REIS(ref_cell_n(ref_grid_tri(seq_grid)),
           ref_cell_n(ref_grid_tri(para_grid)), "ntri mis-match");
      REIS(ref_cell_n(ref_grid_tet(seq_grid)),
           ref_cell_n(ref_grid_tet(para_grid)), "ntet mis-match");
      RSS(ref_grid_free(para_grid), "free");
      RSS(ref_grid_free(seq_grid), "free");
      REIS(0, remove(seq_file), "test clean up");
      REIS(0, remove(para_file), "test clean up");
    }
  }

  { /* export part gather .meshb tet with cad_model, version 2 */
    REF_GRID ref_grid;
    REF_INT version = 2;
    char file[] = "ref_gather_test_ver2.meshb";

    RSS(ref_gather_meshb_fixture(ref_mpi, file, version), "fixture");
    if (transmesh && ref_mpi_once(ref_mpi))
      REIS(0,
           system("transmesh ref_gather_test_ver2.meshb "
                  "ref_gather_test_ver2_import.mesh"),
           "mesh");

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, file), "gather");
    RSS(ref_gather_by_extension(ref_grid, file), "gather");
    RSS(ref_grid_free(ref_grid), "free");

    if (transmesh && ref_mpi_once(ref_mpi))
      REIS(0,
           system("transmesh ref_gather_test_ver2.meshb "
                  "ref_gather_test_ver2_gather.mesh"),
           "mesh");

    if (ref_mpi_once(ref_mpi)) {
      if (!transmesh) REIS(0, remove(file), "test clean up");
    }
  }

  { /* export part gather .meshb tet with cad_model, version 3 */
    REF_GRID ref_grid;
    REF_INT version = 3;
    char file[] = "ref_gather_test_ver3.meshb";

    RSS(ref_gather_meshb_fixture(ref_mpi, file, version), "fixture");
    if (transmesh && ref_mpi_once(ref_mpi))
      REIS(0,
           system("transmesh ref_gather_test_ver3.meshb "
                  "ref_gather_test_ver3_import.mesh"),
           "mesh");

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, file), "gather");
    RSS(ref_gather_by_extension(ref_grid, file), "gather");
    RSS(ref_grid_free(ref_grid), "free");

    if (transmesh && ref_mpi_once(ref_mpi))
      REIS(0,
           system("transmesh ref_gather_test_ver3.meshb "
                  "ref_gather_test_ver3_gather.mesh"),
           "mesh");

    if (ref_mpi_once(ref_mpi)) {
      if (!transmesh) REIS(0, remove(file), "test clean up");
    }
  }

  { /* export part gather .meshb tet with cad_model, version 4 */
    REF_GRID ref_grid;
    REF_INT version = 4;
    char file[] = "ref_gather_test_ver4.meshb";

    RSS(ref_gather_meshb_fixture(ref_mpi, file, version), "fixture");
    if (transmesh && ref_mpi_once(ref_mpi))
      REIS(0,
           system("transmesh ref_gather_test_ver4.meshb "
                  "ref_gather_test_ver4_import.mesh"),
           "mesh");

    RSS(ref_part_by_extension(&ref_grid, ref_mpi, file), "gather");
    RSS(ref_gather_by_extension(ref_grid, file), "gather");
    RSS(ref_grid_free(ref_grid), "free");

    if (transmesh && ref_mpi_once(ref_mpi))
      REIS(0,
           system("transmesh ref_gather_test_ver4.meshb "
                  "ref_gather_test_ver4_gather.mesh"),
           "mesh");

    if (ref_mpi_once(ref_mpi)) {
      if (!transmesh) REIS(0, remove(file), "test clean up");
    }
  }

  { /* gather .solb by extension */
    REF_GRID ref_grid;
    REF_INT ldim;
    REF_DBL *scalar;
    const char **scalar_names = NULL;
    char filename[] = "ref_gather_test.solb";

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    ldim = 2;
    ref_malloc_init(scalar, ldim * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 1.0);
    scalar_names = NULL;

    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, scalar, scalar_names,
                                       filename),
        "gather");

    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "free");

    if (ref_mpi_once(ref_mpi)) REIS(0, remove(filename), "test clean up");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
