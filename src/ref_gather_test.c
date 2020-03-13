
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

int main(int argc, char *argv[]) {
  REF_INT pos;
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  RXS(ref_args_find(argc, argv, "--subset", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos == 1 && argc == 12) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CELL ref_cell;
    char *filename;
    REF_INT i, node, ldim, group, cell, cell_node;
    REF_INT nodes[REF_CELL_MAX_SIZE_PER];
    REF_DBL bbox[6];
    REF_DBL *solution;
    REF_BOOL keep;

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
    each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        keep = REF_FALSE;
        each_ref_cell_cell_node(ref_cell, cell_node) {
          node = nodes[cell_node];
          if (bbox[0] <= ref_node_xyz(ref_node, 0, node) &&
              ref_node_xyz(ref_node, 0, node) <= bbox[3] &&
              bbox[1] <= ref_node_xyz(ref_node, 1, node) &&
              ref_node_xyz(ref_node, 1, node) <= bbox[4] &&
              bbox[2] <= ref_node_xyz(ref_node, 2, node) &&
              ref_node_xyz(ref_node, 2, node) <= bbox[5]) {
            keep = REF_TRUE;
          }
        }
        if (!keep) RSS(ref_cell_remove(ref_cell, cell), "rm");
      }
    }
    ref_mpi_stopwatch_stop(ref_mpi, "prune cells");

    each_ref_node_valid_node(ref_node, node) {
      keep = REF_FALSE;
      each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
        if (!ref_cell_node_empty(ref_cell, node)) keep = REF_TRUE;
      }
      if (!keep) RSS(ref_node_remove(ref_node, node), "rm");
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

  if (2 == argc) {
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
    RSS(ref_gather_scalar_tec(import_grid, ldim, field, NULL, argv[3]),
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

  { /* export import .meshb tet with cad_model */
    REF_GRID export_grid, import_grid;
    REF_GEOM ref_geom;
    char file[] = "ref_gather_test.meshb";

    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    ref_geom = ref_grid_geom(export_grid);
    ref_geom_cad_data_size(ref_geom) = 3;
    ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                      ref_geom_cad_data_size(ref_geom), REF_BYTE);
    ref_geom_cad_data(ref_geom)[0] = 5;
    ref_geom_cad_data(ref_geom)[1] = 4;
    ref_geom_cad_data(ref_geom)[2] = 3;
    RSS(ref_gather_by_extension(export_grid, file), "gather");
    RSS(ref_grid_free(export_grid), "free");

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
      ref_geom = ref_grid_geom(import_grid);
      REIS(3, ref_geom_cad_data_size(ref_geom), "cad size");
      REIS(5, ref_geom_cad_data(ref_geom)[0], "cad[0]");
      REIS(4, ref_geom_cad_data(ref_geom)[1], "cad[1]");
      REIS(3, ref_geom_cad_data(ref_geom)[2], "cad[2]");
      RSS(ref_grid_free(import_grid), "free");
      REIS(0, remove(file), "test clean up");
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
