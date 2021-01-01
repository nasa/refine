
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

#include "ref_migrate.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_sort.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  if (!ref_mpi_para(ref_mpi)) { /* keep local, lose ghost */
    REF_GRID ref_grid;
    REF_MIGRATE ref_migrate;
    REF_INT keep, lose;
    REF_INT update_local, update_part;
    REF_ADJ ref_adj;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "set up grid");
    /* fake 2 proc */
    ref_node_part(ref_grid_node(ref_grid), 3) = 1;
    ref_node_part(ref_grid_node(ref_grid), 4) = 1;
    ref_node_part(ref_grid_node(ref_grid), 5) = 1;
    RSS(ref_migrate_create(&ref_migrate, ref_grid), "set up mig");

    keep = 0;
    lose = 3;
    RSS(ref_migrate_2d_agglomeration_keep(ref_migrate, keep, lose), "0-3");

    REIS(0, ref_migrate_global(ref_migrate, 0), "mark");
    REIS(REF_EMPTY, ref_migrate_global(ref_migrate, 3), "mark");

    ref_adj = ref_migrate_parent_local(ref_migrate);
    update_local = ref_adj_item_ref(ref_adj, ref_adj_first(ref_adj, keep));
    REIS(3, update_local, "glob");

    ref_adj = ref_migrate_parent_part(ref_migrate);
    update_part = ref_adj_item_ref(ref_adj, ref_adj_first(ref_adj, keep));
    REIS(1, update_part, "part");

    RSS(ref_migrate_free(ref_migrate), "free migrate");
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (!ref_mpi_para(ref_mpi)) { /* keep ghost, lose local */
    REF_GRID ref_grid;
    REF_MIGRATE ref_migrate;
    REF_INT keep, lose;
    REF_ADJ ref_adj;

    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "set up grid");
    /* fake 2 proc */
    ref_node_part(ref_grid_node(ref_grid), 0) = 1;
    ref_node_part(ref_grid_node(ref_grid), 1) = 1;
    ref_node_part(ref_grid_node(ref_grid), 2) = 1;
    RSS(ref_migrate_create(&ref_migrate, ref_grid), "set up mig");

    keep = 0;
    lose = 3;
    RSS(ref_migrate_2d_agglomeration_keep(ref_migrate, keep, lose), "0-3");

    REIS(REF_EMPTY, ref_migrate_global(ref_migrate, 0), "mark");
    REIS(REF_EMPTY, ref_migrate_global(ref_migrate, 3), "mark");

    ref_adj = ref_migrate_parent_local(ref_migrate);
    RAS(ref_adj_empty(ref_adj, keep), "glob");

    ref_adj = ref_migrate_parent_part(ref_migrate);
    RAS(ref_adj_empty(ref_adj, keep), "part");

    RSS(ref_migrate_free(ref_migrate), "free migrate");
    RSS(ref_grid_free(ref_grid), "free grid");
  }

  if (1 == argc) { /* part and migrate tet b8.ugrid, world comm */
    REF_GRID import_grid;
    char grid_file[] = "ref_migrate_test.b8.ugrid";

    if (ref_mpi_once(ref_mpi)) {
      REF_GRID export_grid;
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }

    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");
    RSS(ref_migrate_to_balance(import_grid), "create");

    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  if (1 == argc) {        /* part and migrate tet b8.ugrid, world comm */
    REF_UINT x = 1048576; /* (2^21 >> 1) */
    REF_UINT y = 0;
    REF_UINT z = 0;

    REF_ULONG morton_id = ref_migrate_morton_id(x, y, z);
    REF_ULONG expected = (REF_ULONG)1 << 60;
    RES(expected, morton_id, "Expected morton id mismatch.");
  }

  if (1 == argc) { /* part and migrate tet lb8.ugrid, split comm */
    REF_MPI split_mpi;
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_migrate_test.lb8.ugrid";
    char tec_file[1024];

    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
    }
    RSS(ref_grid_free(export_grid), "free");

    RSS(ref_mpi_half_comm(ref_mpi, &split_mpi), "split");
    RSS(ref_part_by_extension(&import_grid, split_mpi, grid_file), "import");

    if (ref_mpi_para(split_mpi)) {
      snprintf(tec_file, 1024, "ref_migrate_part_%d.tec",
               ref_mpi_rank(ref_mpi));
      RSS(ref_gather_tec_part(import_grid, tec_file), "tec part");
    }

    RSS(ref_migrate_to_balance(import_grid), "migrate");

    if (ref_mpi_para(split_mpi)) {
      snprintf(tec_file, 1024, "ref_migrate_bal_%d.tec", ref_mpi_rank(ref_mpi));
      RSS(ref_gather_tec_part(import_grid, tec_file), "tec part");
    }

    {
      REF_LONG *global;
      REF_INT group;
      REF_CELL ref_cell;
      each_ref_grid_all_ref_cell(import_grid, group, ref_cell) {
        RSS(ref_cell_global(ref_cell, ref_grid_node(import_grid), &global),
            "cell global");
        ref_free(global);
      }
    }

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_mpi_join_comm(split_mpi), "join");
    RSS(ref_mpi_free(split_mpi), "free");

    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  if (1 == argc) { /* part and migrate tet b8.ugrid native rcb */
    REF_GRID import_grid;
    char grid_file[] = "ref_migrate_test.b8.ugrid";

    if (ref_mpi_once(ref_mpi)) {
      REF_GRID export_grid;
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }

    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");
    ref_grid_partitioner(import_grid) = REF_MIGRATE_NATIVE_RCB;
    RSS(ref_migrate_to_balance(import_grid), "create");

    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  if (1 == argc) { /* replicate tet b8.ugrid */
    REF_GRID ref_grid = NULL;

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "set up tet");
    }

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_grid_free(ref_grid), "free");
    }
  }

  if (1 < argc) { /* part and migrate argument, world comm */
    REF_GRID import_grid;

    if (ref_mpi_once(ref_mpi))
      printf("%d procs, read %s\n", ref_mpi_n(ref_mpi), argv[1]);

    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_part_by_extension(&import_grid, ref_mpi, argv[1]), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "read/part grid");

    RSS(ref_migrate_to_balance(import_grid), "new part");
    ref_mpi_stopwatch_stop(ref_mpi, "balance parts");

    RSS(ref_gather_tec_part(import_grid, "ref_migrate_test_world.tec"),
        "part_viz");
    ref_mpi_stopwatch_stop(ref_mpi, "gather ref_migrate_test_world.tec");

    RSS(ref_grid_free(import_grid), "free");
  }

  if (1 < argc) { /* part and migrate argument, world comm, native */
    REF_GRID import_grid;

    if (ref_mpi_once(ref_mpi))
      printf("%d procs, read %s\n", ref_mpi_n(ref_mpi), argv[1]);

    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_part_by_extension(&import_grid, ref_mpi, argv[1]), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "read/part grid");

    ref_grid_partitioner(import_grid) = REF_MIGRATE_NATIVE_RCB;
    RSS(ref_migrate_to_balance(import_grid), "new part");
    ref_mpi_stopwatch_stop(ref_mpi, "balance parts");

    RSS(ref_gather_tec_part(import_grid, "ref_migrate_test_native.tec"),
        "part_viz");
    ref_mpi_stopwatch_stop(ref_mpi, "gather ref_migrate_test_world.tec");

    RSS(ref_grid_free(import_grid), "free");
  }

  if (1 < argc) { /* part and migrate argument, split comm */
    REF_MPI split_mpi;
    REF_GRID import_grid;
    char tec_file[1024];

    RSS(ref_mpi_half_comm(ref_mpi, &split_mpi), "split");
    if (ref_mpi_once(split_mpi))
      printf("%d procs, read %s\n", ref_mpi_n(split_mpi), argv[1]);

    ref_mpi_stopwatch_start(split_mpi);
    RSS(ref_part_by_extension(&import_grid, split_mpi, argv[1]), "import");
    ref_mpi_stopwatch_stop(split_mpi, "read");

    if (ref_mpi_para(split_mpi)) {
      snprintf(tec_file, 1024, "ref_migrate_test_%d.tec",
               ref_mpi_rank(ref_mpi));
      RSS(ref_gather_tec_part(import_grid, tec_file), "tec part");
      ref_mpi_stopwatch_stop(split_mpi, "tec");
    }

    RSS(ref_migrate_to_balance(import_grid), "migrate");
    ref_mpi_stopwatch_stop(split_mpi, "migrate");

    if (ref_mpi_para(split_mpi)) {
      snprintf(tec_file, 1024, "ref_migrate_test_bal_%d.tec",
               ref_mpi_rank(ref_mpi));
      RSS(ref_gather_tec_part(import_grid, tec_file), "tec part");
      ref_mpi_stopwatch_stop(split_mpi, "tec");
    }

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_mpi_join_comm(split_mpi), "join");
    RSS(ref_mpi_free(split_mpi), "free");
  }

  {
    REF_INT n = 2;
    REF_DBL xyz[] = {0, 0, 0, 5, 2, 3};
    REF_INT dir;
    RSS(ref_migrate_split_dir(ref_mpi, n, xyz, &dir), "dir");
    REIS(0, dir, "expects x");
  }

  {
    REF_INT n = 2;
    REF_DBL xyz[] = {0, 0, 0, 1, 7, 3};
    REF_INT dir;
    RSS(ref_migrate_split_dir(ref_mpi, n, xyz, &dir), "dir");
    REIS(1, dir, "expects y");
  }

  {
    REF_INT n = 2;
    REF_DBL xyz[] = {0, 0, 0, 1, 2, 3};
    REF_INT dir;
    RSS(ref_migrate_split_dir(ref_mpi, n, xyz, &dir), "dir");
    REIS(2, dir, "expects z");
  }

  {
    REF_INT n = 2;
    REF_DBL ratio;
    RSS(ref_migrate_split_ratio(n, &ratio), "ratio");
    RWDS(0.5, ratio, -1.0, "equal");
  }

  {
    REF_INT n = 5;
    REF_DBL ratio;
    RSS(ref_migrate_split_ratio(n, &ratio), "ratio");
    RWDS(0.4, ratio, -1.0, "equal");
  }

  {
    REF_INT n = 0;
    REF_DBL ratio;
    REIS(REF_DIV_ZERO, ref_migrate_split_ratio(n, &ratio), "ratio");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
