
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

#include "ref_part.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adapt.h"
#include "ref_adj.h"
#include "ref_cavity.h"
#include "ref_cell.h"
#include "ref_clump.h"
#include "ref_collapse.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_subdiv.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT ngeom;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  if (1 < argc) { /* part */
    REF_GRID import_grid;
    char viz_file[256];

    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_part_by_extension(&import_grid, ref_mpi, argv[1]), "import");
    ref_mpi_stopwatch_stop(ref_mpi, "import");

    sprintf(viz_file, "ref_part_test_n%d_p%d.tec",
            ref_mpi_n(ref_grid_mpi(import_grid)),
            ref_mpi_rank(ref_grid_mpi(import_grid)));
    RSS(ref_export_by_extension(import_grid, viz_file), "export");
    ref_mpi_stopwatch_stop(ref_mpi, "export");

    RSS(ref_gather_tec_part(import_grid, "ref_part_test.tec"), "part_viz");
    ref_mpi_stopwatch_stop(ref_mpi, "gather");

    RSS(ref_gather_ngeom(ref_grid_node(import_grid), ref_grid_geom(import_grid),
                         REF_GEOM_FACE, &ngeom),
        "count ngeom");
    if (ngeom > 0) {
      RSS(ref_geom_tec_para_shard(import_grid, "ref_part_geom"), "shard geom");
      ref_mpi_stopwatch_stop(ref_mpi, "geom");
    }

    RSS(ref_grid_free(import_grid), "free");

    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* one thing per part, equal things and parts */
    REIS(1, ref_part_large_part_size(4, 4), "large");
    REIS(0, ref_part_small_part_size(4, 4), "large");
    REIS(4, ref_part_n_large_part(4, 4), "large");

    REIS(0, ref_part_first(4, 4, 0), "first");
    REIS(1, ref_part_first(4, 4, 1), "first");
    REIS(2, ref_part_first(4, 4, 2), "first");
    REIS(3, ref_part_first(4, 4, 3), "first");
    REIS(4, ref_part_first(4, 4, 4), "first");

    REIS(0, ref_part_implicit(4, 4, 0), "part");
    REIS(1, ref_part_implicit(4, 4, 1), "part");
    REIS(2, ref_part_implicit(4, 4, 2), "part");
    REIS(3, ref_part_implicit(4, 4, 3), "part");
  }

  { /* one thing per part, less things than parts */
    REIS(1, ref_part_large_part_size(2, 4), "large");
    REIS(0, ref_part_small_part_size(2, 4), "large");
    REIS(2, ref_part_n_large_part(2, 4), "large");

    REIS(0, ref_part_first(2, 4, 0), "first");
    REIS(1, ref_part_first(2, 4, 1), "first");
    REIS(2, ref_part_first(2, 4, 2), "first");
    REIS(2, ref_part_first(2, 4, 3), "first");
    REIS(2, ref_part_first(2, 4, 4), "first");

    REIS(0, ref_part_implicit(2, 4, 0), "part");
    REIS(1, ref_part_implicit(2, 4, 1), "part");
  }

  { /* two things per part, equal things and parts */
    REIS(2, ref_part_large_part_size(4, 2), "large");
    REIS(1, ref_part_small_part_size(4, 2), "large");
    REIS(2, ref_part_n_large_part(4, 2), "large");

    REIS(0, ref_part_first(4, 2, 0), "first");
    REIS(2, ref_part_first(4, 2, 1), "first");
    REIS(4, ref_part_first(4, 2, 2), "first");

    REIS(0, ref_part_implicit(4, 2, 0), "part");
    REIS(0, ref_part_implicit(4, 2, 1), "part");
    REIS(1, ref_part_implicit(4, 2, 2), "part");
    REIS(1, ref_part_implicit(4, 2, 3), "part");
  }

  { /* two things per part, less things than parts */
    REIS(2, ref_part_large_part_size(7, 4), "large");
    REIS(1, ref_part_small_part_size(7, 4), "large");
    REIS(3, ref_part_n_large_part(7, 4), "large");

    REIS(0, ref_part_first(7, 4, 0), "first");
    REIS(2, ref_part_first(7, 4, 1), "first");
    REIS(4, ref_part_first(7, 4, 2), "first");
    REIS(6, ref_part_first(7, 4, 3), "first");
    REIS(7, ref_part_first(7, 4, 4), "first");

    REIS(0, ref_part_implicit(7, 4, 0), "part");
    REIS(0, ref_part_implicit(7, 4, 1), "part");
    REIS(1, ref_part_implicit(7, 4, 2), "part");
    REIS(1, ref_part_implicit(7, 4, 3), "part");
    REIS(2, ref_part_implicit(7, 4, 4), "part");
    REIS(2, ref_part_implicit(7, 4, 5), "part");
    REIS(3, ref_part_implicit(7, 4, 6), "part");
  }

  { /* one part */
    REIS(4, ref_part_large_part_size(4, 1), "large");
    REIS(3, ref_part_small_part_size(4, 1), "large");
    REIS(1, ref_part_n_large_part(4, 1), "large");

    REIS(0, ref_part_first(4, 1, 0), "first");
    REIS(4, ref_part_first(4, 1, 1), "first");

    REIS(0, ref_part_implicit(4, 1, 0), "part");
    REIS(0, ref_part_implicit(4, 1, 1), "part");
    REIS(0, ref_part_implicit(4, 1, 2), "part");
    REIS(0, ref_part_implicit(4, 1, 3), "part");
  }

  { /* four parts: two with 3, two with 2 */
    REIS(3, ref_part_large_part_size(10, 4), "large");
    REIS(2, ref_part_small_part_size(10, 4), "large");
    REIS(2, ref_part_n_large_part(10, 4), "large");

    REIS(0, ref_part_first(10, 4, 0), "first");
    REIS(3, ref_part_first(10, 4, 1), "first");
    REIS(6, ref_part_first(10, 4, 2), "first");
    REIS(8, ref_part_first(10, 4, 3), "first");
    REIS(10, ref_part_first(10, 4, 4), "first");

    REIS(0, ref_part_implicit(10, 4, 0), "part");
    REIS(0, ref_part_implicit(10, 4, 1), "part");
    REIS(0, ref_part_implicit(10, 4, 2), "part");

    REIS(1, ref_part_implicit(10, 4, 3), "part");
    REIS(1, ref_part_implicit(10, 4, 4), "part");
    REIS(1, ref_part_implicit(10, 4, 5), "part");

    REIS(2, ref_part_implicit(10, 4, 6), "part");
    REIS(2, ref_part_implicit(10, 4, 7), "part");

    REIS(3, ref_part_implicit(10, 4, 8), "part");
    REIS(3, ref_part_implicit(10, 4, 9), "part");
  }

  { /* four parts: one with 3, three with 2 */
    REIS(3, ref_part_large_part_size(9, 4), "large");
    REIS(2, ref_part_small_part_size(9, 4), "large");
    REIS(1, ref_part_n_large_part(9, 4), "large");

    REIS(0, ref_part_first(9, 4, 0), "first");
    REIS(3, ref_part_first(9, 4, 1), "first");
    REIS(5, ref_part_first(9, 4, 2), "first");
    REIS(7, ref_part_first(9, 4, 3), "first");
    REIS(9, ref_part_first(9, 4, 4), "first");

    REIS(0, ref_part_implicit(9, 4, 0), "part");
    REIS(0, ref_part_implicit(9, 4, 1), "part");
    REIS(0, ref_part_implicit(9, 4, 2), "part");

    REIS(1, ref_part_implicit(9, 4, 3), "part");
    REIS(1, ref_part_implicit(9, 4, 4), "part");

    REIS(2, ref_part_implicit(9, 4, 5), "part");
    REIS(2, ref_part_implicit(9, 4, 6), "part");

    REIS(3, ref_part_implicit(9, 4, 7), "part");
    REIS(3, ref_part_implicit(9, 4, 8), "part");
  }
  { /* 800 things on 64 parts */
    REIS(13, ref_part_large_part_size(800, 64), "large");
    REIS(12, ref_part_small_part_size(800, 64), "large");
    REIS(32, ref_part_n_large_part(800, 64), "large");

    REIS(63, ref_part_implicit(800, 64, 799), "part");
  }

  { /* part tet brick lb8.ugrid */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.lb8.ugrid";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }
    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");
    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part tet brick b8.ugrid */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.b8.ugrid";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }
    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");
    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part tet brick lb8l.ugrid */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.lb8l.ugrid";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }
    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");
    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part tet brick b8l.ugrid */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.b8l.ugrid";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }
    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");
    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part tet brick lb8.ugrid64 */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.lb8.ugrid64";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }
    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");
    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part tet brick b8.ugrid64 */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.b8.ugrid64";
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }
    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");
    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part prism b8.ugrid */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.b8.ugrid";

    RSS(ref_fixture_pri_stack_grid(&export_grid, ref_mpi), "set up tet");
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
    }
    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part prism lb8.ugrid, world comm */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.lb8.ugrid";

    RSS(ref_fixture_pri_stack_grid(&export_grid, ref_mpi), "set up tet");
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
    }
    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part tet lb8.ugrid, split comm */
    REF_MPI split_mpi;
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test_split.lb8.ugrid";
    REF_INT dummy = 0;

    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_export_by_extension(export_grid, grid_file), "export");
    }
    RSS(ref_grid_free(export_grid), "free");

    RSS(ref_mpi_bcast(ref_mpi, &dummy, 1, REF_INT_TYPE),
        "force sync before split");

    RSS(ref_mpi_half_comm(ref_mpi, &split_mpi), "split");
    RSS(ref_part_by_extension(&import_grid, split_mpi, grid_file), "import");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_mpi_join_comm(split_mpi), "join");
    RSS(ref_mpi_free(split_mpi), "free");

    RSS(ref_mpi_bcast(ref_mpi, &dummy, 1, REF_INT_TYPE),
        "force sync after join");

    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part meshb, version 2 */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test_ver2.meshb";

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      ref_grid_meshb_version(export_grid) = 2;
      RSS(ref_export_meshb(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }

    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");

    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part meshb, version 3 */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test_ver3.meshb";

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      ref_grid_meshb_version(export_grid) = 3;
      RSS(ref_export_meshb(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }

    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");

    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part meshb, version 4 */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test_ver4.meshb";

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      ref_grid_meshb_version(export_grid) = 4;
      RSS(ref_export_meshb(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }

    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");

    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* part meshb with cad_data */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.meshb";
    REF_GEOM ref_geom;
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
      ref_geom = ref_grid_geom(export_grid);
      ref_geom_cad_data_size(ref_geom) = 3;
      ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                        ref_geom_cad_data_size(ref_geom), REF_BYTE);
      ref_geom_cad_data(ref_geom)[0] = 5;
      ref_geom_cad_data(ref_geom)[1] = 4;
      ref_geom_cad_data(ref_geom)[2] = 3;
      RSS(ref_export_meshb(export_grid, grid_file), "export");
      RSS(ref_grid_free(export_grid), "free");
    }

    RSS(ref_part_by_extension(&import_grid, ref_mpi, grid_file), "import");

    ref_geom = ref_grid_geom(import_grid);
    REIS(3, ref_geom_cad_data_size(ref_geom), "cad size");
    REIS(5, ref_geom_cad_data(ref_geom)[0], "cad[0]");
    REIS(4, ref_geom_cad_data(ref_geom)[1], "cad[1]");
    REIS(3, ref_geom_cad_data(ref_geom)[2], "cad[2]");

    RSS(ref_grid_free(import_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(grid_file), "test clean up");
  }

  { /* metric */
    REF_GRID ref_grid;
    char metric_file[] = "ref_part_test.metric";

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");

    if (ref_mpi_once(ref_mpi)) {
      REF_INT node;
      FILE *file;
      file = fopen(metric_file, "w");
      RNS(file, "unable to open file");

      for (node = 0; node < 4; node++)
        fprintf(file, " %f %f %f %f %f %f\n", 1.0, 0.0, 0.0, 2.0, 0.0, 4.0);

      fclose(file);
    }

    RSS(ref_part_metric(ref_grid_node(ref_grid), metric_file), "metric");

    RSS(ref_grid_free(ref_grid), "free");

    if (ref_mpi_once(ref_mpi)) REIS(0, remove(metric_file), "test clean up");
  }

  { /* gather/part .solb by extension, version 2 */
    REF_GRID ref_grid;
    REF_INT ldim;
    REF_DBL *scalar;
    const char **scalar_names = NULL;
    char meshb[] = "ref_part_test.meshb";
    char solb[] = "ref_part_test.solb";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(ref_grid) = 2;
    RSS(ref_gather_by_extension(ref_grid, meshb), "gather meshb");
    ldim = 2;
    ref_malloc_init(scalar, ldim * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 1.0);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, scalar, scalar_names,
                                       solb),
        "gather solb");
    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, meshb), "part meshb");
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, solb),
        "part solb");
    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(solb), "test clean up");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(meshb), "test clean up");
  }

  { /* gather/part .solb by extension, version 3 */
    REF_GRID ref_grid;
    REF_INT ldim;
    REF_DBL *scalar;
    const char **scalar_names = NULL;
    char meshb[] = "ref_part_test.meshb";
    char solb[] = "ref_part_test.solb";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(ref_grid) = 3;
    RSS(ref_gather_by_extension(ref_grid, meshb), "gather meshb");
    ldim = 2;
    ref_malloc_init(scalar, ldim * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 1.0);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, scalar, scalar_names,
                                       solb),
        "gather solb");
    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, meshb), "part meshb");
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, solb),
        "part solb");
    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(solb), "test clean up");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(meshb), "test clean up");
  }

  { /* gather/part .solb by extension, version 4 */
    REF_GRID ref_grid;
    REF_INT ldim;
    REF_DBL *scalar;
    const char **scalar_names = NULL;
    char meshb[] = "ref_part_test.meshb";
    char solb[] = "ref_part_test.solb";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(ref_grid) = 4;
    RSS(ref_gather_by_extension(ref_grid, meshb), "gather meshb");
    ldim = 2;
    ref_malloc_init(scalar, ldim * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 1.0);
    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, scalar, scalar_names,
                                       solb),
        "gather solb");
    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, meshb), "part meshb");
    RSS(ref_part_scalar(ref_grid_node(ref_grid), &ldim, &scalar, solb),
        "part solb");
    ref_free(scalar);
    RSS(ref_grid_free(ref_grid), "free");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(solb), "test clean up");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(meshb), "test clean up");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
