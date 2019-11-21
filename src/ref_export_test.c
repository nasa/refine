
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

#include "ref_export.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_fixture.h"
#include "ref_geom.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  if (2 == argc) {
    REF_GRID ref_grid;
    char file[] = "ref_export_test.tec";
    ref_mpi_stopwatch_start(ref_mpi);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "examine header");
    ref_mpi_stopwatch_stop(ref_mpi, "import");
    RSS(ref_export_tec_surf(ref_grid, file), "export");
    ref_mpi_stopwatch_stop(ref_mpi, "export");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (3 == argc) {
    REF_GRID ref_grid;
    RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi), "set up brick");
    RSS(ref_metric_olympic_node(ref_grid_node(ref_grid), 0.0001), "oly");
    RSS(ref_metric_twod_node(ref_grid_node(ref_grid)), "2d");

    RSS(ref_export_twod_msh(ref_grid, argv[1]), "export grid");
    RSS(ref_export_metric2d(ref_grid, argv[2]), "export m2d");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* export .vtk tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.vtk";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_vtk(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .tec tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.tec";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_tec(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .tec hex */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.tec";
    RSS(ref_fixture_hex_grid(&ref_grid, ref_mpi), "set up hex");
    RSS(ref_export_tec(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .fgrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.fgrid";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_fgrid(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.ugrid";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .b8.ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.b8.ugrid";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .lb8.ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.lb8.ugrid";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .b8l.ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.b8l.ugrid";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .lb8l.ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.lb8l.ugrid";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .b8.ugrid64 tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.b8.ugrid64";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export .lb8.ugrid64 tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.lb8.ugrid64";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* face id flag range */
    REF_GRID ref_grid;
    REF_INT min_faceid, max_faceid;
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");
    REIS(10, min_faceid, "min");
    REIS(10, max_faceid, "max");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* export .cogsg tet */
    REF_GRID ref_grid;
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_cogsg(ref_grid, "ref_export_test.cogsg"), "export");
    RSS(ref_grid_free(ref_grid), "free");
    REIS(0, remove("ref_export_test.cogsg"), "test clean up");
    REIS(0, remove("ref_export_test.bc"), "test clean up");
  }

  if (REF_FALSE) /* removes gnuplot dependency from unit tests */
  {              /* export .eps pri */
    REF_GRID ref_grid;
    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_eps(ref_grid, "ref_export_test.eps"), "export");
    REIS(0, remove("ref_export_test.eps"), "test clean up");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* export .html */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.html";
    RSS(ref_fixture_pri_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_html(ref_grid, file), "export");
    REIS(0, remove(file), "test clean up");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* export .meshb */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.meshb";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_meshb(ref_grid, file), "export");
    REIS(0, remove(file), "test clean up");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* export twod .msh */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.msh";
    RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi), "set up pri brick");
    RSS(ref_export_twod_msh(ref_grid, file), "export");
    REIS(0, remove(file), "test clean up");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* export twod .meshb */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.meshb";
    RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi), "set up pri brick");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    REIS(0, remove(file), "test clean up");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* export binary tecplot .plt */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.plt";
    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_plt(ref_grid, file), "export");
    REIS(0, remove(file), "test clean up");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* export cubic edge .meshb */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.meshb";
    RSS(ref_fixture_twod_cubic_edge(&ref_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(ref_grid, file), "export");
    REIS(0, remove(file), "test clean up");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
