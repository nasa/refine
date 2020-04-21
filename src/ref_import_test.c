
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

#include "ref_import.h"

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
#include "ref_geom.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT pos;
  REF_BOOL transmesh = REF_FALSE;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--transmesh", &pos), REF_NOT_FOUND,
      "arg search");
  if (REF_EMPTY != pos) transmesh = REF_TRUE;

  if (2 == argc && !transmesh) {
    RSS(ref_import_examine_header(argv[1]), "examine header");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--prop", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos && pos == 1 && argc == 3) {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CELL ref_cell;
    REF_INT i, cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_DBL total_area, total_normal[3], center[3];
    REF_DBL area, normal[3], centroid[3];
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "import");
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tri(ref_grid);
    total_area = 0;
    for (i = 0; i < 3; i++) center[i] = 0;
    for (i = 0; i < 3; i++) total_normal[i] = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
      total_area += area;
      RSS(ref_node_tri_centroid(ref_node, nodes, centroid), "area");
      for (i = 0; i < 3; i++) center[i] += area * centroid[i];
      RSS(ref_node_tri_normal(ref_node, nodes, normal), "area");
      for (i = 0; i < 3; i++) total_normal[i] += normal[i];
    }
    for (i = 0; i < 3; i++) center[i] /= total_area;
    RSS(ref_math_normalize(total_normal), "norm");
    printf("%d tri\n", ref_cell_n(ref_cell));
    printf("%f %f %f center\n", center[0], center[1], center[2]);
    printf("%f %f %f normal\n", total_normal[0], total_normal[1],
           total_normal[2]);
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* export import twod .msh brick */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.msh";
    RSS(ref_fixture_twod_brick_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_twod_msh(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import twod .meshb brick */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test-2d.meshb";
    RSS(ref_fixture_twod_brick_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import twod cubic edge brick */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test-ed3.meshb";
    RSS(ref_fixture_twod_cubic_edge(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_ed3(export_grid)),
         ref_cell_n(ref_grid_ed3(import_grid)), "ed3 count");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .meshb tet brick, version 2 */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test_ver2.meshb";
    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(export_grid) = 2;
    RSS(ref_export_by_extension(export_grid, file), "export");
    if (transmesh) {
      REIS(
          0,
          system(
              "transmesh ref_import_test_ver2.meshb ref_import_test_ver2.mesh"),
          "mesh");
      REIS(
          0,
          system(
              "transmesh ref_import_test_ver2.mesh ref_import_test_ver2.meshb"),
          "meshb");
    }
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    if (!transmesh) REIS(0, remove(file), "test clean up");
  }

  { /* export import .meshb tet brick, version 3 */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test_ver3.meshb";
    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(export_grid) = 3;
    RSS(ref_export_by_extension(export_grid, file), "export");
    if (transmesh) {
      REIS(
          0,
          system(
              "transmesh ref_import_test_ver3.meshb ref_import_test_ver3.mesh"),
          "mesh");
      REIS(
          0,
          system(
              "transmesh ref_import_test_ver3.mesh ref_import_test_ver3.meshb"),
          "meshb");
    }
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    if (!transmesh) REIS(0, remove(file), "test clean up");
  }

  { /* export import .meshb tet brick, version 4 */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test_ver4.meshb";
    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(export_grid) = 4;
    RSS(ref_export_by_extension(export_grid, file), "export");
    if (transmesh) {
      REIS(
          0,
          system(
              "transmesh ref_import_test_ver4.meshb ref_import_test_ver4.mesh"),
          "mesh");
      REIS(
          0,
          system(
              "transmesh ref_import_test_ver4.mesh ref_import_test_ver4.meshb"),
          "meshb");
    }
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    if (!transmesh) REIS(0, remove(file), "test clean up");
  }

  { /* export import .meshb tet brick with cad_model, default */
    REF_GRID export_grid, import_grid;
    REF_GEOM ref_geom;
    char file[] = "ref_import_test.meshb";
    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    ref_geom = ref_grid_geom(export_grid);
    ref_geom_cad_data_size(ref_geom) = 5;
    ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                      ref_geom_cad_data_size(ref_geom), REF_BYTE);
    ref_geom_cad_data(ref_geom)[0] = 5;
    ref_geom_cad_data(ref_geom)[1] = 4;
    ref_geom_cad_data(ref_geom)[2] = 3;
    ref_geom_cad_data(ref_geom)[3] = 2;
    ref_geom_cad_data(ref_geom)[4] = 1;
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    ref_geom = ref_grid_geom(import_grid);
    REIS(5, ref_geom_cad_data_size(ref_geom), "cad size");
    REIS(5, ref_geom_cad_data(ref_geom)[0], "cad[0]");
    REIS(4, ref_geom_cad_data(ref_geom)[1], "cad[1]");
    REIS(3, ref_geom_cad_data(ref_geom)[2], "cad[2]");
    REIS(2, ref_geom_cad_data(ref_geom)[3], "cad[3]");
    REIS(1, ref_geom_cad_data(ref_geom)[4], "cad[4]");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .meshb tet brick with cad_model, version 4  */
    REF_GRID export_grid, import_grid;
    REF_GEOM ref_geom;
    char file[] = "ref_import_test.meshb";
    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(export_grid) = 4;
    ref_geom = ref_grid_geom(export_grid);
    ref_geom_cad_data_size(ref_geom) = 5;
    ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                      ref_geom_cad_data_size(ref_geom), REF_BYTE);
    ref_geom_cad_data(ref_geom)[0] = 5;
    ref_geom_cad_data(ref_geom)[1] = 4;
    ref_geom_cad_data(ref_geom)[2] = 3;
    ref_geom_cad_data(ref_geom)[3] = 2;
    ref_geom_cad_data(ref_geom)[4] = 1;
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    ref_geom = ref_grid_geom(import_grid);
    REIS(5, ref_geom_cad_data_size(ref_geom), "cad size");
    REIS(5, ref_geom_cad_data(ref_geom)[0], "cad[0]");
    REIS(4, ref_geom_cad_data(ref_geom)[1], "cad[1]");
    REIS(3, ref_geom_cad_data(ref_geom)[2], "cad[2]");
    REIS(2, ref_geom_cad_data(ref_geom)[3], "cad[3]");
    REIS(1, ref_geom_cad_data(ref_geom)[4], "cad[4]");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .meshb tet brick with geom, default */
    REF_GRID export_grid, import_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_GEOM ref_geom;
    REF_INT type, id, node;
    REF_DBL param[2];
    char file[] = "ref_import_test_geom2.meshb";
    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    ref_geom = ref_grid_geom(export_grid);
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 15;
    RSS(ref_cell_add(ref_grid_edg(export_grid), nodes, &cell), "add edge");
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

    RSS(ref_export_by_extension(export_grid, file), "export");
    if (transmesh) {
      REIS(0,
           system("transmesh ref_import_test_geom2.meshb "
                  "ref_import_test_geom2.mesh"),
           "mesh");
      REIS(0,
           system("transmesh ref_import_test_geom2.mesh "
                  "ref_import_test_geom2.meshb"),
           "meshb");
    }
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");

    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_edg(export_grid)),
         ref_cell_n(ref_grid_edg(import_grid)), "edg count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");

    REIS(ref_geom_n(ref_grid_geom(export_grid)),
         ref_geom_n(ref_grid_geom(import_grid)), "tet count");

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    if (!transmesh) REIS(0, remove(file), "test clean up");
  }

  { /* export import .meshb tet brick with geom, version 4 */
    REF_GRID export_grid, import_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    REF_GEOM ref_geom;
    REF_INT type, id, node;
    REF_DBL param[2];
    char file[] = "ref_import_test_geom4.meshb";
    RSS(ref_fixture_tet_brick_grid(&export_grid, ref_mpi), "set up tet");
    ref_grid_meshb_version(export_grid) = 4;
    ref_geom = ref_grid_geom(export_grid);
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 15;
    RSS(ref_cell_add(ref_grid_edg(export_grid), nodes, &cell), "add edge");
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

    RSS(ref_export_by_extension(export_grid, file), "export");
    if (transmesh) {
      REIS(0,
           system("transmesh ref_import_test_geom4.meshb "
                  "ref_import_test_geom4.mesh"),
           "mesh");
      REIS(0,
           system("transmesh ref_import_test_geom4.mesh "
                  "ref_import_test_geom4.meshb"),
           "meshb");
    }
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");

    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_edg(export_grid)),
         ref_cell_n(ref_grid_edg(import_grid)), "edg count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");

    REIS(ref_geom_n(ref_grid_geom(export_grid)),
         ref_geom_n(ref_grid_geom(import_grid)), "tet count");

    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    if (!transmesh) REIS(0, remove(file), "test clean up");
  }

  { /* export import .fgrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.fgrid";
    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.ugrid";
    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .lb8.ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.lb8.ugrid";
    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RWDS(ref_node_xyz(ref_grid_node(export_grid), 0, 1),
         ref_node_xyz(ref_grid_node(import_grid), 0, 1), 1e-15, "x 1");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 0, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 0, 0), "tet node0");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 1, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 1, 0), "tet node 1");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .b8.ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.b8.ugrid";
    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RWDS(ref_node_xyz(ref_grid_node(export_grid), 0, 1),
         ref_node_xyz(ref_grid_node(import_grid), 0, 1), 1e-15, "x 1");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 0, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 0, 0), "tet node0");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 1, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 1, 0), "tet node 1");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .lb8l.ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.lb8l.ugrid";
    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RWDS(ref_node_xyz(ref_grid_node(export_grid), 0, 1),
         ref_node_xyz(ref_grid_node(import_grid), 0, 1), 1e-15, "x 1");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 0, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 0, 0), "tet node0");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 1, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 1, 0), "tet node 1");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .b8l.ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.b8l.ugrid";
    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RWDS(ref_node_xyz(ref_grid_node(export_grid), 0, 1),
         ref_node_xyz(ref_grid_node(import_grid), 0, 1), 1e-15, "x 1");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 0, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 0, 0), "tet node0");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 1, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 1, 0), "tet node 1");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .lb8.ugrid64 tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.lb8.ugrid64";
    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RWDS(ref_node_xyz(ref_grid_node(export_grid), 0, 1),
         ref_node_xyz(ref_grid_node(import_grid), 0, 1), 1e-15, "x 1");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 0, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 0, 0), "tet node0");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 1, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 1, 0), "tet node 1");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* export import .b8.ugrid64 tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.b8.ugrid64";
    RSS(ref_fixture_tet_grid(&export_grid, ref_mpi), "set up tet");
    RSS(ref_export_by_extension(export_grid, file), "export");
    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");
    REIS(ref_node_n(ref_grid_node(export_grid)),
         ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(ref_cell_n(ref_grid_tri(export_grid)),
         ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(ref_cell_n(ref_grid_qua(export_grid)),
         ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    REIS(ref_cell_n(ref_grid_tet(export_grid)),
         ref_cell_n(ref_grid_tet(import_grid)), "tet count");
    RWDS(ref_node_xyz(ref_grid_node(export_grid), 0, 1),
         ref_node_xyz(ref_grid_node(import_grid), 0, 1), 1e-15, "x 1");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 0, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 0, 0), "tet node0");
    REIS(ref_cell_c2n(ref_grid_tet(export_grid), 1, 0),
         ref_cell_c2n(ref_grid_tet(import_grid), 1, 0), "tet node 1");
    RSS(ref_grid_free(import_grid), "free");
    RSS(ref_grid_free(export_grid), "free");
    REIS(0, remove(file), "test clean up");
  }

  { /* read AFLR3 < 14 surf */
    REF_GRID import_grid;
    char file[] = "ref_import_test.surf";
    FILE *f;
    f = fopen(file, "w");
    fprintf(f, "1 0 3\n");
    fprintf(f, "12.0 14.0 16.0\n");
    fprintf(f, "22.0 24.0 26.0\n");
    fprintf(f, "32.0 34.0 36.0\n");
    fprintf(f, "1 2 3 1 0 1\n");
    fclose(f);

    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");

    REIS(3, ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(1, ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(0, ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    RWDS(12.0, ref_node_xyz(ref_grid_node(import_grid), 0, 0), 1e-15, "x 0");
    RWDS(36.0, ref_node_xyz(ref_grid_node(import_grid), 2, 2), 1e-15, "z 2");
    RSS(ref_grid_free(import_grid), "free");

    REIS(0, remove(file), "test clean up");
  }

  { /* read AFLR3 > 15 surf */
    REF_GRID import_grid;
    char file[] = "ref_import_test.surf";
    FILE *f;
    f = fopen(file, "w");
    fprintf(f, "1 0 3\n");
    fprintf(f, "12.0 14.0 16.0 0 0\n");
    fprintf(f, "22.0 24.0 26.0 0 0\n");
    fprintf(f, "32.0 34.0 36.0 0 0\n");
    fprintf(f, "1 2 3 1 0 0\n");
    fclose(f);

    RSS(ref_import_by_extension(&import_grid, ref_mpi, file), "import");

    REIS(3, ref_node_n(ref_grid_node(import_grid)), "node count");
    REIS(1, ref_cell_n(ref_grid_tri(import_grid)), "tri count");
    REIS(0, ref_cell_n(ref_grid_qua(import_grid)), "qua count");
    RWDS(12.0, ref_node_xyz(ref_grid_node(import_grid), 0, 0), 1e-15, "x 0");
    RWDS(36.0, ref_node_xyz(ref_grid_node(import_grid), 2, 2), 1e-15, "z 2");
    RSS(ref_grid_free(import_grid), "free");

    REIS(0, remove(file), "test clean up");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
