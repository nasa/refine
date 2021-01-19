
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
#include <string.h>

#include "ref_adj.h"
#include "ref_args.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_face.h"
#include "ref_fixture.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_sort.h"

static REF_STATUS ref_validation_lb8_ugrid_volume(const char *filename) {
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_FILEPOS conn_offset, ibyte = 4;
  REF_DBL *xyz, *xyzs[4];
  REF_INT cell, nodes[4];
  REF_DBL volume, min_volume, max_volume;
  printf("open %s\n", filename);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  REIS(1, fread(&nnode, sizeof(REF_INT), 1, file), "nnode");
  REIS(1, fread(&ntri, sizeof(REF_INT), 1, file), "ntri");
  REIS(1, fread(&nqua, sizeof(REF_INT), 1, file), "nqua");
  REIS(1, fread(&ntet, sizeof(REF_INT), 1, file), "ntet");
  REIS(1, fread(&npyr, sizeof(REF_INT), 1, file), "npyr");
  REIS(1, fread(&npri, sizeof(REF_INT), 1, file), "npri");
  REIS(1, fread(&nhex, sizeof(REF_INT), 1, file), "nhex");

  printf("nnode %d\n", nnode);

  ref_malloc(xyz, 3 * nnode, REF_DBL);
  REIS(3 * nnode, fread(xyz, sizeof(REF_DBL), (size_t)(3 * nnode), file),
       "xyz");

  printf("tet %d\n", ntet);
  conn_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                (REF_FILEPOS)ntri * 4 * ibyte + (REF_FILEPOS)nqua * 5 * ibyte;
  REIS(0, fseeko(file, conn_offset, SEEK_SET), "seek tet failed");
  max_volume = REF_DBL_MIN;
  min_volume = REF_DBL_MAX;
  for (cell = 0; cell < ntet; cell++) {
    REIS(4, fread(nodes, sizeof(REF_INT), 4, file), "tet");
    xyzs[0] = &(xyz[3 * (nodes[0] - 1)]);
    xyzs[1] = &(xyz[3 * (nodes[1] - 1)]);
    xyzs[2] = &(xyz[3 * (nodes[2] - 1)]);
    xyzs[3] = &(xyz[3 * (nodes[3] - 1)]);
    RSS(ref_node_xyz_vol(xyzs, &volume), "vol");
    min_volume = MIN(min_volume, volume);
    max_volume = MAX(max_volume, volume);
  }
  printf("volume range %e %e\n", min_volume, max_volume);
  ref_free(xyz);
  fclose(file);
  return REF_SUCCESS;
}

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_INT pos = REF_EMPTY;

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");

  RXS(ref_args_find(argc, argv, "--const", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) {
    REF_GRID ref_grid;
    REF_INT ldim = 4;
    REF_DBL *counts;
    REF_INT *total, *four, *three;
    REF_INT part, hits, node, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
    REF_NODE ref_node;
    REF_CELL ref_cell;
    const char *varnames[] = {"degree", "three", "four", "threeratio"};

    REIS(1, pos, " ref_validation_test --const input.meshb");
    REIS(3, argc, " ref_validation_test --const input.meshb");

    if (ref_mpi_para(ref_mpi)) {
      if (ref_mpi_once(ref_mpi)) printf("part %s\n", argv[2]);
      RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[2]), "part");
      ref_mpi_stopwatch_stop(ref_mpi, "part");
    } else {
      if (ref_mpi_once(ref_mpi)) printf("import %s\n", argv[2]);
      RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "import");
      ref_mpi_stopwatch_stop(ref_mpi, "import");
    }
    ref_node = ref_grid_node(ref_grid);
    ref_cell = ref_grid_tet(ref_grid);
    ref_malloc_init(counts, ldim * ref_node_max(ref_grid_node(ref_grid)),
                    REF_DBL, 0);
    ref_malloc_init(total, ref_node_max(ref_grid_node(ref_grid)), REF_INT, 0);
    ref_malloc_init(four, ref_node_max(ref_grid_node(ref_grid)), REF_INT, 0);
    ref_malloc_init(three, ref_node_max(ref_grid_node(ref_grid)), REF_INT, 0);

    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) == part) {
        hits = 0;
        each_ref_cell_cell_node(ref_cell, cell_node) {
          if (!ref_cell_node_empty(ref_grid_tri(ref_grid), nodes[cell_node]))
            hits++;
        }
        if (4 == hits) {
          each_ref_cell_cell_node(ref_cell, cell_node) {
            four[nodes[cell_node]] += 1;
          }
        }
        if (3 == hits) {
          each_ref_cell_cell_node(ref_cell, cell_node) {
            three[nodes[cell_node]] += 1;
          }
        }
        each_ref_cell_cell_node(ref_cell, cell_node) {
          total[nodes[cell_node]] += 1;
        }
      }
    }

    RSS(ref_node_localize_ghost_int(ref_node, total), "total");
    RSS(ref_node_localize_ghost_int(ref_node, three), "total");
    RSS(ref_node_localize_ghost_int(ref_node, four), "total");

    each_ref_node_valid_node(ref_node, node) {
      if (ref_node_owned(ref_node, node)) {
        counts[0 + ldim * node] = (REF_DBL)total[node];
        counts[1 + ldim * node] = (REF_DBL)three[node];
        counts[2 + ldim * node] = (REF_DBL)four[node];
        if (0 < total[node]) {
          counts[3 + ldim * node] = (REF_DBL)four[node] / (REF_DBL)total[node];
        }
      }
    }

    RSS(ref_node_ghost_dbl(ref_grid_node(ref_grid), counts, ldim),
        "update ghosts");

    ref_mpi_stopwatch_stop(ref_mpi, "count");

    RSS(ref_gather_scalar_by_extension(ref_grid, ldim, counts, varnames,
                                       "ref_validation_counts.plt"),
        "plt");

    ref_mpi_stopwatch_stop(ref_mpi, "gather");

    ref_free(four);
    ref_free(three);
    ref_free(total);
    ref_free(counts);
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--vol", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) {
    REIS(1, pos, " ref_validation_test --vol grid.lb8.ugrid");
    REIS(3, argc, " ref_validation_test --vol grid.lb8.ugrid");
    RSS(ref_validation_lb8_ugrid_volume(argv[2]), "ugrid vol");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--repair", &pos), REF_NOT_FOUND, "arg search");
  if (REF_EMPTY != pos) {
    REF_GRID ref_grid;
    REIS(1, pos, " ref_validation_test --repair grid.ext");
    REIS(3, argc, " ref_validation_test --repair grid.ext");
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[2]), "ugrid vol");
    RSS(ref_validation_repair(ref_grid), "ugrid vol");
    RSS(ref_export_by_extension(ref_grid, "fixed.tec"), "ugrid vol");
    RSS(ref_grid_free(ref_grid), "free");
    RSS(ref_mpi_free(ref_mpi), "free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  if (argc > 1) {
    REF_GRID ref_grid;
    printf("validating\n");

    printf("reading %s\n", argv[1]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "from ugrid");
    printf("complete.\n");

    RSS(ref_grid_inspect(ref_grid), "inspection");

    if (ref_grid_twod(ref_grid)) {
      RSS(ref_validation_twod_orientation(ref_grid), "twod tri orientation");
    } else {
      printf("validate.\n");
      RSS(ref_validation_volume_status(ref_grid), "tet volume grid");
      RSS(ref_validation_all(ref_grid), "invalid grid");
    }
    RSS(ref_grid_free(ref_grid), "free");
    printf("done.\n");
  }

  if (!ref_mpi_para(ref_mpi)) {
    REF_GRID ref_grid;
    RSS(ref_fixture_twod_brick_grid(&ref_grid, ref_mpi, 4), "twod brick");
    RSS(ref_validation_twod_orientation(ref_grid), "twod tri orientation");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
