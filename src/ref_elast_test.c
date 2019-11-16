
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

#include "ref_elast.h"

#include <math.h>
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
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_subdiv.h"
#include "ref_twod.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  { /* init */
    REF_GRID ref_grid;
    REF_ELAST ref_elast;

    RSS(ref_grid_create(&ref_grid, ref_mpi), "create");
    RSS(ref_elast_create(&ref_elast, ref_grid), "create");

    RSS(ref_elast_assemble(ref_elast), "elast");

    RSS(ref_elast_free(ref_elast), "elast");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* tet */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_ELAST ref_elast;
    REF_INT node;
    REF_DBL dxyz[3];
    REF_DBL l2norm;
    REF_INT sweep;

    RSS(ref_fixture_tet_grid(&ref_grid, ref_mpi), "create");
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_elast_create(&ref_elast, ref_grid), "create");

    dxyz[0] = 0.0;
    dxyz[1] = 0.0;
    dxyz[2] = 1.0;
    if (REF_SUCCESS == ref_node_local(ref_node, 0, &node)) {
      if (ref_node_owned(ref_node, node))
        RSS(ref_elast_displace(ref_elast, node, dxyz), "create");
    }
    if (REF_SUCCESS == ref_node_local(ref_node, 1, &node)) {
      if (ref_node_owned(ref_node, node))
        RSS(ref_elast_displace(ref_elast, node, dxyz), "create");
    }
    if (REF_SUCCESS == ref_node_local(ref_node, 2, &node)) {
      if (ref_node_owned(ref_node, node))
        RSS(ref_elast_displace(ref_elast, node, dxyz), "create");
    }

    RSS(ref_elast_assemble(ref_elast), "elast");
    for (sweep = 0; sweep < 2; sweep++) {
      RSS(ref_elast_relax(ref_elast, &l2norm), "elast");
    }
    RWDS(0.0, l2norm, -1.0, "not coverged on to steps");
    if (REF_SUCCESS == ref_node_local(ref_node, 3, &node)) {
      RWDS(ref_elast->displacement[0 + 3 * node], 0.0, -1.0, "x");
      RWDS(ref_elast->displacement[1 + 3 * node], 0.0, -1.0, "y");
      RWDS(ref_elast->displacement[2 + 3 * node], 1.0, -1.0, "z");
    }

    RSS(ref_elast_free(ref_elast), "elast");
    RSS(ref_grid_free(ref_grid), "free");
  }

  { /* bricks */
    REF_GRID ref_grid;
    REF_ELAST ref_elast;
    REF_INT node;
    REF_DBL dxyz[3];
    REF_DBL l2norm;
    REF_INT sweep;
    char file[] = "ref_elast_test.meshb";

    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_fixture_tet_brick_grid(&ref_grid, ref_mpi), "brick");
      RSS(ref_export_by_extension(ref_grid, file), "export");
      RSS(ref_grid_free(ref_grid), "free");
    }
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, file), "import");
    if (ref_mpi_once(ref_mpi)) REIS(0, remove(file), "test clean up");

    RSS(ref_elast_create(&ref_elast, ref_grid), "create");

    each_ref_node_valid_node(
        ref_grid_node(ref_grid),
        node) if ((-0.01 < ref_node_xyz(ref_grid_node(ref_grid), 2, node) &&
                   0.01 > ref_node_xyz(ref_grid_node(ref_grid), 2, node))) {
      dxyz[0] = 0.0;
      dxyz[1] = 0.0;
      dxyz[2] = 1.0;
      RSS(ref_elast_displace(ref_elast, node, dxyz), "create");
    }

    RSS(ref_elast_assemble(ref_elast), "elast");
    for (sweep = 0; sweep < 1000; sweep++) {
      RSS(ref_elast_relax(ref_elast, &l2norm), "elast");
    }
    RWDS(0.0, l2norm, -1.0, "not coverged on to steps");
    each_ref_node_valid_node(
        ref_grid_node(ref_grid),
        node) if ((0.99 < ref_node_xyz(ref_grid_node(ref_grid), 2, node) &&
                   1.01 > ref_node_xyz(ref_grid_node(ref_grid), 2, node))) {
      RWDS(0.0, ref_elast->displacement[0 + 3 * node], -1.0, "x");
      RWDS(0.0, ref_elast->displacement[1 + 3 * node], -1.0, "y");
      RWDS(1.0, ref_elast->displacement[2 + 3 * node], -1.0, "z");
    }

    RSS(ref_elast_free(ref_elast), "elast");
    RSS(ref_grid_free(ref_grid), "free");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
