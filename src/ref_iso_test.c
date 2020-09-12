
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

#include "ref_iso.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adapt.h"
#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_collapse.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_fixture.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_metric.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_smooth.h"
#include "ref_sort.h"
#include "ref_split.h"
#include "ref_swap.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  {
    REF_GRID ref_grid, iso_grid;
    REF_NODE ref_node;
    REF_DBL *field;
    REF_INT node;
    REF_DBL offset = 0.2;

    RSS(ref_fixture_tri_grid(&ref_grid, ref_mpi), "tri");
    ref_node = ref_grid_node(ref_grid);
    ref_malloc(field, ref_node_max(ref_node), REF_DBL);
    each_ref_node_valid_node(ref_node, node) {
      field[node] = ref_node_xyz(ref_node, 0, node) - offset;
    }
    RSS(ref_iso_insert(&iso_grid, ref_grid, field), "iso");
    ref_grid_free(iso_grid);
    ref_free(field);
    ref_grid_free(ref_grid);
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
