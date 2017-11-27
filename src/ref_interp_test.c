
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_interp.h"

#include "ref_mpi.h"
#include "ref_fixture.h"
#include "ref_grid.h"
#include  "ref_cell.h"
#include   "ref_adj.h"
#include  "ref_node.h"
#include   "ref_matrix.h"
#include   "ref_math.h"
#include   "ref_sort.h"
#include   "ref_list.h"
#include  "ref_geom.h"
#include   "ref_dict.h"
#include   "ref_export.h"
#include  "ref_gather.h"
#include  "ref_adapt.h"
#include   "ref_collapse.h"
#include    "ref_edge.h"
#include   "ref_split.h"
#include    "ref_cavity.h"
#include    "ref_twod.h"
#include   "ref_smooth.h"
#include    "ref_clump.h"

int main( int argc, char *argv[] )
{
  REF_MPI ref_mpi;
  RSS( ref_mpi_start( argc, argv ), "start" );
  RSS( ref_mpi_create( &ref_mpi ), "make mpi" );

  {
    REF_INTERP ref_interp;
    RSS( ref_interp_create( &ref_interp ), "make interp" );
    REIS(0, ref_interp_steps(ref_interp), "steps" );
    RSS( ref_interp_free( ref_interp ), "interp free" );
  }

  {
    REF_GRID from, to;
    RSS( ref_fixture_twod_brick_grid( &from, ref_mpi ), "brick" );
    RSS( ref_fixture_twod_brick_grid( &to, ref_mpi ), "brick" );
    RSS( ref_grid_free(to),"free");
    RSS( ref_grid_free(from),"free");
  }


  RSS( ref_mpi_free( ref_mpi ), "mpi free" );
  RSS( ref_mpi_stop( ), "stop" );
  return 0;
}
