
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
#include <math.h>

#include "ref_elast.h"
#include "ref_mpi.h"
#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_sort.h"
#include "ref_matrix.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_geom.h"
#include "ref_dict.h"
#include "ref_gather.h"
#include "ref_export.h"
#include "ref_adj.h"
#include "ref_adapt.h"
#include "ref_edge.h"
#include "ref_split.h"
#include "ref_collapse.h"
#include "ref_smooth.h"
#include "ref_cavity.h"
#include "ref_subdiv.h"
#include "ref_twod.h"
#include "ref_clump.h"

#include "ref_fixture.h"

int main( int argc, char *argv[] )
{
  REF_MPI ref_mpi;
  RSS( ref_mpi_start( argc, argv ), "start" );
  RSS( ref_mpi_create( &ref_mpi ), "make mpi" );
 
  {  /* init */
    REF_GRID ref_grid;
    REF_ELAST ref_elast;
    REF_INT node;
    REF_DBL dxyz[3];
    REF_DBL l2norm;
    REF_INT sweep;
    
    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");
    RSS(ref_elast_create(&ref_elast,ref_grid),"create");

    node=0; dxyz[0] = 0.0; dxyz[1] = 0.0; dxyz[2] = 1.0;
    RSS(ref_elast_displace(ref_elast,node,dxyz),"create");
    node=1; dxyz[0] = 0.0; dxyz[1] = 0.0; dxyz[2] = 1.0;
    RSS(ref_elast_displace(ref_elast,node,dxyz),"create");
    node=2; dxyz[0] = 0.0; dxyz[1] = 0.0; dxyz[2] = 1.0;
    RSS(ref_elast_displace(ref_elast,node,dxyz),"create");
    
    RSS(ref_elast_assemble(ref_elast),"elast");
    for (sweep=0;sweep<0;sweep++)
      {
        RSS(ref_elast_relax(ref_elast,&l2norm),"elast");
        printf("res %e\n",l2norm);
      }
    RSS(ref_elast_free(ref_elast),"elast");
    RSS(ref_grid_free(ref_grid),"free");
  }

  RSS( ref_mpi_free( ref_mpi ), "mpi free" );
  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
