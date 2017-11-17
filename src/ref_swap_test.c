
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

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_sort.h"
#include "ref_matrix.h"

#include "ref_mpi.h"

#include "ref_swap.h"

#include "ref_fixture.h"

int main( void )
{
  REF_MPI ref_mpi;
  RSS( ref_mpi_create( &ref_mpi ), "create" );
  
  { /* leave single faces alone */
    REF_GRID ref_grid;

    RSS(ref_fixture_tet_grid(&ref_grid,ref_mpi),"set up");

    REIS(REF_INVALID,ref_swap_remove_two_face_cell(ref_grid,0),"cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* leave different 2 faces alone */
    REF_GRID ref_grid;
    REF_INT nodes[4] = {0,3,1,50}, cell;

    RSS(ref_fixture_tet_grid(&ref_grid,ref_mpi),"set up");
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");
    
    REIS(REF_INVALID,ref_swap_remove_two_face_cell(ref_grid,0),"cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* remove same 2 faces */
    REF_GRID ref_grid;
    REF_INT nodes[4] = {0,3,1,10}, cell;

    RSS(ref_fixture_tet_grid(&ref_grid,ref_mpi),"set up");
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");
    
    RSS(ref_swap_remove_two_face_cell(ref_grid,0),"removal failed");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* leave different 3 faces alone */
    REF_GRID ref_grid;
    REF_INT nodes[4], cell;

    RSS(ref_fixture_tet_grid(&ref_grid,ref_mpi),"set up");

    nodes[0]=0;nodes[1]=3;nodes[2]=1;nodes[3]=50;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");

    nodes[0]=0;nodes[1]=2;nodes[2]=3;nodes[3]=50;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");
    
    REIS(REF_INVALID,ref_swap_remove_three_face_cell(ref_grid,0),"cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* remove same 3 faces */
    REF_GRID ref_grid;
    REF_INT nodes[4], cell;

    RSS(ref_fixture_tet_grid(&ref_grid,ref_mpi),"set up");
    nodes[0]=0;nodes[1]=3;nodes[2]=1;nodes[3]=10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");

    nodes[0]=0;nodes[1]=2;nodes[2]=3;nodes[3]=10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");

    RSS(ref_swap_remove_three_face_cell(ref_grid,0),"removal failed");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)),"tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    REIS(3, ref_node_n(ref_grid_node(ref_grid)),"nodes");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }
  
  RSS( ref_mpi_free( ref_mpi ), "free" );
  return 0;
}
