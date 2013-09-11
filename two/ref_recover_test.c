#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_recover.h"
#include "ref_grid.h"
#include  "ref_node.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include   "ref_sort.h"
#include   "ref_list.h"
#include  "ref_cell.h"
#include   "ref_adj.h"
#include "ref_fixture.h"

int main( void )
{

  { /* init */
    REF_GRID ref_grid;
    REF_RECOVER ref_recover;
    REIS(REF_NULL, ref_recover_free(NULL),"dont free NULL");
    RSS(ref_grid_create(&ref_grid),"create");
    RSS(ref_recover_create(&ref_recover,ref_grid),"create");
    REIS( 0, ref_recover_n(ref_recover), "init no recover");
    RSS(ref_recover_free(ref_recover),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  SKIP_BLOCK("implement")  
  { /* insert a node */
    REF_GRID ref_grid;
    REF_RECOVER ref_recover;
    REF_DBL xz[2];
    REF_INT node;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri fixture" );
    RSS(ref_recover_create(&ref_recover,ref_grid),"create");

    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    xz[0] = 0.3; xz[1] = 0.3;
    RSS(ref_recover_insert_twod(ref_recover,xz,&node),"create");
 
    REIS(8, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    RSS(ref_recover_free(ref_recover),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* triangle enclosing a node */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_RECOVER ref_recover;
    REF_DBL xz[2];
    REF_INT global, node, cell;
    REF_DBL bary[3];
    REF_DBL tol = -1.0;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri fixture" );
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_recover_create(&ref_recover,ref_grid),"create");

    xz[0] = 1.0/3.0; xz[1] = 1.0/3.0;
    RSS( ref_node_next_global( ref_node, &global ), "next global");
    RSS( ref_node_add( ref_node, global, &node ), "add node");
    ref_node_xyz(ref_node,0,node) = xz[0];
    ref_node_xyz(ref_node,1,node) = 0.0;
    ref_node_xyz(ref_node,2,node) = xz[1];

    RSS(ref_recover_enclosing_triangle(ref_recover,node,&cell,bary),"create");
 
    REIS(1, cell, "tri");
    RWDS(1.0/3.0, bary[0], tol, "b[0]");
    RWDS(1.0/3.0, bary[1], tol, "b[1]");
    RWDS(1.0/3.0, bary[2], tol, "b[2]");

    RSS(ref_recover_free(ref_recover),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
