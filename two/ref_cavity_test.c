#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_cavity.h"

#include "ref_grid.h"
#include  "ref_node.h"
#include   "ref_sort.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include   "ref_list.h"
#include  "ref_cell.h"
#include   "ref_adj.h"
#include "ref_fixture.h"


int main( void )
{

  { /* init */
    REF_CAVITY ref_cavity;
    REIS(REF_NULL, ref_cavity_free(NULL),"dont free NULL");
    RSS(ref_cavity_create(&ref_cavity,2),"create");
    REIS( 0, ref_cavity_n(ref_cavity), "init no cavity");
    REIS( 2, ref_cavity_node_per(ref_cavity), "init per");
    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add face increments count */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert");
    REIS( 1, ref_cavity_n(ref_cavity), "init no cavity");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add faces, force realloc */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];
    REF_INT f, n;

    RSS(ref_cavity_create(&ref_cavity,2),"create");

    n = ref_cavity_max(ref_cavity) + 3;
    for (f=0;f<n;f++)
      {
	nodes[0]=f;nodes[1]=f+1;
	RSS(ref_cavity_insert(ref_cavity,nodes),"insert");
	REIS( f+1, ref_cavity_n(ref_cavity), "init no cavity");
      }

    REIS( n, ref_cavity_n(ref_cavity), "count");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add same face, raise error */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");
    REIS(REF_INVALID,ref_cavity_insert(ref_cavity,nodes),"insert second");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add opposite face, mutual destruction */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");
    nodes[0]=2;nodes[1]=1;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert opposite");

    REIS( 0, ref_cavity_n(ref_cavity), "cancel");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* find face */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];
    REF_INT face;
    REF_BOOL reversed;

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");

    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find same");
    REIS(0,face,"found");
    REIS(REF_FALSE,reversed,"not rev");

    nodes[0]=2;nodes[1]=1;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find reversed");
    REIS(0,face,"found");
    REIS(REF_TRUE,reversed,"not rev");

    nodes[0]=3;nodes[1]=4;
    REIS(REF_NOT_FOUND,ref_cavity_find(ref_cavity,nodes,
				      &face,&reversed),"missing");
    REIS(REF_EMPTY,face,"found");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add triangle */
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];
    REF_INT face;
    REF_BOOL reversed;

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0]=1;nodes[1]=2;nodes[2]=3;
    RSS(ref_cavity_add_tri(ref_cavity,nodes),"insert first");

    nodes[0]=1;nodes[1]=2;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find reversed");

    nodes[0]=2;nodes[1]=3;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find reversed");

    nodes[0]=3;nodes[1]=1;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find reversed");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* insert node */
    REF_GRID ref_grid;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri" );

    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
