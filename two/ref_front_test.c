#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_front.h"

int main( void )
{

  { /* init */
    REF_FRONT ref_front;
    REIS(REF_NULL, ref_front_free(NULL),"dont free NULL");
    RSS(ref_front_create(&ref_front,2),"create");
    REIS( 0, ref_front_n(ref_front), "init no front");
    REIS( 2, ref_front_node_per(ref_front), "init per");
    RSS(ref_front_free(ref_front),"free");
  }

  { /* add face increments count */
    REF_FRONT ref_front;
    REF_INT nodes[2];

    RSS(ref_front_create(&ref_front,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_front_insert(ref_front,nodes),"insert");
    REIS( 1, ref_front_n(ref_front), "init no front");

    RSS(ref_front_free(ref_front),"free");
  }

  { /* add faces, force realloc */
    REF_FRONT ref_front;
    REF_INT nodes[2];
    REF_INT f, n;

    RSS(ref_front_create(&ref_front,2),"create");

    n = ref_front_max(ref_front) + 3;
    for (f=0;f<n;f++)
      {
	nodes[0]=f;nodes[1]=f+1;
	RSS(ref_front_insert(ref_front,nodes),"insert");
	REIS( f+1, ref_front_n(ref_front), "init no front");
      }

    REIS( n, ref_front_n(ref_front), "count");

    RSS(ref_front_free(ref_front),"free");
  }

  { /* add same face, raise error */
    REF_FRONT ref_front;
    REF_INT nodes[2];

    RSS(ref_front_create(&ref_front,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_front_insert(ref_front,nodes),"insert first");
    REIS(REF_INVALID,ref_front_insert(ref_front,nodes),"insert second");

    RSS(ref_front_free(ref_front),"free");
  }

  { /* add opposite face, mutual destruction */
    REF_FRONT ref_front;
    REF_INT nodes[2];

    RSS(ref_front_create(&ref_front,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_front_insert(ref_front,nodes),"insert first");
    nodes[0]=2;nodes[1]=1;
    RSS(ref_front_insert(ref_front,nodes),"insert opposite");

    REIS( 0, ref_front_n(ref_front), "cancel");

    RSS(ref_front_free(ref_front),"free");
  }

  { /* find face */
    REF_FRONT ref_front;
    REF_INT nodes[2];
    REF_INT face;
    REF_BOOL reversed;

    RSS(ref_front_create(&ref_front,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_front_insert(ref_front,nodes),"insert first");

    RSS(ref_front_find(ref_front,nodes,&face,&reversed),"find same");
    REIS(0,face,"found");
    REIS(REF_FALSE,reversed,"not rev");

    nodes[0]=2;nodes[1]=1;
    RSS(ref_front_find(ref_front,nodes,&face,&reversed),"find reversed");
    REIS(0,face,"found");
    REIS(REF_TRUE,reversed,"not rev");

    nodes[0]=3;nodes[1]=4;
    REIS(REF_NOT_FOUND,ref_front_find(ref_front,nodes,
				      &face,&reversed),"missing");
    REIS(REF_EMPTY,face,"found");

    RSS(ref_front_free(ref_front),"free");
  }

  return 0;
}
