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

  { /* add face */
    REF_FRONT ref_front;
    REF_INT nodes[2];

    RSS(ref_front_create(&ref_front,2),"create");
    nodes[0]=1;nodes[1]=2;
    RSS(ref_front_insert(ref_front,nodes),"insert");
    REIS( 1, ref_front_n(ref_front), "init no front");

    RSS(ref_front_free(ref_front),"free");
  }

  { /* add faces, fore realloc */
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

    RSS(ref_front_free(ref_front),"free");
  }

  return 0;
}
