#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_front.h"

int main( void )
{

  {
    REF_FRONT ref_front;
    REIS(REF_NULL, ref_front_free(NULL),"dont free NULL");
    RSS(ref_front_create(&ref_front,2),"create");
    REIS( 0, ref_front_n(ref_front), "init no front");
    REIS( 2, ref_front_node_per(ref_front), "init per");
    RSS(ref_front_free(ref_front),"free");
  }

  return 0;
}
