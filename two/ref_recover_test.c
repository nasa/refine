#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_recover.h"

int main( void )
{

  { /* init */
    REF_RECOVER ref_recover;
    REIS(REF_NULL, ref_recover_free(NULL),"dont free NULL");
    RSS(ref_recover_create(&ref_recover),"create");
    REIS( 0, ref_recover_n(ref_recover), "init no recover");
    RSS(ref_recover_free(ref_recover),"free");
  }

  return 0;
}
