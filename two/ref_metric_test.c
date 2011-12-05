#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_metric.h"
#include "ref_test.h"

int main( void )
{
  REF_METRIC ref_metric;

  {
    TFS(ref_metric_free(NULL),"dont free NULL");
    TSS(ref_metric_create(&ref_metric),"create");
    TSS(ref_metric_free(ref_metric),"free");
  }

  { /* set */
    REF_DBL m[6]={ 1.0, 0.0, 0.0, 
                        1.0, 0.0,
                             1.0};
    REF_INT i;
    TSS(ref_metric_create(&ref_metric),"create");

    TSS(ref_metric_set(ref_metric,0,m),"set");

    for (i=0;i<6;i++)
      TWDS( m[i], ref_metric_m(ref_metric,i,0), -1.0, "m" );

    TSS(ref_metric_free(ref_metric),"free");
  }

  SKIP_TEST("realloc")
  {
  }

  return 0;
}
