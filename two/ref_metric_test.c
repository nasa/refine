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
    RSS(ref_metric_create(&ref_metric),"create");
    RSS(ref_metric_free(ref_metric),"free");
  }

  { /* set */
    REF_DBL m[6]={ 1.0, 0.0, 0.0, 
                        1.0, 0.0,
                             1.0};
    REF_INT i;
    RSS(ref_metric_create(&ref_metric),"create");

    RSS(ref_metric_set(ref_metric,0,m),"set");

    for (i=0;i<6;i++)
      RWDS( m[i], ref_metric_m(ref_metric,i,0), -1.0, "m" );

    RSS(ref_metric_free(ref_metric),"free");
  }

  { /* realloc */
    REF_DBL m[6]={ 1.0, 0.0, 0.0, 
                        1.0, 0.0,
                             1.0};
    REF_INT max;
    REF_INT i;
    RSS(ref_metric_create(&ref_metric),"create");

    max = ref_metric_max(ref_metric);

    RSS(ref_metric_set(ref_metric,max,m),"set");
    TAS(ref_metric_max(ref_metric)>max,"max increased");

    for (i=0;i<6;i++)
      RWDS( m[i], ref_metric_m(ref_metric,i,max), -1.0, "m" );

    RSS(ref_metric_free(ref_metric),"free");
  }

  return 0;
}
