#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_histogram.h"


int main( void )
{

  {
    REF_HISTOGRAM ref_histogram;
    REIS(REF_NULL, ref_histogram_free(NULL),"dont free NULL");
    RSS(ref_histogram_create(&ref_histogram),"create");
    REIS(10, ref_histogram_n(ref_histogram),"bins");
    RSS(ref_histogram_free(ref_histogram),"free");
  }

  return 0;
}
