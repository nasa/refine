
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_histogram.h"
#include "ref_malloc.h"

REF_STATUS ref_histogram_create( REF_HISTOGRAM *ref_histogram_ptr )
{
  REF_HISTOGRAM ref_histogram;

  ref_malloc( *ref_histogram_ptr, 1, REF_HISTOGRAM_STRUCT );
  ref_histogram = (*ref_histogram_ptr);

  ref_histogram_n(ref_histogram) = 10;

  ref_malloc_init( ref_histogram->bins, 
		   ref_histogram_n(ref_histogram), REF_INT, 0 );

  ref_histogram_max(ref_histogram) = -20.0;
  ref_histogram_min(ref_histogram) =  20.0;

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_free( REF_HISTOGRAM ref_histogram )
{
  if ( NULL == (void *)ref_histogram ) return REF_NULL;
  ref_free( ref_histogram->bins );
  ref_free( ref_histogram );
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_add( REF_HISTOGRAM ref_histogram, REF_DBL observation )
{
  REF_DBL r;
  REF_INT i;

  if ( observation <= 0.0 ) return REF_INVALID;

  r = log10(observation);
  ref_histogram_max(ref_histogram) = MAX(ref_histogram_max(ref_histogram),r);
  ref_histogram_min(ref_histogram) = MIN(ref_histogram_min(ref_histogram),r);

  i = (REF_INT)r;
  i = i+ref_histogram_n(ref_histogram)/2;
  i = MIN(i,ref_histogram_n(ref_histogram)-1);
  i = MAX(i,0);

  ref_histogram_bin( ref_histogram, i )++;

  return REF_SUCCESS;
}

