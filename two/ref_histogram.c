
#include <stdlib.h>
#include <stdio.h>

#include "ref_histogram.h"
#include "ref_malloc.h"

REF_STATUS ref_histogram_create( REF_HISTOGRAM *ref_histogram_ptr )
{
  REF_HISTOGRAM ref_histogram;

  ref_malloc( *ref_histogram_ptr, 1, REF_HISTOGRAM_STRUCT );
  ref_histogram = (*ref_histogram_ptr);

  ref_histogram_n(ref_histogram) = 10;

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_free( REF_HISTOGRAM ref_histogram )
{
  if ( NULL == (void *)ref_histogram ) return REF_NULL;
  ref_free( ref_histogram );
  return REF_SUCCESS;
}
