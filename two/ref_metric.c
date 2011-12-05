
#include <stdlib.h>
#include <stdio.h>

#include "ref_metric.h"

REF_STATUS ref_metric_create( REF_METRIC *ref_metric_ptr )
{
  REF_METRIC ref_metric;

  (*ref_metric_ptr) = NULL;
  (*ref_metric_ptr) = (REF_METRIC)malloc( sizeof(REF_METRIC_STRUCT) );
  RNS(*ref_metric_ptr,"malloc ref_metric NULL");
  ref_metric = (*ref_metric_ptr);

  ref_metric_max(ref_metric) = 10;

  ref_metric->m = (REF_DBL *)malloc( 6*ref_metric_max(ref_metric) * 
				     sizeof(REF_DBL) );
  RNS(ref_metric->m,"malloc ref_metric->m NULL");

  return REF_SUCCESS;
}

REF_STATUS ref_metric_free( REF_METRIC ref_metric )
{
  if ( NULL == (void *)ref_metric ) return REF_NULL;
  ref_cond_free( ref_metric->m );
  ref_cond_free( ref_metric );
  return REF_SUCCESS;
}

REF_STATUS ref_metric_set( REF_METRIC ref_metric, REF_INT node, REF_DBL *m )
{
  REF_INT entry;
  REF_INT orig, chunk;

  if ( node < 0 ) return REF_INVALID;

  if ( node >= ref_metric_max(ref_metric) )
    {
      orig = ref_metric_max( ref_metric );
      chunk = 100 + MAX( 0, node-orig );
      ref_metric_max(ref_metric) = orig + chunk;
      ref_metric->m = (REF_DBL *)realloc( ref_metric->m,
					  6*ref_metric_max(ref_metric) * 
					  sizeof(REF_DBL) );
      RNS(ref_metric->m,"realloc ref_metric->m NULL");
    }

  for( entry=0;entry<6;entry++ )
    ref_metric_m(ref_metric,entry,node) = m[entry];

  return REF_SUCCESS;
}
