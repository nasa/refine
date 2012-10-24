
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_histogram.h"
#include "ref_malloc.h"
#include "ref_edge.h"
#include "ref_mpi.h"

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

REF_STATUS ref_histogram_ratio( REF_GRID ref_grid )
{
  REF_HISTOGRAM ref_histogram;
  REF_EDGE ref_edge;
  REF_INT edge, part;
  REF_DBL ratio;

  RSS( ref_histogram_create(&ref_histogram),"create");
  RSS( ref_edge_create( &ref_edge, ref_grid ), "make edges" );

  for (edge=0;edge< ref_edge_n(ref_edge);edge++)
    {
      RSS( ref_edge_part( ref_edge, edge, &part ), "edge part");
      if ( part == ref_mpi_id )
	{
	  RSS( ref_node_ratio( ref_grid_node(ref_grid), 
			       ref_edge_e2n(ref_edge, 0, edge),
			       ref_edge_e2n(ref_edge, 1, edge), 
			       &ratio ), "rat");
	  RSS( ref_histogram_add( ref_histogram, ratio ), "add");
	}
    }

  RSS( ref_histogram_gather( ref_histogram ), "gather");

  RSS( ref_edge_free(ref_edge), "free edge" );
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_gather( REF_HISTOGRAM ref_histogram )
{
  REF_INT *bins;
  REF_INT i;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_malloc( bins, ref_histogram_n(ref_histogram), REF_INT );

  RSS( ref_mpi_sum( ref_histogram->bins, bins, 
		    ref_histogram_n(ref_histogram), REF_INT_TYPE ), "sum" );


  for ( i=0;i<ref_histogram_n(ref_histogram);i++ )
    if ( ref_mpi_master )
      {
	ref_histogram->bins[i] = bins[i];
      }
    else
      {
	ref_histogram->bins[i] = 0;
      }

  ref_free( bins );

  if ( ref_mpi_master ) printf("imp min max\n");

  return REF_SUCCESS;
}
