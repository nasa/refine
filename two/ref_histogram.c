
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_histogram.h"
#include "ref_malloc.h"
#include "ref_edge.h"
#include "ref_mpi.h"

#include "ref_adapt.h"

REF_STATUS ref_histogram_create( REF_HISTOGRAM *ref_histogram_ptr )
{
  REF_HISTOGRAM ref_histogram;

  ref_malloc( *ref_histogram_ptr, 1, REF_HISTOGRAM_STRUCT );
  ref_histogram = (*ref_histogram_ptr);

  ref_histogram_n(ref_histogram) = 21;

  ref_malloc_init( ref_histogram->bins, 
		   ref_histogram_n(ref_histogram), REF_INT, 0 );

  ref_histogram_max(ref_histogram) = -1.0e20;
  ref_histogram_min(ref_histogram) =  1.0e20;

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
  REF_INT i;

  if ( observation <= 0.0 ) return REF_INVALID;

  ref_histogram_max(ref_histogram) = 
    MAX(ref_histogram_max(ref_histogram),observation);
  ref_histogram_min(ref_histogram) = 
    MIN(ref_histogram_min(ref_histogram),observation);

  i = ref_histogram_to_bin(observation);
  i = MIN(i,ref_histogram_n(ref_histogram)-1);
  i = MAX(i,0);

  ref_histogram_bin( ref_histogram, i )++;

  /*
    printf("%f:%f:%f\n",
    ref_histogram_to_obs(i),
    observation,
    ref_histogram_to_obs(i-1));
  */

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_ratio( REF_GRID ref_grid )
{
  REF_HISTOGRAM ref_histogram;
  REF_EDGE ref_edge;
  REF_INT edge, part;
  REF_DBL ratio;
  REF_BOOL active;

  RSS( ref_histogram_create(&ref_histogram),"create");
  RSS( ref_edge_create( &ref_edge, ref_grid ), "make edges" );

  for (edge=0;edge< ref_edge_n(ref_edge);edge++)
    {
      RSS( ref_edge_part( ref_edge, edge, &part ), "edge part");
      RSS( ref_node_edge_twod( ref_grid_node(ref_grid), 
			       ref_edge_e2n(ref_edge, 0, edge),
			       ref_edge_e2n(ref_edge, 1, edge), 
			       &active ), "twod edge");
      active = ( active || !ref_grid_twod(ref_grid) );
      if ( part == ref_mpi_id && active )
	{
	  RSS( ref_node_ratio( ref_grid_node(ref_grid), 
			       ref_edge_e2n(ref_edge, 0, edge),
			       ref_edge_e2n(ref_edge, 1, edge), 
			       &ratio ), "rat");
	  RSS( ref_histogram_add( ref_histogram, ratio ), "add");
	}
    }

  RSS( ref_histogram_gather( ref_histogram ), "gather");
  if ( ref_mpi_master ) RSS( ref_histogram_print( ref_histogram ), "print");

  RSS( ref_edge_free(ref_edge), "free edge" );
  RSS( ref_histogram_free(ref_histogram), "free gram" );
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_quality( REF_GRID ref_grid )
{
  REF_HISTOGRAM ref_histogram;
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality;

  RSS( ref_histogram_create(&ref_histogram),"create");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      if ( ref_node_part(ref_grid_node(ref_grid),nodes[0]) == ref_mpi_id )
	{
	  RSS( ref_node_tri_quality( ref_grid_node(ref_grid),
				     nodes,&quality ), "qual");
	  RSS( ref_histogram_add( ref_histogram, quality ), "add");
	}
    }

  RSS( ref_histogram_gather( ref_histogram ), "gather");
  if ( ref_mpi_master ) RSS( ref_histogram_print( ref_histogram ), "print");

  RSS( ref_histogram_free(ref_histogram), "free gram" );
  return REF_SUCCESS;
}

REF_STATUS ref_histogram_gather( REF_HISTOGRAM ref_histogram )
{
  REF_INT *bins;
  REF_DBL min, max;
  REF_INT i;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_malloc( bins, ref_histogram_n(ref_histogram), REF_INT );

  RSS( ref_mpi_sum( ref_histogram->bins, bins, 
		    ref_histogram_n(ref_histogram), REF_INT_TYPE ), "sum" );
  RSS( ref_mpi_max( &ref_histogram_max( ref_histogram ), &max, 
		    REF_DBL_TYPE ), "max" );
  RSS( ref_mpi_min( &ref_histogram_min( ref_histogram ), &min, 
		    REF_DBL_TYPE ), "min" );

  if ( ref_mpi_master )
    {
      ref_histogram_max( ref_histogram ) = max;
      ref_histogram_min( ref_histogram ) = min;
      for ( i=0;i<ref_histogram_n(ref_histogram);i++ )
	ref_histogram->bins[i] = bins[i];
    }
  else
    {
      ref_histogram_max( ref_histogram ) = -1.0e20;
      ref_histogram_min( ref_histogram ) =  1.0e20;
      for ( i=0;i<ref_histogram_n(ref_histogram);i++ )
	ref_histogram->bins[i] = 0;
    }

  ref_free( bins );

  return REF_SUCCESS;
}

REF_STATUS ref_histogram_print( REF_HISTOGRAM ref_histogram )
{
  REF_INT i, sum;

  sum = 0;
  for (i=0;i<ref_histogram_n(ref_histogram);i++)
    sum += ref_histogram_bin( ref_histogram, i );

  printf("%10.3f\n", ref_histogram_min( ref_histogram ));
  for (i=0;i<ref_histogram_n(ref_histogram);i++)
    if ( ( ref_histogram_to_obs(i)   > ref_adapt_split_ratio ||
	   ref_histogram_to_obs(i-1) < ref_adapt_collapse_ratio ) &&
	 ref_histogram_bin( ref_histogram, i ) > 0 )
      {
	printf("%2d:%7.3f:%10d *\n", i, 
	       ref_histogram_to_obs(i),ref_histogram_bin( ref_histogram, i ));
      }
    else
      {
	printf("%2d:%7.3f:%10d\n", i, 
	       ref_histogram_to_obs(i),ref_histogram_bin( ref_histogram, i ));
       }
  printf("%10.3f:%10d\n", ref_histogram_max( ref_histogram ), sum);

  return REF_SUCCESS;
}
