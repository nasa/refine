
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_adapt.h"
#include "ref_edge.h"

#include "ref_malloc.h"
#include "ref_mpi.h"

#include "ref_collapse.h"
#include "ref_split.h"
#include "ref_smooth.h"

#include "ref_gather.h"

REF_DBL ref_adapt_split_ratio = 1.5;
REF_DBL ref_adapt_split_quality_absolute = 1.0e-3;
REF_DBL ref_adapt_split_quality_relative = 0.6;

REF_DBL ref_adapt_collapse_ratio = 0.6;
REF_DBL ref_adapt_collapse_quality_absolute = 1.0e-3;
REF_DBL ref_adapt_collapse_ratio_limit = 3.0;

REF_DBL ref_adapt_smooth_min_quality = 1.0e-3;

REF_STATUS ref_adapt_pass( REF_GRID ref_grid )
{
  if (ref_grid_twod(ref_grid))
    {
      RSS( ref_adapt_twod_pass( ref_grid ), "pass");
    }
  else
    {
      RSS( ref_adapt_threed_pass( ref_grid ), "pass");
    }
  return REF_SUCCESS;
}

REF_STATUS ref_adapt_threed_pass( REF_GRID ref_grid )
{

  RSS( ref_geom_verify_topo( ref_grid ), "adapt preflight check");
  ref_gather_blocking_frame( ref_grid, "threed pass" );
  RSS( ref_collapse_pass( ref_grid ), "col pass");
  RSS( ref_geom_verify_topo( ref_grid ), "collapse geom typo check");
  ref_gather_blocking_frame( ref_grid, "collapse" );
  RSS( ref_split_pass( ref_grid ), "split pass");
  RSS( ref_geom_verify_topo( ref_grid ), "split geom typo check");
  ref_gather_blocking_frame( ref_grid, "split" );
  RSS( ref_smooth_threed_pass( ref_grid ), "smooth pass");
  ref_gather_blocking_frame( ref_grid, "smooth" );

  return REF_SUCCESS;
}

REF_STATUS ref_adapt_twod_pass( REF_GRID ref_grid )
{

  ref_gather_blocking_frame( ref_grid, "twod pass" );
  RSS( ref_collapse_twod_pass( ref_grid ), "col pass");
  ref_gather_blocking_frame( ref_grid, "collapse" );
  RSS( ref_split_twod_pass( ref_grid ), "split pass");
  ref_gather_blocking_frame( ref_grid, "split" );
  RSS( ref_smooth_twod_pass( ref_grid ), "smooth pass");
  ref_gather_blocking_frame( ref_grid, "smooth" );

  return REF_SUCCESS;
}
