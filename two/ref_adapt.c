
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_adapt.h"
#include "ref_edge.h"

#include "ref_malloc.h"
#include "ref_mpi.h"

#include "ref_collapse.h"
#include "ref_split.h"

REF_STATUS ref_adapt_pass( REF_GRID ref_grid )
{

  RSS( ref_collapse_pass( ref_grid ), "col pass");
  RSS( ref_split_pass( ref_grid ), "split pass");

  return REF_SUCCESS;
}
