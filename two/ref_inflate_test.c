#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_inflate.h"

#include "ref_math.h"

#include "ref_import.h"
#include "ref_export.h"

#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_sort.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_dict.h"
#include "ref_list.h"

#include "ref_edge.h"

int main(  )
{

  {
    REF_DBL H, h0, rate;
    REF_INT n;
    H = 10.0;
    h0 = 1.0;
    n = 10;
    RSS( ref_inflate_rate(n,h0,H,&rate),"rate");
    RWDS( 1.0, rate, -1.0, "uniform" );
  }

  {
    REF_DBL H, h0, rate;
    REF_INT n;
    h0 = 0.5;
    n = 10;
    rate = 1.1469;
    RSS( ref_inflate_total_thickness(n,h0,rate,&H),"total");
    RWDS( 10.0, H, 1.0e-3, "not total" );
  }

  return 0;
}
