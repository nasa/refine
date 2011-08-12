#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid_import.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  /* aflr3 */

  TSS(ref_grid_import_ugrid("/ump/fldmd/home/mikepark/FUN3D/GnuTestCase/Grids/gbumpf_MX.ugrid",&ref_grid),"from ugrid");

  TSS(ref_grid_free(ref_grid),"free");

  return 0;
}
