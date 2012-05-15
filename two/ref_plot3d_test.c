#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_plot3d.h"

int main( int argc, char *argv[] )
{
  REF_PLOT3D ref_plot3d;

  if (2 == argc) 
    {
      RSS( ref_plot3d_from_file( &ref_plot3d, argv[1] ), "from file" );
      RSS( ref_plot3d_tec( ref_plot3d, "ref_plot3d_test.tec" ), "tec" );
      RSS( ref_plot3d_free( ref_plot3d ), "free" );
      return 0;
    }

  return 0;
}
