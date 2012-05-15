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
      return 0;
    }

  return 0;
}
