#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_html.h"

int main( int argc, char *argv[] )
{

  if (2 == argc) 
    {
      REF_HTML ref_html;
      REF_DBL origin[3];
      REF_DBL sys[12];
      RSS(ref_html_create( &ref_html, argv[1] ), "open html" );

      origin[0] = 0.0;origin[1] = 0.0;origin[2] = 0.0;
      sys[0] = 2.0; sys[1] = 1.0; sys[2] = 1.0;
      sys[ 3] = 1.0; sys[ 4] = 0.0; sys[ 5] = 0.0;
      sys[ 6] = 0.0; sys[ 7] = 1.0; sys[ 8] = 0.0;
      sys[ 9] = 0.0; sys[10] = 0.0; sys[11] = 1.0;      
      RSS(ref_html_diagonal_system( ref_html, origin, sys ), "sys html" );

      RSS(ref_html_free( ref_html ), "close html" );
      return 0;
    }

  { /* html .meshb */
    REF_HTML ref_html;
    char file[] = "ref_html_test.html";
    RSS(ref_html_create( &ref_html, file ), "open html" );
    RSS(ref_html_free( ref_html ), "close html" );
    REIS(0, remove( file ), "test clean up");
  }

  return 0;
}
