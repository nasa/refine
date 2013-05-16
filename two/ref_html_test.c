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
      RSS(ref_html_create( &ref_html, argv[1] ), "open html" );
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
