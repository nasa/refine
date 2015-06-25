#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_args.h"

int main( int argc, char *argv[] )
{

  if ( argc > 1 )
    {
      RSS(ref_args_inspect( argc, argv ), "echo");
    }

  {
    REF_INT n = 3;
    char *a0 = "program";
    char *a1 = "-1";
    char *a2 = "-2";
    char *as[3];
    REF_INT pos;
    as[0]=a0; as[1]=a1;as[2]=a2;

    RSS(ref_args_find( n, as, "-1", &pos ), "echo");
    REIS( 1, pos, "location" );
    
    REIS(REF_NOT_FOUND, ref_args_find( n, as, "-h", &pos ), "not found");
    REIS( REF_EMPTY, pos, "location" );

    REIS(REF_NOT_FOUND, ref_args_find( n, as, "--long", &pos ), "not found");
    REIS( REF_EMPTY, pos, "location" );
  }

  return 0;
}
