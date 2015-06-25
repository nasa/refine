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
    as[0]=a0; as[1]=a1;as[2]=a2;

    RSS(ref_args_inspect( n, as ), "echo");
  }


  
  return 0;
}
