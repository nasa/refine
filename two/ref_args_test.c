#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_args.h"

int main( int argc, char *argv[] )
{

  if ( argc > 1 )
    {
      REF_INT i;
      for (i=0;i<argc;i++)
	printf("%d : '%s'\n",i,argv[i]);
    }

  return 0;
}
