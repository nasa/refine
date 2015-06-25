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

  {
    REF_INT n = 3;
    char *a0 = "program";
    char *a1 = "-1";
    char *a2 = "-2";
    char *as[3];
    REF_INT i;
    as[0]=a0; as[1]=a1;as[2]=a2;
    
    for (i=0;i<n;i++)
      printf("%d : '%s'\n",i,as[i]);
  }


  
  return 0;
}
