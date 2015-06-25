
#include <stdlib.h>
#include <stdio.h>

#include "ref_args.h"

REF_STATUS ref_args_inspect( REF_INT n, char **args )
{
  REF_INT i;
  
  for (i=0;i<n;i++)
    printf("%d : '%s'\n",i,args[i]);

  return REF_SUCCESS;
}
