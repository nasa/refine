
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_args.h"

REF_STATUS ref_args_inspect( REF_INT n, char **args )
{
  REF_INT i;
  
  for (i=0;i<n;i++)
    printf("%d : '%s'\n",i,args[i]);

  return REF_SUCCESS;
}

REF_STATUS ref_args_find( REF_INT n, char **args, char *target, REF_INT *pos )
{
  REF_INT i;

  *pos = REF_EMPTY;
  
  for (i=0;i<n;i++)
    {
      if ( 0 == strcmp( target, args[i] ) )
	{
	  *pos = i;
	  return REF_SUCCESS;
	}
    }
  
  return REF_NOT_FOUND;
}
  
