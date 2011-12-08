
#include <stdlib.h>
#include <stdio.h>

#include "ref_sort.h"

REF_STATUS ref_sort_insertion( REF_INT n, REF_INT *original, REF_INT *sorted )
{
  REF_INT i, j, smallest, temp;

  for(i=0;i<n;i++)
    sorted[i] = original[i];

  for(i=0;i<n;i++)
    {
      smallest = i;
      for(j=i+1;j<n;j++)
	{
	  if ( sorted[j] < sorted[smallest] )
	    smallest = j;
	}
      temp = sorted[i];
      sorted[i] = sorted[smallest];
      sorted[smallest] = temp;
    }

  return REF_SUCCESS;
}
