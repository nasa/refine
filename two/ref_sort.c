
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

REF_STATUS ref_sort_unique( REF_INT n, REF_INT *original, 
			    REF_INT *nunique, REF_INT *unique )
{
  REF_INT i, j;

  *nunique = REF_EMPTY;

  RSS( ref_sort_insertion( n, original, unique ), "sort in unique");

  j=0;
  for(i=1;i<n;i++)
    {
      if ( unique[j] != unique[i] ) j++;
      if ( j != i ) unique[j] = unique[i];
    }

  *nunique = j+1;

  return REF_SUCCESS;
}

REF_STATUS ref_sort_search( REF_INT n, REF_INT *ascending_list, 
			    REF_INT target, REF_INT *position )
{
  int lower, upper, mid;

  *position = REF_EMPTY;

  if (n<1) return REF_INVALID;

  if ( target < ascending_list[0] || target > ascending_list[n-1] ) 
    return REF_NOT_FOUND;

  lower = 0;
  upper = n-1;
  mid = n >> 1; /* fast divide by two */

  if (target==ascending_list[lower]) 
    {
      *position = lower;
      return REF_SUCCESS;
    }
  if (target==ascending_list[upper])
    {
      *position = upper;
      return REF_SUCCESS;
    }

  while ( (lower < mid) && (mid < upper) ) {
    if ( target >= ascending_list[mid] ) {
      if ( target == ascending_list[mid] )
	{
	  *position = mid;
	  return REF_SUCCESS;
	}
      lower = mid;
    } else {
      upper = mid;
    }
    mid = (lower+upper) >> 1;
  }

  return REF_NOT_FOUND;
}
