
/* Sort - for doing sorts and searches
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "sort.h"

void sortHeap( int length, int *arrayInput, int *sortedIndex  )
{
  int i, l, ir, indxt, q;
  unsigned int n, j;

  for(i=0;i<length;i++) sortedIndex[i] = i;

  if (length < 2) return;

  n = length;
  l=(n >> 1)+1; 
  ir=n-1;
  for (;;) {
    if (l > 1) {
      l--;
      indxt=sortedIndex[l-1]; 
      q=arrayInput[indxt];
    } else {
      indxt=sortedIndex[ir];  
      q=arrayInput[indxt];
      sortedIndex[ir]=sortedIndex[0];
      if (--ir == 0) {
	sortedIndex[0]=indxt;
	break; 
      } 
    }
    i=l-1;
    j=l+i;

    while (j <= ir) { 
      if ( j < ir ) {
	if (arrayInput[sortedIndex[j]] < arrayInput[sortedIndex[j+1]]) j++; 
      }
      if (q < arrayInput[sortedIndex[j]]) { 
	sortedIndex[i]=sortedIndex[j];
	i=j; 

	j++;
	j <<= 1; 
	j--;

      } else break; 
    }
    sortedIndex[i]=indxt;
  } 
}

int sortSearch( int length, int *sortednodes, int targetnode  )
{
  int lower, upper, mid;

  if (length<1) return EMPTY;

  if ( targetnode < sortednodes[0] || targetnode > sortednodes[length-1] ) 
    return EMPTY;

  lower = 0;
  upper = length-1;
  mid = length >> 1;

  if (targetnode==sortednodes[lower]) return lower;
  if (targetnode==sortednodes[upper]) return upper;

  while ( (lower < mid) && (mid < upper) ) {
    if ( targetnode >= sortednodes[mid] ) {
      if ( targetnode == sortednodes[mid] ) return mid;
      lower = mid;
    } else {
      upper = mid;
    }
    mid = (lower+upper) >> 1;
  }

  return EMPTY;
}
