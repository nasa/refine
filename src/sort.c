
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
#include "sort.h"

void sortHeap( int length, int *arrayInput, int *sortedIndex  )
{
  int i,l,j,ir, indxt,q;

  if(0==length) return;
  for(i=0;i<length;i++) sortedIndex[i] = i;
  if(1==length) return;

  l=(length >> 1)+1; 
  ir=length;
  for (;;) {
    if (l > 1) {
      indxt=sortedIndex[--l]; 
      q=arrayInput[indxt];
    } else {
      indxt=sortedIndex[ir];  
      q=arrayInput[indxt];
      sortedIndex[ir]=sortedIndex[1];
      if (--ir == 1) {
	sortedIndex[1]=indxt;
	break; 
      } 
    }
    i=l;
    j=l+l;

    while (j <= ir) { 
      if ( j < ir && 
	   arrayInput[sortedIndex[j]] < arrayInput[sortedIndex[j+1]]) j++; 
      if (q < arrayInput[sortedIndex[j]]) { 
	sortedIndex[i]=sortedIndex[j];
	i=j; 
	j <<= 1; 
      } else break; 
    }
    sortedIndex[i]=indxt;
  } 
}

int sortSearch( int length, int *arrayInput, int index  )
{
  return EMPTY;
}
