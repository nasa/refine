
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
  int i,l,j,ir, indxt,q,n;

  if(0==length) return;
  for(i=0;i<length;i++) sortedIndex[i] = i;
  if(1==length) return;

  //printf("length: %d\n",length);
  //printf("arrin:  %d %d %d\n",arrayInput[0],arrayInput[1],arrayInput[2]);
  //printf("arrout: %d %d %d\n",sortedIndex[0],sortedIndex[1],sortedIndex[2]);

  n = length;

  if (n < 2) return;

  l=(n >> 1)+1; 
  ir=n;
  //printf("l ir: %d %d\n",l,ir);
  for (;;) {
    if (l > 1) {
      indxt=sortedIndex[(--l)-1]; 
      q=arrayInput[indxt];
    } else {
      indxt=sortedIndex[ir-1];  
      q=arrayInput[indxt];
      sortedIndex[ir-1]=sortedIndex[0];
      if (--ir == 1) {
	sortedIndex[0]=indxt;
	break; 
      } 
    }
    i=l;
    j=l+l;

    while (j <= ir) { 
      //printf("j ir: %d %d\n",j,ir);
      if ( j < ir ) {
	if (arrayInput[sortedIndex[j-1]] < arrayInput[sortedIndex[j]]) j++; 
      }
      if (q < arrayInput[sortedIndex[j-1]]) { 
	sortedIndex[i-1]=sortedIndex[j-1];
	i=j; 
	j <<= 1; 
      } else break; 
    }
    sortedIndex[i-1]=indxt;
  } 
}

int sortSearch( int length, int *arrayInput, int index  )
{
  return EMPTY;
}
