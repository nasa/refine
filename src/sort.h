
/* Sort - for doing sorts and searches
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef SORT_H
#define SORT_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

void sortHeap( int length, int *arrayInput, int *sortedIndex );
void sortDoubleHeap( int length, double *arrayInput, int *sortedIndex );
int sortSearch( int length, int *arrayInput, int index );

END_C_DECLORATION

#endif /* QUEUE_H */
