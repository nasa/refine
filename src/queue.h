
/* Queue for storing grid transformations that can be shared to off-processors
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef QUEUE_H
#define QUEUE_H

#include "master_header.h"

BEGIN_C_DECLORATION

typedef struct Queue Queue;

struct Queue {
  int transactions;
  int maxTransactions;
  int *removed;
  int maxRemoved;
  int nRemoved;
  int *removedNodes;
  int *added;
  int maxAdded;
  int nAdded;
  int *addedNodes;
  double *addedXYZs;
};

Queue *queueCreate(  );
void queueFree( Queue * );
Queue *queueReset( Queue * );
int queueTransactions( Queue * );
Queue *queueNewTransaction( Queue * );
Queue *queueRemoveCell( Queue *, int *nodes );
int queueRemovedCells( Queue *, int transaction );
Queue *queueRemovedCellNodes( Queue *, int index, int *nodes );
Queue *queueAddCell( Queue *, int *nodes, double *xyzs );
int queueAddedCells( Queue *, int transaction );
Queue *queueAddedCellNodes( Queue *, int index, int *nodes );
Queue *queueAddedCellXYZs( Queue *, int index, double *xyz );

END_C_DECLORATION

#endif /* QUEUE_H */
