
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
  int *transactionNodes;
};

Queue *queueCreate(  );
void queueFree( Queue * );
Queue *queueReset( Queue * );
int queueTransactions( Queue * );
Queue *queueNewTransaction( Queue * );
int queueTransactionNodes( Queue *, int transaction );
Queue *queueAddNode( Queue *, int node );

END_C_DECLORATION

#endif /* QUEUE_H */
