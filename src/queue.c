
/* Queue for storing grid transformations that can be shared to off-processors
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "queue.h"

Queue* queueCreate(  )
{
  Queue *queue;
  queue = malloc( sizeof(Queue) );
  queue->maxTransactions = 100;
  queue->removed = malloc( queue->maxTransactions * sizeof(int) );
  queue->added = malloc( queue->maxTransactions * sizeof(int) );
  queue->maxRemoved = 100;
  queue->removedNodes = malloc( 4 * queue->maxRemoved * sizeof(int) );
  queue->maxAdded = 100;
  queue->addedNodes = malloc( 5 * queue->maxAdded * sizeof(int) );
  queue->addedXYZs = malloc( 12 * queue->maxAdded * sizeof(double) );
  return queueReset(queue);
}

void queueFree( Queue *queue )
{
  free( queue->addedXYZs );
  free( queue->addedNodes );
  free( queue->removedNodes );
  free( queue->added );
  free( queue->removed );
  free( queue );
}

Queue *queueReset( Queue *queue )
{
  queue->transactions = 1;
  queue->nRemoved = 0;
  queue->added[0] = 0;
  queue->removed[0] = 0;
  return queue;
}

int queueTransactions( Queue *queue )
{
  return queue->transactions;
}

Queue *queueNewTransaction( Queue *queue )
{
  if (queue->transactions>=queue->maxTransactions) {
    queue->maxTransactions += 100;
    queue->removed = realloc( queue->removed, 
			      queue->maxTransactions * sizeof(int) );
    queue->added   = realloc( queue->added, 
			      queue->maxTransactions * sizeof(int) );
  }
  queue->removed[queue->transactions] = 0;
  queue->added[queue->transactions] = 0;
  queue->transactions++;
  return queue;
}

Queue *queueRemove( Queue *queue, int *nodes )
{
  int i;
  queue->removed[queue->transactions-1]++;
  if (queue->nRemoved>=queue->maxRemoved) {
    queue->maxRemoved += 100;
    queue->removedNodes = realloc( queue->removedNodes, 
				   4 * queue->maxRemoved * sizeof(int) );
  }
  for (i=0;i<4;i++) queue->removedNodes[i+4*queue->nRemoved] = nodes[i];
  queue->nRemoved++;
  return queue;
}

int queueRemoved( Queue *queue, int transaction )
{
  if ( transaction<0 || transaction>=queueTransactions(queue) ) return EMPTY;
  return queue->removed[transaction];
}

Queue *queueRemovedNodes( Queue *queue, int index, int *nodes )
{
  int i;
  if ( index<0 || index>queue->nRemoved ) return NULL;
  for (i=0;i<4;i++) nodes[i] = queue->removedNodes[i+4*index];
  return queue;
}

Queue *queueAdd( Queue *queue, int *nodes, double *xyzs )
{
  int i;
  queue->added[queue->transactions-1]++;
  if (queue->nAdded>=queue->maxAdded) {
    queue->maxAdded += 100;
    queue->addedNodes = realloc( queue->addedNodes, 
				 5 * queue->maxAdded * sizeof(int) );
    queue->addedXYZs  = realloc( queue->addedXYZs, 
				 12 * queue->maxAdded * sizeof(double) );
  }
  for (i=0;i<5 ;i++) queue->addedNodes[i+5*queue->nAdded] = nodes[i];
  for (i=0;i<12;i++) queue->addedXYZs[i+12*queue->nAdded] = xyzs[i];
  queue->nAdded++;
  return queue;
}

int queueAdded( Queue *queue, int transaction )
{
  if ( transaction<0 || transaction>=queueTransactions(queue) ) return EMPTY;
  return queue->added[transaction];
}

Queue *queueAddedNodes( Queue *queue, int index, int *nodes )
{
  int i;
  if ( index<0 || index>queue->nAdded ) return NULL;
  for (i=0;i<5;i++) nodes[i] = queue->addedNodes[i+5*index];
  return queue;
}

Queue *queueAddedXYZs( Queue *queue, int index, double *xyzs )
{
  int i;
  if ( index<0 || index>queue->nAdded ) return NULL;
  for (i=0;i<12;i++) xyzs[i] = queue->addedXYZs[i+12*index];
  return queue;
}
