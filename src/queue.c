
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
  int i;
  queue = malloc( sizeof(Queue) );
  queue->transactions = 0;
  queue->maxTransactions = 10;
  queue->transactionNodes = malloc(queue->maxTransactions * sizeof(int) );
  for (i=0;i<queue->maxTransactions;i++) queue->transactionNodes[i]=0;
  queue->addedNodes = 0;
  queue->maxAddedNode = 10;
  queue->addedNode = malloc(queue->maxAddedNode * sizeof(int) );
  for (i=0;i<queue->maxAddedNode;i++) queue->addedNode[i]=0;
  return queue;
}

void queueFree( Queue *queue )
{
  free( queue->transactionNodes );
  free( queue );
}

Queue *queueReset( Queue *queue )
{
  queue->transactions = 0;
  return queue;
}

int queueTransactions( Queue *queue )
{
  return queue->transactions;
}

Queue *queueNewTransaction( Queue *queue )
{
  queue->transactions++;
  return queue;
}

int queueTransactionNodes( Queue *queue, int transaction )
{
  if ( transaction<0 || transaction>queueTransactions(queue) ) return EMPTY;
  return queue->transactionNodes[transaction];
}

Queue *queueAddNode( Queue *queue, int node )
{
  queue->transactionNodes[queue->transactions]++;
  queue->addedNode[queue->addedNodes] = node;
  queue->addedNodes++;
  return queue;
}

int queueAddedNode( Queue *queue, int index )
{
  if ( index<0 || index>=queue->addedNodes ) return EMPTY;
  return queue->addedNode[index];
}
