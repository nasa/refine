
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

Queue* queueCreate( int nodeSize )
{
  Queue *queue;
  queue = malloc( sizeof(Queue) );
  queue->nodeSize = nodeSize;
  queue->maxTransactions = 100;
  queue->removedCells = malloc( queue->maxTransactions * sizeof(int) );
  queue->addedCells = malloc( queue->maxTransactions * sizeof(int) );
  queue->removedFaces = malloc( queue->maxTransactions * sizeof(int) );
  queue->addedFaces = malloc( queue->maxTransactions * sizeof(int) );
  queue->maxRemovedCells = 100;
  queue->removedCellNodes = malloc( 4 * queue->maxRemovedCells * sizeof(int) );
  queue->maxAddedCells = 100;
  queue->addedCellNodes = malloc( 9 * queue->maxAddedCells * sizeof(int) );
  queue->addedCellXYZs = malloc( 4*queue->nodeSize * queue->maxAddedCells 
				 * sizeof(double) );
  queue->maxRemovedFaces = 100;
  queue->removedFaceNodes = malloc( 3 * queue->maxRemovedFaces * sizeof(int) );
  queue->maxAddedFaces = 100;
  queue->addedFaceNodes = malloc( 4 * queue->maxAddedFaces * sizeof(int) );
  queue->addedFaceUVs = malloc( 6 * queue->maxAddedFaces * sizeof(double) );
  return queueReset(queue);
}

void queueFree( Queue *queue )
{
  if (NULL==queue) return;
  free( queue->addedFaceUVs );
  free( queue->addedFaceNodes );
  free( queue->removedFaceNodes );
  free( queue->addedFaces );
  free( queue->removedFaces );
  free( queue->addedCellXYZs );
  free( queue->addedCellNodes );
  free( queue->removedCellNodes );
  free( queue->addedCells );
  free( queue->removedCells );
  free( queue );
}

Queue *queueReset( Queue *queue )
{
  if (NULL==queue) return NULL;
  queue->transactions = 1;
  queue->nAddedCells = 0;
  queue->nRemovedCells = 0;
  queue->addedCells[0] = 0;
  queue->removedCells[0] = 0;
  queue->nAddedFaces = 0;
  queue->nRemovedFaces = 0;
  queue->addedFaces[0] = 0;
  queue->removedFaces[0] = 0;
  return queue;
}

int queueNodeSize( Queue *queue )
{
  if (NULL==queue) return EMPTY;
  return queue->nodeSize;
}

int queueTransactions( Queue *queue )
{
  if (NULL==queue) return EMPTY;
  return queue->transactions;
}

Queue *queueNewTransaction( Queue *queue )
{
  if (NULL==queue) return NULL;
  if (queue->transactions>=queue->maxTransactions) {
    queue->maxTransactions += 100;
    queue->removedCells = realloc( queue->removedCells, 
				   queue->maxTransactions * sizeof(int) );
    queue->addedCells   = realloc( queue->addedCells, 
				   queue->maxTransactions * sizeof(int) );
    queue->removedFaces = realloc( queue->removedFaces, 
				   queue->maxTransactions * sizeof(int) );
    queue->addedFaces   = realloc( queue->addedFaces, 
				   queue->maxTransactions * sizeof(int) );
  }
  queue->removedCells[queue->transactions] = 0;
  queue->addedCells[queue->transactions] = 0;
  queue->removedFaces[queue->transactions] = 0;
  queue->addedFaces[queue->transactions] = 0;
  queue->transactions++;
  return queue;
}

Queue *queueRemoveCell( Queue *queue, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  queue->removedCells[queue->transactions-1]++;
  if (queue->nRemovedCells>=queue->maxRemovedCells) {
    queue->maxRemovedCells += 100;
    queue->removedCellNodes = realloc( queue->removedCellNodes, 
				       4 * queue->maxRemovedCells* sizeof(int));
  }
  for (i=0;i<4;i++) queue->removedCellNodes[i+4*queue->nRemovedCells]= nodes[i];
  queue->nRemovedCells++;
  return queue;
}

int queueRemovedCells( Queue *queue, int transaction )
{
  if (NULL==queue) return EMPTY;
  if ( transaction<0 || transaction>=queueTransactions(queue) ) return EMPTY;
  return queue->removedCells[transaction];
}

Queue *queueRemovedCellNodes( Queue *queue, int index, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nRemovedCells ) return NULL;
  for (i=0;i<4;i++) nodes[i] = queue->removedCellNodes[i+4*index];
  return queue;
}

Queue *queueAddCell( Queue *queue, int *nodes, double *xyzs )
{
  int i;
  if (NULL==queue) return NULL;
  queue->addedCells[queue->transactions-1]++;
  if (queue->nAddedCells>=queue->maxAddedCells) {
    queue->maxAddedCells += 100;
    queue->addedCellNodes = realloc( queue->addedCellNodes, 
				     9 * queue->maxAddedCells * sizeof(int) );
    queue->addedCellXYZs  = realloc( queue->addedCellXYZs, 
				     4*queue->nodeSize * queue->maxAddedCells 
				     * sizeof(double));
  }
  for (i=0;i<9 ;i++) queue->addedCellNodes[i+9*queue->nAddedCells] = nodes[i];
  for (i=0;i<4*queue->nodeSize;i++) 
    queue->addedCellXYZs[i+4*queue->nodeSize*queue->nAddedCells] = xyzs[i];
  queue->nAddedCells++;
  return queue;
}

int queueAddedCells( Queue *queue, int transaction )
{
  if (NULL==queue) return EMPTY;
  if ( transaction<0 || transaction>=queueTransactions(queue) ) return EMPTY;
  return queue->addedCells[transaction];
}

Queue *queueAddedCellNodes( Queue *queue, int index, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedCells ) return NULL;
  for (i=0;i<9;i++) nodes[i] = queue->addedCellNodes[i+9*index];
  return queue;
}

Queue *queueAddedCellXYZs( Queue *queue, int index, double *xyzs )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedCells ) return NULL;
  for (i=0;i<4*queue->nodeSize;i++) 
    xyzs[i] = queue->addedCellXYZs[i+4*queue->nodeSize*index];
  return queue;
}

int queueTotalRemovedCells( Queue *queue )
{
  if (NULL==queue) return EMPTY;
  return queue->nRemovedCells;
}

Queue *queueRemoveFace( Queue *queue, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  queue->removedFaces[queue->transactions-1]++;
  if (queue->nRemovedFaces>=queue->maxRemovedFaces) {
    queue->maxRemovedFaces += 100;
    queue->removedFaceNodes = realloc( queue->removedFaceNodes, 
				       3 * queue->maxRemovedFaces* sizeof(int));
  }
  for (i=0;i<3;i++) queue->removedFaceNodes[i+3*queue->nRemovedFaces]= nodes[i];
  queue->nRemovedFaces++;
  return queue;
}

int queueRemovedFaces( Queue *queue, int transaction )
{
  if (NULL==queue) return EMPTY;
  if ( transaction<0 || transaction>=queueTransactions(queue) ) return EMPTY;
  return queue->removedFaces[transaction];
}

Queue *queueRemovedFaceNodes( Queue *queue, int index, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nRemovedFaces ) return NULL;
  for (i=0;i<3;i++) nodes[i] = queue->removedFaceNodes[i+3*index];
  return queue;
}

Queue *queueAddFace( Queue *queue, int *nodes, double *uvs )
{
  int i;
  if (NULL==queue) return NULL;
  queue->addedFaces[queue->transactions-1]++;
  if (queue->nAddedFaces>=queue->maxAddedFaces) {
    queue->maxAddedFaces += 100;
    queue->addedFaceNodes = realloc( queue->addedFaceNodes, 
				     4 * queue->maxAddedFaces * sizeof(int) );
    queue->addedFaceUVs   = realloc( queue->addedFaceUVs, 
				     6 * queue->maxAddedFaces* sizeof(double));
  }
  for (i=0;i<4;i++) queue->addedFaceNodes[i+4*queue->nAddedFaces] = nodes[i];
  for (i=0;i<6;i++) queue->addedFaceUVs[i+6*queue->nAddedFaces] = uvs[i];
  queue->nAddedFaces++;
  return queue;
}

Queue *queueAddFaceScalar( Queue *queue, 
			   int n0, double u0, double v0,
			   int n1, double u1, double v1,
			   int n2, double u2, double v2, int faceId)
{
  int nodes[4];
  double uv[6];
  nodes[0] = n0;
  nodes[1] = n1;
  nodes[2] = n2;
  nodes[3] = faceId;
  uv[0] = u0;
  uv[1] = v0;
  uv[2] = u1;
  uv[3] = v1;
  uv[4] = u2;
  uv[5] = v2;
  return queueAddFace( queue, nodes, uv );
}


int queueAddedFaces( Queue *queue, int transaction )
{
  if (NULL==queue) return EMPTY;
  if ( transaction<0 || transaction>=queueTransactions(queue) ) return EMPTY;
  return queue->addedFaces[transaction];
}

Queue *queueAddedFaceNodes( Queue *queue, int index, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedFaces ) return NULL;
  for (i=0;i<4;i++) nodes[i] = queue->addedFaceNodes[i+4*index];
  return queue;
}

Queue *queueAddedFaceUVs( Queue *queue, int index, double *xyzs )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedFaces ) return NULL;
  for (i=0;i<6;i++) xyzs[i] = queue->addedFaceUVs[i+6*index];
  return queue;
}

Queue *queueDumpSize( Queue *queue, int *nInt, int *nDouble )
{
  *nInt = EMPTY;
  *nDouble = EMPTY;
  if (NULL==queue) return NULL;
  *nInt 
    = 6
    + 4 * queue->transactions
    + 4 * queue->nRemovedCells
    + 9 * queue->nAddedCells
    + 3 * queue->nRemovedFaces
    + 4 * queue->nAddedFaces;
  *nDouble = 4*queue->nodeSize * queue->nAddedCells + 6 * queue->nAddedFaces ;
  return queue;
}

Queue *queueDump( Queue *queue, int *ints, double *doubles )
{
  int i, d, node, size, transaction, removed, added;
  if (NULL==queue) return NULL;
  ints[0] = queue->nodeSize;
  ints[1] = queue->transactions;
  ints[2] = queue->nRemovedCells;
  ints[3] = queue->nAddedCells;
  ints[4] = queue->nRemovedFaces;
  ints[5] = queue->nAddedFaces;
  i = 6;
  d = 0;

  for(transaction=0;transaction<queue->transactions;transaction++){
    ints[i] = queue->removedCells[transaction]; i++;
  }
  for(removed=0;removed<queue->nRemovedCells;removed++){
    size = 4;
    for (node=0;node<size;node++) { 
      ints[i] = queue->removedCellNodes[node+size*removed]; i++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    ints[i] = queue->addedCells[transaction]; i++;
  }
  for(added=0;added<queue->nAddedCells;added++){
    size = 9;
    for (node=0;node<size;node++) { 
      ints[i] = queue->addedCellNodes[node+size*added]; i++;
    }
    size = 4*queue->nodeSize;
    for (node=0;node<size;node++) { 
      doubles[d] = queue->addedCellXYZs[node+size*added]; d++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    ints[i] = queue->removedFaces[transaction]; i++;
  }
  for(removed=0;removed<queue->nRemovedFaces;removed++){
    size = 3;
    for (node=0;node<size;node++) { 
      ints[i] = queue->removedFaceNodes[node+size*removed]; i++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    ints[i] = queue->addedFaces[transaction]; i++;
  }
  for(added=0;added<queue->nAddedFaces;added++){
    size = 4;
    for (node=0;node<size;node++) { 
      ints[i] = queue->addedFaceNodes[node+size*added]; i++;
    }
    size = 6;
    for (node=0;node<size;node++) { 
      doubles[d] = queue->addedFaceUVs[node+size*added]; d++;
    }
  }

  return queue;
}

Queue *queueLoad( Queue *queue, int *ints, double *doubles )
{
  int i, d, node, size, transaction, removed, added;

  if (NULL==queue) return NULL;

  if ( queue->nodeSize != ints[0] ) {
    printf("ERROR: %s: %d: queueLoad: incompatable nodeSize: %d %d\n",
	   __FILE__, __LINE__, queue->nodeSize, ints[0] );
    return NULL;
  }			    
  queue->transactions  = ints[1];
  queue->nRemovedCells = ints[2];
  queue->nAddedCells   = ints[3];
  queue->nRemovedFaces = ints[4];
  queue->nAddedFaces   = ints[5];
  i = 6;
  d = 0;

  if (queue->transactions > queue->maxTransactions) {
    queue->maxTransactions = queue->transactions;
    queue->removedCells = realloc( queue->removedCells, 
				   queue->maxTransactions * sizeof(int) );
    queue->addedCells   = realloc( queue->addedCells, 
				   queue->maxTransactions * sizeof(int) );
    queue->removedFaces = realloc( queue->removedFaces, 
				   queue->maxTransactions * sizeof(int) );
    queue->addedFaces   = realloc( queue->addedFaces, 
				   queue->maxTransactions * sizeof(int) );
  }
  if (queue->nRemovedCells > queue->maxRemovedCells) {
    queue->maxRemovedCells = queue->nRemovedCells;
    queue->removedCellNodes = realloc( queue->removedCellNodes, 
				       4 * queue->maxRemovedCells* sizeof(int));
  }

  if (queue->nAddedCells > queue->maxAddedCells) {
    queue->maxAddedCells = queue->nAddedCells;
    queue->addedCellNodes = realloc( queue->addedCellNodes, 
				     9 * queue->maxAddedCells * sizeof(int) );
    queue->addedCellXYZs  = realloc( queue->addedCellXYZs, 
				     4*queue->nodeSize * queue->maxAddedCells
				     * sizeof(double));
  }
  if (queue->nRemovedFaces > queue->maxRemovedFaces) {
    queue->maxRemovedFaces = queue->nRemovedFaces;
    queue->removedFaceNodes = realloc( queue->removedFaceNodes, 
				       3 * queue->maxRemovedFaces* sizeof(int));
  }
  if (queue->nAddedFaces > queue->maxAddedFaces) {
    queue->maxAddedFaces = queue->nAddedFaces;
    queue->addedFaceNodes = realloc( queue->addedFaceNodes, 
				     4 * queue->maxAddedFaces * sizeof(int) );
    queue->addedFaceUVs   = realloc( queue->addedFaceUVs, 
				     6 * queue->maxAddedFaces* sizeof(double));
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    queue->removedCells[transaction] = ints[i]; i++;
  }
  for(removed=0;removed<queue->nRemovedCells;removed++){
    size = 4;
    for (node=0;node<size;node++) { 
      queue->removedCellNodes[node+size*removed] = ints[i]; i++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    queue->addedCells[transaction] = ints[i]; i++;
  }
  for(added=0;added<queue->nAddedCells;added++){
    size = 9;
    for (node=0;node<size;node++) { 
      queue->addedCellNodes[node+size*added] = ints[i]; i++;
    }
    size = 4*queue->nodeSize;
    for (node=0;node<size;node++) { 
      queue->addedCellXYZs[node+size*added] = doubles[d]; d++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    queue->removedFaces[transaction] = ints[i]; i++;
  }
  for(removed=0;removed<queue->nRemovedFaces;removed++){
    size = 3;
    for (node=0;node<size;node++) { 
      queue->removedFaceNodes[node+size*removed] = ints[i]; i++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    queue->addedFaces[transaction] = ints[i]; i++;
  }
  for(added=0;added<queue->nAddedFaces;added++){
    size = 4;
    for (node=0;node<size;node++) { 
      queue->addedFaceNodes[node+size*added] = ints[i]; i++;
    }
    size = 6;
    for (node=0;node<size;node++) { 
      queue->addedFaceUVs[node+size*added] = doubles[d]; d++;
    }
  }

  return queue;
}

Queue *queueContents(Queue *queue, FILE *f)
{
  int transaction, removed, removedcell;
  int globalnodes[9];

  if (NULL==queue) return NULL;

  fprintf(f,"transactions %d\n",queue->transactions);
  fprintf(f,"total removed cells %d\n",queue->nRemovedCells);

  removedcell = 0;
  for (transaction=0;transaction<queueTransactions(queue);transaction++){
    fprintf(f,"transaction %d has %d removed cells\n",
	    transaction,queueRemovedCells(queue,transaction));
    for (removed=0;removed<queueRemovedCells(queue,transaction);removed++) {
      queueRemovedCellNodes( queue, removedcell, globalnodes );
      fprintf(f,"cell %d %d: %d %d %d %d\n",removed, removedcell,
	      globalnodes[0], globalnodes[1], globalnodes[2], globalnodes[3]);
      removedcell++;
    }
  }
  return queue;  
}
