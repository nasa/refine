
/* Queue for storing grid transformations that can be shared to off-processors
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  


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

  queue->maxRemovedCells = 100;
  queue->removedCellNodes = malloc( 8 * queue->maxRemovedCells * sizeof(int) );
  queue->maxAddedCells = 100;
  queue->addedCellNodes = malloc( 9 * queue->maxAddedCells * sizeof(int) );
  queue->addedCellXYZs = malloc( 4*queue->nodeSize * queue->maxAddedCells 
				 * sizeof(double) );

  queue->removedFaces = malloc( queue->maxTransactions * sizeof(int) );
  queue->addedFaces = malloc( queue->maxTransactions * sizeof(int) );

  queue->maxRemovedFaces = 100;
  queue->removedFaceNodes = malloc( 6 * queue->maxRemovedFaces * sizeof(int) );
  queue->maxAddedFaces = 100;
  queue->addedFaceNodes = malloc( 7 * queue->maxAddedFaces * sizeof(int) );
  queue->addedFaceUVs = malloc( 6 * queue->maxAddedFaces * sizeof(double) );

  queue->removedEdges = malloc( queue->maxTransactions * sizeof(int) );
  queue->addedEdges = malloc( queue->maxTransactions * sizeof(int) );

  queue->maxRemovedEdges = 100;
  queue->removedEdgeNodes = malloc( 4 * queue->maxRemovedEdges * sizeof(int) );
  queue->maxAddedEdges = 100;
  queue->addedEdgeNodes = malloc( 5 * queue->maxAddedEdges * sizeof(int) );
  queue->addedEdgeTs = malloc( 2 * queue->maxAddedEdges * sizeof(double) );

  return queueReset(queue);
}

void queueFree( Queue *queue )
{
  if (NULL==queue) return;
  free( queue->addedEdgeTs );
  free( queue->addedEdgeNodes );
  free( queue->removedEdgeNodes );
  free( queue->addedEdges );
  free( queue->removedEdges );

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

  queue->nAddedEdges = 0;
  queue->nRemovedEdges = 0;
  queue->addedEdges[0] = 0;
  queue->removedEdges[0] = 0;

  return queue;
}

Queue *queueResetCurrentTransaction( Queue *queue )
{
  if (NULL==queue) return NULL;

  queue->nAddedCells -= queue->addedCells[queue->transactions-1];
  queue->nRemovedCells -= queue->removedCells[queue->transactions-1];
  queue->addedCells[queue->transactions-1] = 0;
  queue->removedCells[queue->transactions-1] = 0;

  queue->nAddedFaces -= queue->addedFaces[queue->transactions-1];
  queue->nRemovedFaces -= queue->removedFaces[queue->transactions-1];
  queue->addedFaces[queue->transactions-1] = 0;
  queue->removedFaces[queue->transactions-1] = 0;

  queue->nAddedEdges -= queue->addedEdges[queue->transactions-1];
  queue->nRemovedEdges -= queue->removedEdges[queue->transactions-1];
  queue->addedEdges[queue->transactions-1] = 0;
  queue->removedEdges[queue->transactions-1] = 0;

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
    queue->removedEdges = realloc( queue->removedEdges, 
				   queue->maxTransactions * sizeof(int) );
    queue->addedEdges   = realloc( queue->addedEdges, 
				   queue->maxTransactions * sizeof(int) );
  }
  queue->removedCells[queue->transactions] = 0;
  queue->addedCells[queue->transactions] = 0;
  queue->removedFaces[queue->transactions] = 0;
  queue->addedFaces[queue->transactions] = 0;
  queue->removedEdges[queue->transactions] = 0;
  queue->addedEdges[queue->transactions] = 0;
  queue->transactions++;
  return queue;
}

Queue *queueRemoveCell( Queue *queue, int *nodes, int *nodeParts )
{
  int i;
  if (NULL==queue) return NULL;
  queue->removedCells[queue->transactions-1]++;
  if (queue->nRemovedCells>=queue->maxRemovedCells) {
    queue->maxRemovedCells += 100;
    queue->removedCellNodes = realloc( queue->removedCellNodes, 
				       8 * queue->maxRemovedCells* sizeof(int));
  }
  for (i=0;i<4;i++) {
    queue->removedCellNodes[i+8*queue->nRemovedCells]= nodes[i];
    queue->removedCellNodes[4+i+8*queue->nRemovedCells]= nodeParts[i];
  }
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
  for (i=0;i<4;i++) nodes[i] = queue->removedCellNodes[i+8*index];
  return queue;
}

Queue *queueRemovedCellNodeParts( Queue *queue, int index, int *nodeParts )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nRemovedCells ) return NULL;
  for (i=0;i<4;i++) nodeParts[i] = queue->removedCellNodes[4+i+8*index];
  return queue;
}

Queue *queueAddCell( Queue *queue, int *nodes, int cellId, int *nodeParts,
		     double *xyzs )
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
  for (i=0;i<4 ;i++) queue->addedCellNodes[i+9*queue->nAddedCells] = nodes[i];
  queue->addedCellNodes[4+9*queue->nAddedCells] = cellId;
  for (i=0;i<4 ;i++) 
    queue->addedCellNodes[5+i+9*queue->nAddedCells] = nodeParts[i];
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
  for (i=0;i<4;i++) nodes[i] = queue->addedCellNodes[i+9*index];
  return queue;
}

Queue *queueAddedCellId( Queue *queue, int index, int *cellId )
{
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedCells ) return NULL;
  *cellId = queue->addedCellNodes[4+9*index];
  return queue;
}

Queue *queueAddedCellNodeParts( Queue *queue, int index, int *nodeParts )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedCells ) return NULL;
  for (i=0;i<4;i++) nodeParts[i] = queue->addedCellNodes[5+i+9*index];
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

Queue *queueRemoveFace( Queue *queue, int *nodes, int *nodeParts )
{
  int i;
  if (NULL==queue) return NULL;
  queue->removedFaces[queue->transactions-1]++;
  if (queue->nRemovedFaces>=queue->maxRemovedFaces) {
    queue->maxRemovedFaces += 100;
    queue->removedFaceNodes = realloc( queue->removedFaceNodes, 
				       6 * queue->maxRemovedFaces* sizeof(int));
  }
  for (i=0;i<3;i++) {
    queue->removedFaceNodes[i+6*queue->nRemovedFaces]= nodes[i];
    queue->removedFaceNodes[3+i+6*queue->nRemovedFaces]= nodeParts[i];
  }
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
  for (i=0;i<3;i++) nodes[i] = queue->removedFaceNodes[i+6*index];
  return queue;
}

Queue *queueRemovedFaceNodeParts( Queue *queue, int index, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nRemovedFaces ) return NULL;
  for (i=0;i<3;i++) nodes[i] = queue->removedFaceNodes[3+i+6*index];
  return queue;
}

Queue *queueAddFace( Queue *queue, int *nodes, int faceId, int *nodeParts,
		     double *uvs )
{
  int i;
  if (NULL==queue) return NULL;
  queue->addedFaces[queue->transactions-1]++;
  if (queue->nAddedFaces>=queue->maxAddedFaces) {
    queue->maxAddedFaces += 100;
    queue->addedFaceNodes = realloc( queue->addedFaceNodes, 
				     7 * queue->maxAddedFaces * sizeof(int) );
    queue->addedFaceUVs   = realloc( queue->addedFaceUVs, 
				     6 * queue->maxAddedFaces* sizeof(double));
  }
  for (i=0;i<3;i++) {
    queue->addedFaceNodes[i+7*queue->nAddedFaces] = nodes[i];
    queue->addedFaceNodes[4+i+7*queue->nAddedFaces] = nodeParts[i];
  }
  queue->addedFaceNodes[3+7*queue->nAddedFaces] = faceId;
  for (i=0;i<6;i++) queue->addedFaceUVs[i+6*queue->nAddedFaces] = uvs[i];
  queue->nAddedFaces++;
  return queue;
}

Queue *queueAddFaceScalar( Queue *queue, 
			   int n0, int p0, double u0, double v0,
			   int n1, int p1, double u1, double v1,
			   int n2, int p2, double u2, double v2, int faceId)
{
  int nodes[3];
  int nodeParts[3];
  double uv[6];
  nodes[0] = n0;
  nodes[1] = n1;
  nodes[2] = n2;
  nodeParts[0] = p0;
  nodeParts[1] = p1;
  nodeParts[2] = p2;
  uv[0] = u0;
  uv[1] = v0;
  uv[2] = u1;
  uv[3] = v1;
  uv[4] = u2;
  uv[5] = v2;
  return queueAddFace( queue, nodes, faceId, nodeParts, uv );
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
  for (i=0;i<3;i++) nodes[i] = queue->addedFaceNodes[i+7*index];
  return queue;
}

Queue *queueAddedFaceNodeParts( Queue *queue, int index, int *nodeParts )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedFaces ) return NULL;
  for (i=0;i<3;i++) nodeParts[i] = queue->addedFaceNodes[4+i+7*index];
  return queue;
}

Queue *queueAddedFaceId( Queue *queue, int index, int *faceId )
{
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedFaces ) return NULL;
  *faceId = queue->addedFaceNodes[3+7*index];
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

Queue *queueRemoveEdge( Queue *queue, int *nodes, int *nodeParts )
{
  int i;
  if (NULL==queue) return NULL;
  queue->removedEdges[queue->transactions-1]++;
  if (queue->nRemovedEdges>=queue->maxRemovedEdges) {
    queue->maxRemovedEdges += 100;
    queue->removedEdgeNodes = realloc( queue->removedEdgeNodes, 
				       4 * queue->maxRemovedEdges* sizeof(int));
  }
  for (i=0;i<2;i++) {
    queue->removedEdgeNodes[i+4*queue->nRemovedEdges]= nodes[i];
    queue->removedEdgeNodes[2+i+4*queue->nRemovedEdges]= nodeParts[i];
  }
  queue->nRemovedEdges++;
  return queue;
}

int queueRemovedEdges( Queue *queue, int transaction )
{
  if (NULL==queue) return EMPTY;
  if ( transaction<0 || transaction>=queueTransactions(queue) ) return EMPTY;
  return queue->removedEdges[transaction];
}

Queue *queueRemovedEdgeNodes( Queue *queue, int index, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nRemovedEdges ) return NULL;
  for (i=0;i<2;i++) nodes[i] = queue->removedEdgeNodes[i+4*index];
  return queue;
}

Queue *queueRemovedEdgeNodeParts( Queue *queue, int index, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nRemovedEdges ) return NULL;
  for (i=0;i<2;i++) nodes[i] = queue->removedEdgeNodes[2+i+4*index];
  return queue;
}

Queue *queueAddEdge( Queue *queue, int *nodes, int edgeId, int *nodeParts,
		     double *uvs )
{
  int i;
  if (NULL==queue) return NULL;
  queue->addedEdges[queue->transactions-1]++;
  if (queue->nAddedEdges>=queue->maxAddedEdges) {
    queue->maxAddedEdges += 100;
    queue->addedEdgeNodes = realloc( queue->addedEdgeNodes, 
				     5 * queue->maxAddedEdges * sizeof(int) );
    queue->addedEdgeTs   = realloc( queue->addedEdgeTs, 
				     2 * queue->maxAddedEdges* sizeof(double));
  }
  for (i=0;i<2;i++) {
    queue->addedEdgeNodes[i+5*queue->nAddedEdges] = nodes[i];
    queue->addedEdgeNodes[3+i+5*queue->nAddedEdges] = nodeParts[i];
  }
  queue->addedEdgeNodes[2+5*queue->nAddedEdges] = edgeId;
  for (i=0;i<2;i++) queue->addedEdgeTs[i+2*queue->nAddedEdges] = uvs[i];
  queue->nAddedEdges++;
  return queue;
}

int queueAddedEdges( Queue *queue, int transaction )
{
  if (NULL==queue) return EMPTY;
  if ( transaction<0 || transaction>=queueTransactions(queue) ) return EMPTY;
  return queue->addedEdges[transaction];
}

Queue *queueAddedEdgeNodes( Queue *queue, int index, int *nodes )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedEdges ) return NULL;
  for (i=0;i<2;i++) nodes[i] = queue->addedEdgeNodes[i+5*index];
  return queue;
}

Queue *queueAddedEdgeNodeParts( Queue *queue, int index, int *nodeParts )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedEdges ) return NULL;
  for (i=0;i<2;i++) nodeParts[i] = queue->addedEdgeNodes[3+i+5*index];
  return queue;
}

Queue *queueAddedEdgeId( Queue *queue, int index, int *edgeId )
{
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedEdges ) return NULL;
  *edgeId = queue->addedEdgeNodes[2+5*index];
  return queue;
}

Queue *queueAddedEdgeTs( Queue *queue, int index, double *xyzs )
{
  int i;
  if (NULL==queue) return NULL;
  if ( index<0 || index>queue->nAddedEdges ) return NULL;
  for (i=0;i<2;i++) xyzs[i] = queue->addedEdgeTs[i+2*index];
  return queue;
}

Queue *queueDumpSize( Queue *queue, int *nInt, int *nDouble )
{
  *nInt = EMPTY;
  *nDouble = EMPTY;
  if (NULL==queue) return NULL;
  *nInt 
    = 8
    + 6 * queue->transactions
    + 8 * queue->nRemovedCells
    + 9 * queue->nAddedCells
    + 6 * queue->nRemovedFaces
    + 7 * queue->nAddedFaces
    + 4 * queue->nRemovedEdges
    + 5 * queue->nAddedEdges;
  *nDouble = 
    4 * queue->nodeSize * queue->nAddedCells + 
    6 * queue->nAddedFaces +
    2 * queue->nAddedEdges;
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
  ints[6] = queue->nRemovedEdges;
  ints[7] = queue->nAddedEdges;
  i = 8;
  d = 0;

  for(transaction=0;transaction<queue->transactions;transaction++){
    ints[i] = queue->removedCells[transaction]; i++;
  }
  for(removed=0;removed<queue->nRemovedCells;removed++){
    size = 8;
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
    size = 6;
    for (node=0;node<size;node++) { 
      ints[i] = queue->removedFaceNodes[node+size*removed]; i++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    ints[i] = queue->addedFaces[transaction]; i++;
  }
  for(added=0;added<queue->nAddedFaces;added++){
    size = 7;
    for (node=0;node<size;node++) { 
      ints[i] = queue->addedFaceNodes[node+size*added]; i++;
    }
    size = 6;
    for (node=0;node<size;node++) { 
      doubles[d] = queue->addedFaceUVs[node+size*added]; d++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    ints[i] = queue->removedEdges[transaction]; i++;
  }
  for(removed=0;removed<queue->nRemovedEdges;removed++){
    size = 4;
    for (node=0;node<size;node++) { 
      ints[i] = queue->removedEdgeNodes[node+size*removed]; i++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    ints[i] = queue->addedEdges[transaction]; i++;
  }
  for(added=0;added<queue->nAddedEdges;added++){
    size = 5;
    for (node=0;node<size;node++) { 
      ints[i] = queue->addedEdgeNodes[node+size*added]; i++;
    }
    size = 2;
    for (node=0;node<size;node++) { 
      doubles[d] = queue->addedEdgeTs[node+size*added]; d++;
    }
  }

  return queue;
}

Queue *queueLoad( Queue *queue, int *ints, double *doubles )
{
  int i, d, node, size, transaction, removed, added;

  if (NULL==queue) return NULL;

  if ( queue->nodeSize != ints[0] ) {
    /* printf("ERROR: %s: %d: queueLoad: incompatable nodeSize: %d %d\n",
              __FILE__, __LINE__, queue->nodeSize, ints[0] ); */
    return NULL;
  }			    
  queue->transactions  = ints[1];
  queue->nRemovedCells = ints[2];
  queue->nAddedCells   = ints[3];
  queue->nRemovedFaces = ints[4];
  queue->nAddedFaces   = ints[5];
  queue->nRemovedEdges = ints[6];
  queue->nAddedEdges   = ints[7];
  i = 8;
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
    queue->removedEdges = realloc( queue->removedEdges, 
				   queue->maxTransactions * sizeof(int) );
    queue->addedEdges   = realloc( queue->addedEdges, 
				   queue->maxTransactions * sizeof(int) );
  }
  if (queue->nRemovedCells > queue->maxRemovedCells) {
    queue->maxRemovedCells = queue->nRemovedCells;
    queue->removedCellNodes = realloc( queue->removedCellNodes, 
				       8 * queue->maxRemovedCells* sizeof(int));
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
				       6 * queue->maxRemovedFaces* sizeof(int));
  }
  if (queue->nAddedFaces > queue->maxAddedFaces) {
    queue->maxAddedFaces = queue->nAddedFaces;
    queue->addedFaceNodes = realloc( queue->addedFaceNodes, 
				     7 * queue->maxAddedFaces * sizeof(int) );
    queue->addedFaceUVs   = realloc( queue->addedFaceUVs, 
				     6 * queue->maxAddedFaces* sizeof(double));
  }

  if (queue->nRemovedEdges > queue->maxRemovedEdges) {
    queue->maxRemovedEdges = queue->nRemovedEdges;
    queue->removedEdgeNodes = realloc( queue->removedEdgeNodes, 
				       4 * queue->maxRemovedEdges* sizeof(int));
  }
  if (queue->nAddedEdges > queue->maxAddedEdges) {
    queue->maxAddedEdges = queue->nAddedEdges;
    queue->addedEdgeNodes = realloc( queue->addedEdgeNodes, 
				     5 * queue->maxAddedEdges * sizeof(int) );
    queue->addedEdgeTs    = realloc( queue->addedEdgeTs, 
				     2 * queue->maxAddedEdges* sizeof(double));
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    queue->removedCells[transaction] = ints[i]; i++;
  }
  for(removed=0;removed<queue->nRemovedCells;removed++){
    size = 8;
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
    size = 6;
    for (node=0;node<size;node++) { 
      queue->removedFaceNodes[node+size*removed] = ints[i]; i++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    queue->addedFaces[transaction] = ints[i]; i++;
  }
  for(added=0;added<queue->nAddedFaces;added++){
    size = 7;
    for (node=0;node<size;node++) { 
      queue->addedFaceNodes[node+size*added] = ints[i]; i++;
    }
    size = 6;
    for (node=0;node<size;node++) { 
      queue->addedFaceUVs[node+size*added] = doubles[d]; d++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    queue->removedEdges[transaction] = ints[i]; i++;
  }
  for(removed=0;removed<queue->nRemovedEdges;removed++){
    size = 4;
    for (node=0;node<size;node++) { 
      queue->removedEdgeNodes[node+size*removed] = ints[i]; i++;
    }
  }

  for(transaction=0;transaction<queue->transactions;transaction++){
    queue->addedEdges[transaction] = ints[i]; i++;
  }
  for(added=0;added<queue->nAddedEdges;added++){
    size = 5;
    for (node=0;node<size;node++) { 
      queue->addedEdgeNodes[node+size*added] = ints[i]; i++;
    }
    size = 2;
    for (node=0;node<size;node++) { 
      queue->addedEdgeTs[node+size*added] = doubles[d]; d++;
    }
  }

  return queue;
}

Queue *queueGlobalShiftNode(Queue *queue,
			    int old_nnode_global,
			    int node_offset )
{
  int index, i;

  for ( index = 0 ; index < queue->nAddedCells ; index++ )
    for ( i = 0; i < 4 ; i++ ) 
      if ( queue->addedCellNodes[i+9*index] >= old_nnode_global )
	queue->addedCellNodes[i+9*index] += node_offset;

  for ( index = 0 ; index < queue->nRemovedCells ; index++ )
    for ( i = 0; i < 4 ; i++ ) 
      if ( queue->removedCellNodes[i+8*index] >= old_nnode_global )
	queue->removedCellNodes[i+8*index] += node_offset;

  for ( index = 0 ; index < queue->nAddedFaces ; index++ )
    for ( i = 0; i < 3 ; i++ ) 
      if ( queue->addedFaceNodes[i+7*index] >= old_nnode_global )
	queue->addedFaceNodes[i+7*index] += node_offset;

  for ( index = 0 ; index < queue->nRemovedFaces ; index++ )
    for ( i = 0; i < 3 ; i++ ) 
      if ( queue->removedFaceNodes[i+6*index] >= old_nnode_global )
	queue->removedFaceNodes[i+6*index] += node_offset;

  for ( index = 0 ; index < queue->nAddedEdges ; index++ )
    for ( i = 0; i < 2 ; i++ ) 
      if ( queue->addedEdgeNodes[i+5*index] >= old_nnode_global )
	queue->addedEdgeNodes[i+5*index] += node_offset;

  for ( index = 0 ; index < queue->nRemovedEdges ; index++ )
    for ( i = 0; i < 2 ; i++ ) 
      if ( queue->removedEdgeNodes[i+4*index] >= old_nnode_global )
	queue->removedEdgeNodes[i+4*index] += node_offset;

  return queue;
}

Queue *queueGlobalShiftCell(Queue *queue,
			    int old_ncell_global,
			    int cell_offset )
{
  int index;

  for ( index = 0 ; index < queue->nAddedCells ; index++ )
    if ( queue->addedCellNodes[4+9*index] >= old_ncell_global )
      queue->addedCellNodes[4+9*index] += cell_offset;

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
