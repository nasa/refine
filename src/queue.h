
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

#include "refine_defs.h"
#include <stdio.h>

BEGIN_C_DECLORATION

typedef struct Queue Queue;

struct Queue {
  int nodeSize;
  int transactions;
  int maxTransactions;
  int *removedCells;
  int maxRemovedCells;
  int nRemovedCells;
  int *removedCellNodes;
  int *addedCells;
  int maxAddedCells;
  int nAddedCells;
  int *addedCellNodes;
  double *addedCellXYZs;
  int *removedFaces;
  int maxRemovedFaces;
  int nRemovedFaces;
  int *removedFaceNodes;
  int *addedFaces;
  int maxAddedFaces;
  int nAddedFaces;
  int *addedFaceNodes;
  double *addedFaceUVs;
};

Queue *queueCreate( int nodeSize );
void queueFree( Queue * );
Queue *queueReset( Queue * );
int queueNodeSize( Queue * );
int queueTransactions( Queue * );
Queue *queueNewTransaction( Queue * );

Queue *queueRemoveCell( Queue *, int *nodes );
int queueRemovedCells( Queue *, int transaction );
Queue *queueRemovedCellNodes( Queue *, int index, int *nodes );
Queue *queueAddCell( Queue *, int *nodes, int cellId, int *nodeParts,
		     double *xyzs );
int queueAddedCells( Queue *, int transaction );
Queue *queueAddedCellNodes( Queue *, int index, int *nodes );
Queue *queueAddedCellId( Queue *, int index, int *cellId );
Queue *queueAddedCellNodeParts( Queue *, int index, int *nodeParts );
Queue *queueAddedCellXYZs( Queue *, int index, double *xyzs );
int queueTotalRemovedCells( Queue * );

Queue *queueRemoveFace( Queue *, int *nodes );
int queueRemovedFaces( Queue *, int transaction );
Queue *queueRemovedFaceNodes( Queue *, int index, int *nodes );
Queue *queueAddFace( Queue *, int *nodes, double *uvs );
Queue *queueAddFaceScalar( Queue *, 
			   int n0, double u0, double v0,
			   int n1, double u1, double v1,
			   int n2, double u2, double v2, int faceId);
int queueAddedFaces( Queue *, int transaction );
Queue *queueAddedFaceNodes( Queue *, int index, int *nodes );
Queue *queueAddedFaceUVs( Queue *, int index, double *uvs );

Queue *queueDumpSize( Queue *, int *nInt, int *nDouble );
Queue *queueDump( Queue *, int *ints, double *doubles );
Queue *queueLoad( Queue *, int *ints, double *doubles );

Queue *queueContents(Queue *queue, FILE *file);

END_C_DECLORATION

#endif /* QUEUE_H */
