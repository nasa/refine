
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
Queue *queueAddedCellXYZs( Queue *, int index, double *xyzs );

Queue *queueRemoveFace( Queue *, int *nodes );
int queueRemovedFaces( Queue *, int transaction );
Queue *queueRemovedFaceNodes( Queue *, int index, int *nodes );
Queue *queueAddFace( Queue *, int *nodes, double *uvs );
int queueAddedFaces( Queue *, int transaction );
Queue *queueAddedFaceNodes( Queue *, int index, int *nodes );
Queue *queueAddedFaceUVs( Queue *, int index, double *uvs );

END_C_DECLORATION

#endif /* QUEUE_H */
