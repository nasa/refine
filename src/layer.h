
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef LAYER_H
#define LAYER_H

#include "master_header.h"
#include "grid.h"

BEGIN_C_DECLORATION

typedef struct Layer Layer;

Layer *layerCreate(Grid *);
void layerFree(Layer *);
void layerSortGlobalNodes(void *layer, int *o2n);
int layerNFront(Layer *);
int layerNNormal(Layer *);
int layerMaxNode(Layer *);
Layer *layerMakeFront(Layer *, int nbc, int *bc);
Layer *layerFront(Layer *, int front, int *nodes);
Layer *layerFrontDirection(Layer *, int front, double *direction);
Layer *layerMakeNormal(Layer *);
Layer *layerFrontNormals(Layer *, int front, int *normals);
int layerNormalRoot(Layer *, int normal );
int layerNormalDeg(Layer *, int normal );
Layer *layerNormalFronts(Layer *, int normal, int maxfront, int *fronts);
Layer *layerNormalDirection(Layer *, int normal, double *direction);
Layer *layerVisibleNormals(Layer *);
Layer *layerConstrainNormal(Layer *, int bc );
Layer *layerConstrainFrontSide(Layer *, int normal0, int normal1, int bc );
int layerConstrained(Layer *, int normal );
int layerConstrainedSide(Layer *, int front, int side );
Layer *layerTerminateNormal(Layer *, int normal );
bool layerNormalTerminated(Layer *, int normal );
int layerNActiveNormal(Layer *);
Layer *layerAdvance(Layer *, double height );
Layer *layerWiggle(Layer *, double height );
Layer *layerSmoothLayerNeighbors(Layer * );
Layer *layerTerminateNormalWithSpacing(Layer *, double spacing);

Layer *layerInsertPhantomFront(Layer *, double dz);
Layer *layerVerifyPhantomEdges(Layer *);
Layer *layerVerifyPhantomFaces(Layer *);

END_C_DECLORATION

#endif /* LAYER_H */
