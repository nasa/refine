
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
Layer *formAdvancingFront( Grid *grid, char *project );
Grid *layerGrid(Layer *);
void layerFree(Layer *);
void layerSortGlobalNodes(void *layer, int *o2n);
int layerNFront(Layer *);
int layerNBlend(Layer *);
int layerNNormal(Layer *);
int layerMaxNode(Layer *);
Layer *layerMakeFront(Layer *, int nbc, int *bc);
bool layerParentFace(Layer *, int faceId);
Layer *layerFront(Layer *, int front, int *nodes);
Layer *layerFrontDirection(Layer *, int front, double *direction);
Layer *layerMakeNormal(Layer *);
Layer *layerInitializeNormal(Layer *, int normal );
Layer *layerCopyNormal(Layer *, int originalNormal, int newNormal );
Layer *layerFrontNormals(Layer *, int front, int *normals);
int layerNormalRoot(Layer *, int normal );
int layerNormalDeg(Layer *, int normal );
Layer *layerNormalFronts(Layer *, int normal, int maxfront, int *fronts);
Layer *layerNormalDirection(Layer *, int normal, double *direction);
Layer *layerSetNormalHeight(Layer *, int normal, double height);
Layer *layerGetNormalHeight(Layer *, int normal, double *height);
Layer *layerScaleNormalHeight(Layer *, double scale);
Layer *layerLaminarInitialHeight(Layer *, double Re, double xStart );
Layer *layerVisibleNormals(Layer *);
Layer *layerConstrainNormal(Layer *, int edgeface );
bool layerConstrainingGeometry(Layer *, int edgeface );
int layerConstrained(Layer *, int normal );
Layer *layerConstrainFrontSide(Layer *, int normal0, int normal1, int bc );
int layerConstrainedSide(Layer *, int front, int side );
int layerNConstrainedSides(Layer *, int faceId );
Layer *layerFindParentEdges(Layer *);
Layer *layerSetParentEdge(Layer *, int normal0, int normal1, int edgeId );
int layerParentEdge(Layer *, int front, int side );
int layerNParentEdgeSegments(Layer *, int edgeId );

Layer *layerTerminateNormal(Layer *, int normal );
bool layerNormalTerminated(Layer *, int normal );
int layerNActiveNormal(Layer *);
Layer *layerAdvance(Layer * );
Layer *layerAdvanceConstantHeight(Layer *, double height );
Layer *layerWiggle(Layer *, double height );
Layer *layerBlend(Layer *);

Layer *layerSmoothLayerNeighbors(Layer * );
Layer *layerTerminateNormalWithSpacing(Layer *, double spacing);
Layer *layerTerminateNormalWithX(Layer *, int direction, double x);

Layer *layerInsertPhantomFront(Layer *, double dz);
Layer *layerVerifyPhantomEdges(Layer *);
Layer *layerVerifyPhantomFaces(Layer *);

END_C_DECLORATION

#endif /* LAYER_H */
