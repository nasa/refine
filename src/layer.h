
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
int layerMaxTriangle(Layer *);
int layerNTriangle(Layer *);
int layerMaxNormal(Layer *);
int layerNNormal(Layer *);
int layerMaxNode(Layer *);
Layer *layerPopulateAdvancingFront(Layer *, int nbc, int *bc);
Layer *layerAddParentGeomFace(Layer *, int faceId);
bool layerParentGeomFace(Layer *, int faceId);
Layer *layerAddTriangle(Layer *, int n0, int n1, int n2);
Layer *layerTriangle(Layer *, int triangle, int *nodes);
Layer *layerTriangleDirection(Layer *, int triangle, double *direction);
int layerAddNormal(Layer *, int globalNodeId );
int layerUniqueNormalId(Layer *, int globalNodeId );
Layer *layerInitializeNormal(Layer *, int normal );
Layer *layerTriangleNormals(Layer *, int triangle, int *normals);
int layerNormalRoot(Layer *, int normal );
int layerNormalDeg(Layer *, int normal );
Layer *layerNormalTriangles(Layer *, int normal, int maxtriangle, int *triangles);
Layer *layerNormalDirection(Layer *, int normal, double *direction);
Layer *layerSetHeightOfAllNormals(Layer *, double height);
Layer *layerSetNormalHeightOfFace(Layer *, int faceId, double height);
Layer *layerSetNormalHeight(Layer *, int normal, double height);
Layer *layerGetNormalHeight(Layer *, int normal, double *height);
Layer *layerScaleNormalHeight(Layer *, double scale);
Layer *layerLaminarInitialHeight(Layer *, double Re, double xStart );
Layer *layerLaminarInitialHeightNegZ(Layer *);
Layer *layerVisibleNormals(Layer *, double dotLimit, double radianLimit );
Layer *layerSmoothNormalDirection(Layer *);
Layer *layerProjectNormalsToConstraints(Layer *);

Layer *layerConstrainNormal(Layer *, int edgeface );
bool layerConstrainingGeometry(Layer *, int edgeface );
int layerConstrained(Layer *, int normal );
Layer *layerConstrainTriangleSide(Layer *, int normal0, int normal1, int bc );
int layerConstrainedSide(Layer *, int triangle, int side );
int layerNConstrainedSides(Layer *, int faceId );
Layer *layerFindParentGeomEdges(Layer *);
Layer *layerSetParentGeomEdge(Layer *, int normal0, int normal1, int edgeId );
int layerParentGeomEdge(Layer *, int triangle, int side );
int layerNParentGeomEdgeSegments(Layer *, int edgeId );

Layer *layerTerminateNormal(Layer *, int normal );
bool layerNormalTerminated(Layer *, int normal );
Layer *layerTerminateFaceNormals(Layer *, int faceId );
int layerNActiveNormal(Layer *);

bool layerCellInLayer(Layer *, int cell);
bool layerFaceInLayer(Layer *, int face);
bool layerEdgeInLayer(Layer *, int edge);

Layer *layerReconnectCellUnlessInLayer(Layer *, int oldNode, int newNode );
Layer *layerReconnectEdgeUnlessInLayer(Layer *, int edgeId, 
				       int oldNode, int newNode );
Layer *layerAdvance(Layer * );
Layer *layerAdvanceConstantHeight(Layer *, double height );
Layer *layerWiggle(Layer *, double height );

Layer *layerSmoothLayerNeighbors(Layer * );
Layer *layerTerminateNormalWithSpacing(Layer *, double spacing);
Layer *layerTerminateNormalWithX(Layer *, int direction, double x);

Layer *layerInsertPhantomTriangle(Layer *, double dz);
Layer *layerVerifyPhantomEdges(Layer *);
Layer *layerVerifyPhantomFaces(Layer *);

bool layerTetrahedraOnly(Layer *);
Layer *layerToggleMixedElementMode(Layer *);

END_C_DECLORATION

#endif /* LAYER_H */
