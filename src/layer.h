
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
Grid *layerGrid(Layer *);
void layerFree(Layer *);
Layer *formAdvancingFront( Grid *grid, char *project );
void layerSortGlobalNodes(void *layer, int *o2n);
int layerMaxTriangle(Layer *);
int layerNTriangle(Layer *);
int layerNBlend(Layer *);
int layerMaxNormal(Layer *);
int layerNNormal(Layer *);
int layerMaxNode(Layer *);
Layer *layerPopulateAdvancingFront(Layer *, int nbc, int *bc);
Layer *layerBuildNormalTriangleAdjacency(Layer *);
Layer *layerInitializeTriangleNormalDirection(Layer *);
Layer *layerAddParentGeomFace(Layer *, int faceId);
bool layerParentGeomFace(Layer *, int faceId);
Layer *layerAddTriangle(Layer *, int n0, int n1, int n2);
int layerForceTriangle(Layer *, int normal0, int normal1, int nnormal2);
Layer *layerTriangle(Layer *, int triangle, int *nodes);
Layer *layerTriangleDirection(Layer *, int triangle, double *direction);
Layer *layerTriangleArea(Layer *, int triangle, double *area);
Layer *layerTriangleCenter(Layer *, int triangle, double *center);
Layer *layerTriangleFourthNode(Layer *, int triangle, double *xyz);
Layer *layerTriangleMaxEdgeLength(Layer *, int triangle, double *length);
Layer *layerNormalMaxEdgeLength(Layer *, int normal, double *length);
int layerAddNormal(Layer *, int globalNodeId );
int layerUniqueNormalId(Layer *, int globalNodeId );
int layerDuplicateNormal(Layer *, int normal );
Layer *layerInitializeNormal(Layer *, int normal );
Layer *layerTriangleNormals(Layer *, int triangle, int *normals);
int layerNormalRoot(Layer *, int normal );
int layerNormalDeg(Layer *, int normal );
Layer *layerNormalTriangles(Layer *, int normal, int maxtriangle, int *triangles);
int layerPreviousTriangle(Layer *, int normal, int triangle );
int layerNextTriangle(Layer *, int normal, int triangle );
Layer *layerCommonEdge(Layer *, int triangle0, int triangle1, int *nodes);
double layerEdgeAngle(Layer *, int triangle0, int triangle1 );
Layer *layerNormalDirection(Layer *, int normal, double *direction);
Layer *layerAssignPolynomialNormalHeight(Layer *, double constant, double slope, 
                                         double exponent, double *origin,
					 double *direction);
Layer *layerSetHeightOfAllNormals(Layer *, double height);
Layer *layerLaminarInitialHeight(Layer *, double Re, double xStart );
Layer *layerLaminarInitialHeightNegZ(Layer *);
Layer *layerSetNormalHeightOfFace(Layer *, int faceId, double height);
Layer *layerSetNormalHeight(Layer *, int normal, double height);
Layer *layerGetNormalHeight(Layer *, int normal, double *height);
Layer *layerScaleNormalHeight(Layer *, double scale);
Layer *layerScaleNormalHeightWithPolynomial(Layer *, 
					    double constant, double slope,
				            double exponent, double scale,
					    double *origin, double *direction);
Layer *layerSetNormalMaxLength(Layer *, int normal, double maxLength);
double layerNormalMaxLength(Layer *, int normal);
Layer *layerSetPolynomialMaxLength(Layer *, double constant, double slope, 
				            double exponent, double *origin,
				            double *direction);
Layer *layerSetNormalInitialHeight(Layer *layer, int normal, 
				   double initialHeight);
Layer *layerSaveInitialNormalHeight(Layer *);
double layerNormalInitialHeight(Layer *, int normal);
Layer *layerSetNormalRate(Layer *, int normal, double rate);
double layerNormalRate(Layer *, int normal);
Layer *layerSetAllNormalRate(Layer *, double rate);
Layer *layerSetNormalHeightWithRate(Layer *);
Layer *layerSetNormalHeightWithMaxRate(Layer *, double maxRate);
Layer *layerVisibleNormals(Layer *, double dotLimit, double radianLimit );
Layer *layerSmoothNormalDirection(Layer *);
Layer *layerProjectNormalsToConstraints(Layer *);

Layer *layerAdjustNormalHeightToSmoothFront(Layer *, double maxHeight);

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
bool layerAnyActiveNormals(Layer *);

bool layerCellInLayer(Layer *, int cell);
bool layerFaceInLayer(Layer *, int face);
bool layerEdgeInLayer(Layer *, int edge);
int layerNEdgeInLayer(Layer *, int edgeId);
int layerEdgeEndPoint(Layer *, int edgeId, int startNode);

Layer *layerReconnectCellUnlessInLayer(Layer *, int oldNode, int newNode );
Layer *layerReconnectEdgeUnlessInLayer(Layer *, int edgeId, 
				       int oldNode, int newNode );
Layer *layerReconnectFaceUnlessInLayer(Layer *, int faceId, 
				       int oldNode, int newNode );

Layer *layerAdvance(Layer * );
Layer *layerAdvanceConstantHeight(Layer *, double height );
Layer *layerWiggle(Layer *, double height );

Layer *layerSmoothLayerNeighbors(Layer * );
Layer *layerTerminateNormalWithSpacing(Layer *, double spacing);
Layer *layerTerminateNormalWithX(Layer *, int direction, double x);
int layerTerminateNormalWithLength(Layer *, double ratio);

Layer *layerInsertPhantomTriangle(Layer *, double dz);
Layer *layerVerifyPhantomEdges(Layer *);
Layer *layerVerifyPhantomFaces(Layer *);

Layer *layerThaw(Layer*);

bool layerTetrahedraOnly(Layer *);
Layer *layerToggleMixedElementMode(Layer *);

Adj *layerBuildNormalBlendAdjacency(Layer *layer);
Layer *layerSplitBlend(Layer *);
Layer *layerBlend(Layer *);
Layer *layerAddBlend(Layer *, int normal0, int normal1, int otherNode );
Layer *layerDuplicateAllBlend(Layer *);
Layer *layerBlendNormals(Layer *, int blend, int *normals );
Layer *layerExtrudeBlend(Layer *, double dx, double dy, double dz );

Layer *layerPopulateNormalNearTree(Layer *);
Layer *layerPopulateTriangleNearTree(Layer *);
Layer *layerTerminateCollidingNormals(Layer *);

Layer *layerWriteTecplotFrontGeometry(Layer *);
Layer *layerWriteTecplotFrontWithData(Layer *, int);

END_C_DECLORATION

#endif /* LAYER_H */
