
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDMETRIC_H
#define GRIDMETRIC_H

#include "master_header.h"
#include "grid.h"

BEGIN_C_DECLORATION

double gridEdgeLength(Grid *g, int n0, int n1 );
double gridAverageEdgeLength(Grid *g, int node );
int gridLongestEdge(Grid *g, int node );
Grid *gridLargestRatioEdge(Grid *g, int node, int *edgeNode, double *ratio );
Grid *gridSmallestRatioEdge(Grid *g, int node, int *edgeNode, double *ratio );
double gridSpacing(Grid *g, int node );
Grid *gridResetSpacing(Grid *g );
Grid *gridScaleSpacing(Grid *g, int node, double scale );
Grid *gridScaleSpacingSphere(Grid *g, double x, double y, double z, double r,
			     double scale );
Grid *gridScaleSpacingSphereDirection(Grid *g, 
				 double x, double y, double z, double r,
				 double scalex, double scaley, double scalez );

double gridVolume(Grid *g, int *nodes );
double gridAR(Grid *g, int *nodes );
Grid *gridNodeAR(Grid *g, int node, double *ar );
Grid *gridCellARDerivative(Grid *g, int *nodes, double *ar, double *dARdx );
Grid *gridNodeARDerivative(Grid *g, int node, double *ar, double *dARdx );
double gridMinVolume(Grid *g);
bool gridNegCellAroundNode(Grid *g, int node );
bool gridNegCellAroundNodeExceptGem(Grid *g, int node );
double gridMinAR(Grid *g);

bool gridRightHandedFace(Grid *g, int face );
bool gridRightHandedBoundary(Grid *g );

double gridFaceArea(Grid *g, int n0, int n1, int n2);
double gridFaceAR(Grid *g, int n0, int n1, int n2);
double gridFaceMR(Grid *g, int n0, int n1, int n2);
Grid *gridFaceMRDerivative(Grid *g, int* nodes, double *mr, double *dMRdx );

void FaceMRDerivative(double x1, double y1, double z1,
		      double x2, double y2, double z2,
		      double x3, double y3, double z3,
		      double *mr, double *dMRdx  );

Grid *gridMapMatrix(Grid *g, int node, double *m);
void gridMapXYZWithM( double *m, double *x, double *y, double *z );
void gridMapXYZWithNode( Grid *g, int node, double *x, double *y, double *z );

double gridMappedEdgeLength(Grid *g, int n0, int n1 );

END_C_DECLORATION

#endif /* GRIDMETRIC_H */
