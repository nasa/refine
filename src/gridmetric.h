
/* Computes metrics from faces and tets 
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDMETRIC_H
#define GRIDMETRIC_H

#include "refine_defs.h"
#include "grid.h"
#include "gridmath.h"

BEGIN_C_DECLORATION

Grid *gridSetMapMatrixToAverageOfNodes(Grid *g, int avgNode, int n0, int n1 );
void gridMapXYZWithJ( double *j, double *x, double *y, double *z );

double gridEdgeLength(Grid *g, int n0, int n1 );
double gridEdgeRatio(Grid *g, int n0, int n1 );
double gridAverageEdgeLength(Grid *g, int node );
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
Grid *gridCopySpacing(Grid *g, int originalNode, int newNode);

Grid *gridConvertMetricToJacobian(Grid *g, double *m, double *j);

double gridVolume(Grid *g, int *nodes );
Grid *gridCellVolumeDerivative(Grid *g, int *nodes, 
			       double *volume, double *dVoldx );
Grid *gridNodeVolumeDerivative(Grid *g, int node, 
			       double *volume, double *dVoldx );

double gridAR(Grid *g, int *nodes );
double gridCellAspectRatio( double *n0, double *n1, double *n2, double *n3 );
Grid *gridNodeAR(Grid *g, int node, double *ar );
Grid *gridNodeVolume(Grid *g, int node, double *volume );
Grid *gridGemAR(Grid *g, double *ar);
Grid *gridCellARDerivative(Grid *g, int *nodes, double *ar, double *dARdx );

void gridCellAspectRatioDerivative( double *xyz1, double *xyz2, 
				    double *xyz3, double *xyz4,
				    double *ar, double *dARdx);

Grid *gridNodeARDerivative(Grid *g, int node, double *ar, double *dARdx );
Grid *gridStoreARDerivative(Grid *g, int node );
double gridMinVolume(Grid *g);
GridBool gridNegCellAroundNode(Grid *g, int node );
GridBool gridNegCellAroundNodeExceptGem(Grid *g, int node );
double gridMinARAroundNodeExceptGem(Grid *g, int node );
double gridMinAR(Grid *g);
double gridMinThawedAR(Grid *g);

Grid *gridFreezeSmallARCells(Grid *g, double minAR );

GridBool gridRightHandedFace(Grid *g, int face );
GridBool gridRightHandedBoundary(Grid *g );

double gridFaceArea(Grid *g, int n0, int n1, int n2);
double gridFaceAR(Grid *g, int n0, int n1, int n2);
double gridFaceMR(Grid *g, int n0, int n1, int n2);
double gridMinFaceMR(Grid *g);
double gridMinThawedFaceMR(Grid *g);
Grid *gridFaceMRDerivative(Grid *g, int* nodes, double *mr, double *dMRdx );

void FaceMRDerivative(double x1, double y1, double z1,
		      double x2, double y2, double z2,
		      double x3, double y3, double z3,
		      double *mr, double *dMRdx  );
Grid *gridNodeFaceMR(Grid *g, int node, double *mr );
Grid *gridNodeFaceMRDerivative(Grid *g, int node, double *mr, double *dMRdx );

double gridCellMeanRatio( double *n0, double *n1, double *n2, double *n3 );

void gridCellMeanRatioDerivative( double *xyz0, double *xyz1, 
				  double *xyz2, double *xyz3,
				  double *mr, double *dMRdx);

END_C_DECLORATION

#endif /* GRIDMETRIC_H */
