
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
double gridSpacing(Grid *g, int node );
Grid *gridResetSpacing(Grid *g );
Grid *gridScaleSpacing(Grid *g, int node, double scale );
Grid *gridScaleSpacingSphere(Grid *g, double x, double y, double z, double r,
			     double scale );

double gridVolume(Grid *g, int *nodes );
double gridAR(Grid *g, int *nodes );
Grid *gridNodeAR(Grid *g, int node, double *ar );
Grid *gridCellARDerivative(Grid *g, int *nodes, double *ar, double *dARdx );
Grid *gridNodeARDerivative(Grid *g, int node, double *ar, double *dARdx );
double gridMinVolume(Grid *g);
bool gridNegCellAroundNode(Grid *g, int node );
double gridMinAR(Grid *g);

bool gridRightHandedFace(Grid *g, int face );
bool gridRightHandedBoundary(Grid *g );

END_C_DECLORATION

#endif /* GRIDMETRIC_H */
