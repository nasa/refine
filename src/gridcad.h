
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDCAD_H
#define GRIDCAD_H

#include "master_header.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridProjectNodeToEdge(Grid *g, int node, int edgeId );
Grid *gridProjectNodeToFace(Grid *g, int node, int faceId );
Grid *gridSafeProjectNode(Grid *g, int node );
Grid *gridSafeProjectNodeToFace(Grid *g, int node, int faceId );
Grid *gridSafeProjectNodeToEdge(Grid *g, int node, int edgeId );
Grid *gridProject(Grid *g);

Grid *gridSmooth(Grid *g );
Grid *gridSmoothNode(Grid *g, int node );
Grid *gridOptimizeT(Grid *g, int node, double dt );
Grid *gridOptimizeUV(Grid *g, int node, double *dudv );
Grid *gridOptimizeXYZ(Grid *g, int node, double *dxdydz );

END_C_DECLORATION

#endif /* GRIDCAD_H */
