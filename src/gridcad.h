
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
Grid *gridSafeProjectNode(Grid *g, int node, double ratio );
Grid *gridSafeProjectNodeToFace(Grid *g, int node, int faceId, double ratio );
Grid *gridSafeProjectNodeToEdge(Grid *g, int node, int edgeId, double ratio );
Grid *gridProject(Grid *g);
Grid *gridRobustProjectNode(Grid *g, int node);
Grid *gridRobustProject(Grid *g);

Grid *gridSmooth(Grid *g );
Grid *gridSmoothFaceMR(Grid *g );
Grid *gridSmoothVolume(Grid *g );
Grid *gridSmoothNode(Grid *g, int node );
Grid *gridSmoothNodeFaceMR(Grid *g, int node );
Grid *gridOptimizeT(Grid *g, int node, double dt );
Grid *gridOptimizeUV(Grid *g, int node, double *dudv );
Grid *gridOptimizeFaceUV(Grid *g, int node, double *dudv );
Grid *gridOptimizeXYZ(Grid *g, int node, double *dxdydz );
Grid *gridSmartLaplacian(Grid *g, int node );

END_C_DECLORATION

#endif /* GRIDCAD_H */
