
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDCAD_H
#define GRIDCAD_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridForceNodeToEdge(Grid *g, int node, int edgeId );
Grid *gridForceNodeToFace(Grid *g, int node, int faceId );

Grid *gridProjectNodeToEdge(Grid *g, int node, int edgeId );
Grid *gridProjectNodeToFace(Grid *g, int node, int faceId );

Grid *gridEvaluateEdgeAtT(Grid *g, int node, double t );
Grid *gridEvaluateFaceAtUV(Grid *g, int node, double *uv );
Grid *gridUpdateFaceParameter(Grid *g, int node );

Grid *gridProjectToFace(Grid *g, int faceId, 
			double *xyz, double *uv, double *newxyz );
Grid *gridFaceNormalAtUV(Grid *g, int faceId, 
			 double *uv, double *xyz, double *normal );
Grid *gridFaceNormalAtXYZ(Grid *g, int faceId, double *xyz, double *normal );

Grid *gridSafeProjectNode(Grid *g, int node, double ratio );
Grid *gridSafeProjectNodeToFace(Grid *g, int node, int faceId, double ratio );
Grid *gridSafeProjectNodeToEdge(Grid *g, int node, int edgeId, double ratio );

Grid *gridNodeProjectionDisplacement(Grid *g, int node, double *displacement );

Grid *gridProject(Grid *g);
Grid *gridRobustProjectNode(Grid *g, int node);
Grid *gridRobustProject(Grid *g);

Grid *gridSmooth(Grid *g );
Grid *gridSmoothFaceMR(Grid *g, double optimizationLimit );
Grid *gridSmoothVolume(Grid *g );
Grid *gridSmoothNearNode1(Grid *g, int node );
Grid *gridSmoothNearNode(Grid *g, int node );
Grid *gridSmoothNode(Grid *g, int node, GridBool smoothOnSurface );
Grid *gridSmoothNodeFaceMR(Grid *g, int node );

Grid *gridOptimizeT(Grid *g, int node, double dt );
Grid *gridOptimizeUV(Grid *g, int node, double *dudv );
Grid *gridOptimizeFaceUV(Grid *g, int node, double *dudv );
Grid *gridOptimizeXYZ(Grid *g, int node, double *dxdydz );

Grid *gridSmartLaplacian(Grid *g, int node );
Grid *gridSmartVolumeLaplacian(Grid *g, int node );

Grid *gridSmoothNodeQP(Grid *g, int node );

Grid *gridSmoothNodeVolume(Grid *g, int node );
Grid *gridSmoothNodeVolumeWithSurf(Grid *g, int node );
Grid *gridSmoothNodeVolumeSimplex(Grid *g, int node );

Grid *gridRelaxNegativeCells(Grid *g, GridBool dumpTecplot );
Grid *gridSmoothVolumeNearNode(Grid *grid, int node, 
			       GridBool smoothOnSurface );

GridBool nearestOnEdge(int vol, int edge, double *xyz, double *t,
                       double *xyznew);
GridBool nearestOnFace(int vol, int face, double *xyz, double *uv,
                       double *xyznew);

double gridFaceAreaUV(Grid *g, int face);
Grid *gridMinFaceAreaUV(Grid *g, int node, double *min_area);

Grid *gridSmoothNodeFaceAreaUV(Grid *g, int node );
Grid *gridSmoothNodeFaceAreaUVSimplex( Grid *g, int node );

END_C_DECLORATION

#endif /* GRIDCAD_H */
