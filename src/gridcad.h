
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

Grid *gridUpdateParameters(Grid *g, int node );
Grid *gridUpdateFaceParameters(Grid *g, int node );

Grid *gridProjectToEdge(Grid *g, int edgeId, 
			double *xyz, double *t, double *newxyz );
Grid *gridProjectToFace(Grid *g, int faceId, 
			double *xyz, double *uv, double *newxyz );
Grid *gridEvaluateOnEdge(Grid *g, int edgeId, double t, double *xyz );
Grid *gridEvaluateOnFace(Grid *g, int faceId, double *uv, double *xyz );
Grid *gridFaceNormalAtUV(Grid *g, int faceId, 
			 double *uv, double *xyz, double *normal );
Grid *gridFaceNormalAtXYZ(Grid *g, int faceId, double *xyz, double *normal );

Grid *gridProjectNode(Grid *g, int node );

Grid *gridNodeProjectionDisplacement(Grid *g, int node, double *displacement );

Grid *gridRobustProjectNode(Grid *g, int node);
Grid *gridRobustProject(Grid *g);

/* the {optimzation,lapacian}Limits will be set to a default if < 0.0 */ 
Grid *gridSmooth(Grid *g, double optimizationLimit, double laplacianLimit );
Grid *gridSmoothFaceMR(Grid *g, double optimizationLimit );
Grid *gridSmoothVolume(Grid *g );
Grid *gridSmoothNearNode1(Grid *g, int node );
Grid *gridSmoothNearNode(Grid *g, int node );
Grid *gridSmoothNode(Grid *g, int node, GridBool smoothOnSurface );
Grid *gridSmoothNodeFaceMR(Grid *g, int node );

#define gridOPTIM_COST_FLOOR (0.01)
Grid *gridLineSearchT(Grid *g, int node, double optimized_cost_limit );
Grid *gridLineSearchUV(Grid *g, int node, double *dudv, 
		       double optimized_cost_limit );
Grid *gridOptimizeFaceUV(Grid *g, int node, double *dudv );
Grid *gridLinearProgramUV(Grid *g, int node, GridBool *callAgain );
Grid *gridOptimizeXYZ(Grid *g, int node, double *dxdydz );

Grid *gridSmartLaplacian(Grid *g, int node );
Grid *gridSmartVolumeLaplacian(Grid *g, int node );

Grid *gridStoreVolumeCostDerivatives(Grid *g, int node );
Grid *gridStoreFaceCostParameterDerivatives(Grid *g, int node );

Grid *gridLinearProgramXYZ(Grid *g, int node, GridBool *callAgain );

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
