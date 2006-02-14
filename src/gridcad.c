
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <values.h>
#ifdef HAVE_SDK
#include "CADGeom/CADGeom.h"
#else
#include "FAKEGeom.h"
#endif
#include "plan.h"
#include "tableau.h"
#include "gridmath.h"
#include "gridmetric.h"
#include "gridshape.h"
#include "gridcad.h"

Grid *gridForceNodeToEdge(Grid *grid, int node, int edgeId )
{
  int vol = 1;
  double t, xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  t = DBL_MAX;

  if (!nearestOnEdge( vol, edgeId, xyz, &t, xyznew) ) return NULL;  

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;

  return grid;
}

Grid *gridForceNodeToFace(Grid *grid, int node, int faceId )
{
  int vol = 1;
  double uv[2], xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  uv[0] = DBL_MAX;
  uv[1] = DBL_MAX;

  if (!nearestOnFace( vol, faceId, xyz, uv, xyznew) ) {
    printf("%s: %d: nearestOnFace failed.\n",__FILE__,__LINE__);
    return NULL;  
  }

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;

  return grid;
}

Grid *gridProjectNodeToEdge(Grid *grid, int node, int edgeId )
{
  int vol = 1;
  double t, xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  if ( grid != gridNodeT( grid, node, edgeId, &t ) ) return NULL;

  if (!nearestOnEdge( vol, edgeId, xyz, &t, xyznew) ) return NULL;  

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;
  if ( grid != gridSetNodeT( grid, node, edgeId, t ) ) return NULL;

  return gridUpdateFaceParameters(grid,node);
}

Grid *gridProjectNodeToFace(Grid *grid, int node, int faceId )
{
  int vol = 1;
  double uv[2], xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  if ( grid != gridNodeUV( grid, node, faceId, uv ) ) return NULL;

  if (!nearestOnFace( vol, faceId, xyz, uv, xyznew) ) {
    printf("%s: %d: nearestOnFace failed.\n",__FILE__,__LINE__);
    return NULL;  
  }

  if ( grid != gridSetNodeXYZ( grid, node, xyznew ) ) return NULL;
  if ( grid != gridSetNodeUV( grid, node, faceId, uv[0], uv[1] ) ) return NULL;

  return gridUpdateFaceParameters(grid,node); /* only needed for FAKEGeom */
}

Grid *gridEvaluateEdgeAtT(Grid *grid, int node, double t )
{
  int vol =1;
  int nodes[2];
  int edge, edgeId;
  double xyz[3];

  if (gridGeometryNode(grid,node)) return NULL;
  if ( grid != gridNodeXYZ(grid, node, xyz) ) return NULL;

  edge = adjItem(adjFirst(gridEdgeAdj(grid), node));
  if ( grid != gridEdge(grid,edge,nodes,&edgeId)) return NULL;

  if ( !CADGeom_PointOnEdge( vol, edgeId, t, xyz, 
			     0, NULL, NULL) ) {
    printf ( "ERROR: CADGeom_PointOnEdge, %d: %s\n",__LINE__,__FILE__ );
    return NULL;
  }
  gridSetNodeT(grid, node, edgeId, t);
  gridSetNodeXYZ(grid, node, xyz);

  return gridUpdateFaceParameters(grid,node);
}

Grid *gridEvaluateFaceAtUV(Grid *grid, int node, double *uv )
{
  int vol =1;
  int nodes[3];
  int face, faceId;
  double xyz[3];

  if ( gridGeometryEdge(grid, node) ) return NULL;
  if ( grid != gridNodeXYZ(grid, node, xyz) ) return NULL;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if (EMPTY == face) return NULL;

  gridFace(grid,face,nodes,&faceId);

  if ( !CADGeom_PointOnFace( vol, faceId, uv, xyz, 
			     0, NULL, NULL, NULL, NULL, NULL) ){
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
    return NULL;
  }
  gridSetNodeUV(grid, node, faceId, uv[0], uv[1]);
  gridSetNodeXYZ(grid, node, xyz);
  return grid;
}

Grid *gridResolveTofEdge(Grid *grid, int node, int edgeId )
{
  int vol = 1;
  double t, xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  if ( grid != gridNodeT( grid, node, edgeId, &t ) ) return NULL;

  if (!CADGeom_ResolveOnEdgeWCS( vol, edgeId, xyz, &t, xyznew) ) {
    printf("%s: %d: CADGeom_ResolveOnEdgeWCS failed.\n",__FILE__,__LINE__);
    return NULL;
  }

  if ( grid != gridSetNodeT( grid, node, edgeId, t ) ) return NULL;

  return grid;
}

Grid *gridResolveUVofFace(Grid *grid, int node, int faceId )
{
  int vol = 1;
  double uv[2], xyz[3], xyznew[3];

  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  if ( grid != gridNodeUV( grid, node, faceId, uv ) ) return NULL;

  if (!CADGeom_ResolveOnFaceWCS( vol, faceId, xyz, uv, xyznew) ) {
    printf("%s: %d: CADGeom_ResolveOnFaceWCS failed.\n",__FILE__,__LINE__);
    return NULL;
  }

  if ( grid != gridSetNodeUV( grid, node, faceId, uv[0], uv[1] ) ) return NULL;

  return grid;
}

Grid *gridUpdateParameters(Grid *grid, int node ){

  AdjIterator it;
  int nodes[2];
  int edge, edgeId;

  for ( it = adjFirst(gridEdgeAdj(grid),node); 
	adjValid(it); 
	it = adjNext(it) ){
    edge = adjItem(it);
    if ( grid != gridEdge(grid, edge, nodes, &edgeId) ) return NULL;
    if ( grid != gridResolveTofEdge(grid, node, edgeId) ) return NULL;
  }
  return gridUpdateFaceParameters(grid, node );
}

Grid *gridUpdateFaceParameters(Grid *grid, int node ){

  AdjIterator it;
  int nodes[3];
  int face, faceId;

  for ( it = adjFirst(gridFaceAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ){
    face = adjItem(it);
    if ( grid != gridFace(grid, face, nodes, &faceId) )  return NULL;
    if ( grid != gridResolveUVofFace( grid, node, faceId ) ) return NULL;
  }

  return grid;
}

Grid *gridProjectToEdge(Grid *grid, int edgeId, 
			double *xyz, double *t, double *newxyz )
{
  int vol = 1;

  if (!nearestOnEdge( vol, edgeId, xyz, t, newxyz) ) {
    printf("%s: %d: nearestOnEdge failed.\n",__FILE__,__LINE__);
    return NULL;  
  }

  return grid;
}

Grid *gridProjectToFace(Grid *grid, int faceId, 
			double *xyz, double *uv, double *newxyz )
{
  int vol = 1;

  if (!nearestOnFace( vol, faceId, xyz, uv, newxyz) ) {
    printf("%s: %d: nearestOnFace failed.\n",__FILE__,__LINE__);
    return NULL;  
  }

  return grid;
}

Grid *gridEvaluateOnEdge(Grid *grid, int edgeId, double t, double *xyz )
{
  int vol = 1;

  if (!CADGeom_PointOnEdge( vol, edgeId, t, xyz, 0, NULL, NULL) ) {
    printf("%s: %d: CADGeom_PointOnEdge( failed.\n",__FILE__,__LINE__);
    return NULL;  
  }

  return grid;
}

Grid *gridEvaluateOnFace(Grid *grid, int faceId, double *uv, double *xyz )
{
  int vol = 1;

  if (!CADGeom_PointOnFace( vol, faceId, uv, xyz, 
			     0, NULL, NULL, NULL, NULL, NULL) ) {
    printf("%s: %d: CADGeom_PointOnFace failed.\n",__FILE__,__LINE__);
    return NULL;  
  }

  return grid;
}

Grid *gridResolveOnFace(Grid *grid, int faceId,
			double *uv, double *original_xyz, double *resolved_xyz )
{
  int vol = 1;

  if (!CADGeom_ResolveOnFaceWCS(vol, faceId, original_xyz, uv, resolved_xyz) ) {
    printf("%s: %d: CADGeom_ResolveOnFaceWCS failed.\n",__FILE__,__LINE__);
    return NULL;  
  }

  return grid;
}

Grid *gridFaceNormalAtUV(Grid *grid, int faceId,
			 double *uv, double *xyz, double *normal )

{
  int vol =1;

  if ( !CADGeom_NormalToFace( vol, faceId, uv, xyz, normal) ){
    printf ( "ERROR: CADGeom_NormalToFace, %d: %s\n",__LINE__,__FILE__ );
    return NULL;
  }
  return grid;
}

Grid *gridFaceNormalAtXYZ(Grid *grid, int faceId, double *xyz, double *normal )
{
  double uv[2], newxyz[3];

  uv[0] = DBL_MAX;
  uv[1] = DBL_MAX;
  if (grid != gridProjectToFace(grid,faceId,xyz,uv,newxyz)) return NULL;
  if (grid != gridFaceNormalAtUV(grid, faceId, uv, newxyz, normal)) return NULL;

  return grid;
}

Grid *gridProjectNode(Grid *grid, int node )
{
  int nodes[3];
  int edge, edgeId;
  int face, faceId;

  if ( gridGeometryNode( grid, node ) ) { 
    return gridUpdateParameters(grid,node);
  }
  if ( gridGeometryEdge( grid, node ) ) {
    edge = adjItem(adjFirst(gridEdgeAdj(grid), node));
    gridEdge(grid, edge, nodes, &edgeId );
    return gridProjectNodeToEdge( grid, node, edgeId );
  }
  if ( gridGeometryFace( grid, node ) ) {
    face = adjItem(adjFirst(gridFaceAdj(grid), node));
    gridFace(grid, face, nodes, &faceId );
    return gridProjectNodeToFace( grid, node, faceId );
  }

  return grid;
}

Grid *gridNodeProjectionDisplacement(Grid *grid, int node,
				     double *displacement )
{
  int nodes[3];
  int edge, edgeId;
  int face, faceId;
  double t, uv[2], xyz[3], xyznew[3];
  GridBool evaluate = TRUE;

  displacement[0] = displacement[1] = displacement[2] = 0.0;
  
  if (!gridGeometryFace( grid, node ) ) return grid;
  
  if ( grid != gridNodeXYZ( grid, node, xyz ) ) return NULL;
  
  if (evaluate) {
    if ( gridGeometryNode( grid, node ) ) { 
    } else {
      if ( gridGeometryEdge( grid, node ) ) {
	edge = adjItem(adjFirst(gridEdgeAdj(grid), node));
	gridEdge(grid, edge, nodes, &edgeId );
	gridNodeT(grid, node, edgeId, &t);
	gridEvaluateEdgeAtT(grid, node, t );
      } else {
	face = adjItem(adjFirst(gridFaceAdj(grid), node));
	gridFace(grid, face, nodes, &faceId );
	gridNodeUV(grid, node, faceId, uv);
	gridEvaluateFaceAtUV(grid, node, uv );
      }
    }
  }else{
    gridProjectNode(grid, node);
  }
  gridNodeXYZ( grid, node, xyznew );
  gridSubtractVector(xyznew,xyz,displacement);
  gridSetNodeXYZ( grid, node, xyz );
  return grid;
}

Grid *gridRobustProject(Grid *grid)
{
  int  node;
  int notProjected;
  notProjected = 0;
  for (node=0;node<gridMaxNode(grid);node++)
    if ( gridValidNode( grid, node ) && !gridNodeFrozen(grid, node) ) 
      if (grid != gridProjectNode( grid, node) ) 
	notProjected++;

  if (notProjected > 0){
    printf("gridRobustProject: %d of %d nodes not projected.\n",
	   notProjected,gridNNode(grid));
    return NULL;
  }

  return grid;
}

Grid *gridUntangle(Grid *grid)
{
  double allowedVolume, allowedArea;
  double minVolume, minArea;
  double last_volume;
  int count;
  double volume, area;
  int fix_node;
  int active_nodes;
  int tries;
  int stalled;

  allowedArea = 1.0e-12;

  minArea = gridMinGridFaceAreaUV(grid);
  tries = 0;
  while ( minArea <= allowedArea ) {
    tries++;
    if (tries >50) {
      printf("unable to fix min face UV area %e\n",minArea);
      return NULL;
    }
    if (tries >3) 
      printf("untangle faces %2d, min face UV area %e\n",tries,minArea);
    active_nodes = 0;
    for( fix_node=0; fix_node < gridMaxNode(grid);fix_node++) {
      if ( gridValidNode( grid, fix_node ) &&
	   !gridNodeFrozen( grid, fix_node) &&
	   gridGeometryFace(grid, fix_node) &&
	   !gridGeometryEdge(grid, fix_node) ) {
	gridMinFaceAreaUV(grid, fix_node, &area);
	if ( (area <= allowedArea) || (tries>10) ) {
	  active_nodes++;
	  if ( grid != gridUntangleAreaUV( grid, fix_node, 1 ) ) {
	    printf( "%s: %d: %s: gridUntangleAreaUV NULL\n",
		    __FILE__, __LINE__, "gridUntangle");
	    return NULL;
	  }
	}
      }
    }
    minArea = gridMinGridFaceAreaUV(grid);
  }

  allowedVolume = 1.0e-14;

  gridMinVolumeAndCount( grid, &minVolume, &count );
  tries = 0;
  stalled = 0;
  active_nodes=0;
  while ( minVolume <= allowedVolume ) {
    tries++;
    if (tries >20) {
      printf("unable to fix min volume %e\n",minVolume);
      return NULL;
    }
    if (tries >4)
      printf("untangle tets %3d, min volume %25.15e of %4d %4d\n",
	     tries,minVolume,count,active_nodes);
    if (tries >3) { 
      gridCollapseWedgeCells(grid);
      gridMinVolumeAndCount( grid, &minVolume, &count );
    }
    active_nodes = 0;
    for( fix_node=0; fix_node < gridMaxNode(grid);fix_node++) {
      if ( gridValidNode( grid, fix_node ) &&
	   !gridNodeFrozen( grid, fix_node) &&
	   !gridGeometryFace(grid, fix_node) ) {
	gridNodeVolume(grid, fix_node, &volume );
	if ( (volume <= allowedVolume) ) {
	  active_nodes++;
	  if ( grid != gridUntangleVolume( grid, fix_node, 2 ) ) {
	    printf( "%s: %d: %s: gridUntangleVolume NULL\n",
		    __FILE__, __LINE__, "gridUntangle");	    
	    return NULL;
	  }
	}
      }
    }
    last_volume = minVolume;
    gridMinVolumeAndCount( grid, &minVolume, &count );
    if (!(last_volume < minVolume)) {
      stalled++;
    }else{
      stalled=0;
    }
    if (stalled >= 5) {
      printf("unable to make additional headway\n");
      return NULL;
    }
  }
  if (tries >4)
    printf("untangle tets %3d, min volume %25.15e of %4d %4d\n",
	   tries,minVolume,count,active_nodes);
  return grid;
}

Grid *gridSequentialEvaluation(Grid *grid)
{
  int node;
  
  double displacement[3];
  double original_xyz[3],projected_xyz[3];

  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) && 
	 !gridNodeFrozen( grid, node) &&
	 gridGeometryFace(grid,node) ) {
      gridNodeXYZ(grid,node,original_xyz);
      gridNodeProjectionDisplacement( grid, node, displacement );
      projected_xyz[0] = original_xyz[0] + displacement[0];
      projected_xyz[1] = original_xyz[1] + displacement[1];
      projected_xyz[2] = original_xyz[2] + displacement[2];
      gridSetNodeXYZ(grid,node,projected_xyz);

      if (grid != gridUntangle(grid)) return NULL;
    }
  }


  return grid;
}

Grid *gridWholesaleEvaluation(Grid *grid)
{
  int node;

  double displacement[3];
  double original_xyz[3],projected_xyz[3];

  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) && 
	 !gridNodeFrozen( grid, node) &&
	 gridGeometryFace(grid,node) ) {
      gridNodeXYZ(grid,node,original_xyz);
      gridNodeProjectionDisplacement( grid, node, displacement );
      projected_xyz[0] = original_xyz[0] + displacement[0];
      projected_xyz[1] = original_xyz[1] + displacement[1];
      projected_xyz[2] = original_xyz[2] + displacement[2];
      gridSetNodeXYZ(grid,node,projected_xyz);
    }
  }

  if (grid != gridUntangle(grid)) return NULL;

  return grid;
}

GridBool gridNewGeometryEdgeSiteAllowedAt(Grid *grid, 
					  int node0, int node1, double t )
{
  int edge, nodes[2], edgeId;
  double t0, t1;
  double dt, ratio;
  double xyz[3];
  int face0, face1;
  int faceId0, faceId1;
  int nodes0[3], nodes1[3];
  double uv0[3][2], uv1[3][2];
  double interpolated_uv0[2], interpolated_uv1[2];
  double projected_uv0[2], projected_uv1[2];
  double xyz0[3], xyz1[3];
  double uv[3][2];
  
  double du, dv, distance_squared;

  int i;

  /* require a geometry edge (segment) for support */
  edge = gridFindEdge(grid, node0, node1 );
  if ( EMPTY == edge ) return FALSE;
  gridEdge(grid, edge, nodes, &edgeId );

  /* ensure monotonicity of edge */
  gridNodeT( grid, node0, edgeId, &t0);
  gridNodeT( grid, node1, edgeId, &t1);
  if ( t < MIN(t0,t1) || MAX(t0,t1) < t ) return FALSE;

  /* test for collapsed edge and compute ratio of end points for new t */
  dt = t1 - t0;
  if ( ABS(dt) < 1.0e-12 ) return FALSE;
  ratio = (t - t0) / (t1 - t0); /* t = r*t1 + (1-r)*t0 */

  /* make sure that CAD can sucessfully evaluate at this point */
  if ( grid != gridEvaluateOnEdge( grid, edgeId, t, xyz ) ) return FALSE;

  /* require a triangle on each side of edge */
  face0 = gridFindFaceWithNodesUnless( grid, node0, node1, EMPTY );
  if (EMPTY == face0) return FALSE;
  face1 = gridFindFaceWithNodesUnless( grid, node0, node1, face0 );
  if (EMPTY == face1) return FALSE;
  gridFace( grid, face0, nodes0, &faceId0 );
  gridFace( grid, face1, nodes1, &faceId1 );

  /* and obtain UV values at triangle nodes */
  for ( i=0 ; i < 3 ; i++ ) {
    gridNodeUV( grid, nodes0[i], faceId0, uv0[i] );
    gridNodeUV( grid, nodes1[i], faceId1, uv1[i] );
  }

  /* and obtain a linear guess for new UV location uv = (1-r)*uv + r*uv */
  interpolated_uv0[0] = interpolated_uv0[1] = 0.0;
  for ( i=0 ; i < 3 ; i++ ) {
    if ( nodes0[i] == node0 ) {
      interpolated_uv0[0] += (1-ratio)*uv0[i][0];
      interpolated_uv0[1] += (1-ratio)*uv0[i][1];
    }
    if ( nodes0[i] == node1 ) {
      interpolated_uv0[0] += ratio*uv0[i][0];
      interpolated_uv0[1] += ratio*uv0[i][1];
    }
  }
  interpolated_uv1[0] = interpolated_uv1[1] = 0.0;
  for ( i=0 ; i < 3 ; i++ ) {
    if ( nodes1[i] == node0 ) {
      interpolated_uv1[0] += (1-ratio)*uv1[i][0];
      interpolated_uv1[1] += (1-ratio)*uv1[i][1];
    }
    if ( nodes1[i] == node1 ) {
      interpolated_uv1[0] += ratio*uv1[i][0];
      interpolated_uv1[1] += ratio*uv1[i][1];
    }
  }

  /* make sure that the geometry faces can be projected to */
  projected_uv0[0] = interpolated_uv0[0];
  projected_uv0[1] = interpolated_uv0[1];
  if ( grid != gridProjectToFace( grid, faceId0, xyz, projected_uv0, xyz0 ) ) 
    return FALSE;
  projected_uv1[0] = interpolated_uv1[0];
  projected_uv1[1] = interpolated_uv1[1];
  if ( grid != gridProjectToFace( grid, faceId1, xyz, projected_uv1, xyz1 ) ) 
    return FALSE;

  /* if the t value is near an edge segment endpoint make sure the uv
     is close to the existing trangle uv value */
  if ( ratio < 1.0e-7 ){
    du = projected_uv0[0]-interpolated_uv0[0];
    dv = projected_uv0[1]-interpolated_uv0[1];
    distance_squared = du*du + dv*dv;
    if ( distance_squared > 1.0e-12 ) return FALSE;
    return TRUE;
  }
  if ( ratio > (1.0-1.0e-7) ) {
    du = projected_uv1[0]-interpolated_uv1[0];
    dv = projected_uv1[1]-interpolated_uv1[1];
    distance_squared = du*du + dv*dv;
    if ( distance_squared > 1.0e-12 ) return FALSE;
    return TRUE;
  }

  /* make sure that projected uv values produce positve uv area triangles 
     for the middle of the edges */
  for ( i=0 ; i < 3 ; i++ ) {
    if ( nodes0[i] == node0 ) {
      uv[i][0] = projected_uv0[0];
      uv[i][1] = projected_uv0[1];
    } else {
      uv[i][0] = uv0[i][0];
      uv[i][1] = uv0[i][1];
    }
  }
  if ( gridFaceAreaUVDirect(grid, uv[0], uv[1], uv[2], faceId0 ) < 1.0e-12 )
    return FALSE;

  for ( i=0 ; i < 3 ; i++ ) {
    if ( nodes0[i] == node1 ) {
      uv[i][0] = projected_uv0[0];
      uv[i][1] = projected_uv0[1];
    } else {
      uv[i][0] = uv0[i][0];
      uv[i][1] = uv0[i][1];
    }
  }
  if ( gridFaceAreaUVDirect(grid, uv[0], uv[1], uv[2], faceId0 ) < 1.0e-12 )
    return FALSE;

  for ( i=0 ; i < 3 ; i++ ) {
    if ( nodes1[i] == node0 ) {
      uv[i][0] = projected_uv1[0];
      uv[i][1] = projected_uv1[1];
    } else {
      uv[i][0] = uv1[i][0];
      uv[i][1] = uv1[i][1];
    }
  }
  if ( gridFaceAreaUVDirect(grid, uv[0], uv[1], uv[2], faceId1 ) < 1.0e-12 )
    return FALSE;

  for ( i=0 ; i < 3 ; i++ ) {
    if ( nodes1[i] == node1 ) {
      uv[i][0] = projected_uv1[0];
      uv[i][1] = projected_uv1[1];
    } else {
      uv[i][0] = uv1[i][0];
      uv[i][1] = uv1[i][1];
    }
  }
  if ( gridFaceAreaUVDirect(grid, uv[0], uv[1], uv[2], faceId1 ) < 1.0e-12 )
    return FALSE;

  return TRUE;
}

Grid *gridCurveIntersectsFace(Grid *grid, int *face_nodes, int parent,
			      double *tuv0_start, double *tuv1_start,
			      double *tuv, double *curve, double *bary )
{
  double xyz0[3], xyz1[3], xyz2[3];
  double edge0[3], edge1[3];
  double norm[3];
  double tuv0[2], tuv1[2], ratio;
  GridBool keep_going;
  double curve0[3], curve1[3];
  double dir0[3], dir1[3], dir[3];
  double dot0, dot1, dot;

  if (grid != gridNodeXYZ(grid, face_nodes[0], xyz0) ) return NULL;
  if (grid != gridNodeXYZ(grid, face_nodes[1], xyz1) ) return NULL;
  if (grid != gridNodeXYZ(grid, face_nodes[2], xyz2) ) return NULL;

  gridSubtractVector(xyz1, xyz0, edge0);
  gridSubtractVector(xyz2, xyz0, edge1);
  gridCrossProduct(edge0,edge1,norm);
  gridVectorNormalize(norm);


  tuv0[0] = tuv0_start[0]; tuv0[1] = tuv0_start[1]; 
  tuv1[0] = tuv1_start[0]; tuv1[1] = tuv1_start[1]; 
  if (parent > 0) {
    gridEvaluateOnFace(grid, parent, tuv0, curve0 );
    gridEvaluateOnFace(grid, parent, tuv1, curve1 );
  }else{
    gridEvaluateOnEdge(grid, -parent, tuv0[0], curve0 );
    gridEvaluateOnEdge(grid, -parent, tuv1[0], curve1 );
  }
  gridSubtractVector(curve0, xyz0, dir0);
  dot0 = gridDotProduct(dir0,norm);
  gridSubtractVector(curve1, xyz0, dir1);
  dot1 = gridDotProduct(dir1,norm);

  tuv[0] = tuv[1] = DBL_MAX;
  ratio = 0.5;

  keep_going = TRUE;
  while (keep_going) {

    if (parent > 0) {
      tuv[0] = (1-ratio)*tuv0[0]+ratio*tuv1[0];
      tuv[1] = (1-ratio)*tuv0[1]+ratio*tuv1[1];
      gridEvaluateOnFace(grid, parent, tuv, curve );
    }else{
      tuv[0] = (1-ratio)*tuv0[0]+ratio*tuv1[0];
      gridEvaluateOnEdge(grid, -parent, tuv[0], curve );
    }
    gridSubtractVector(curve, xyz0, dir);
    dot = gridDotProduct(dir,norm);
    
    //printf("ratio%11.8f dots%23.15e%23.15e%23.15e\n",ratio,dot0,dot,dot1);

    if ( dot0 < 0.0 || dot1 > 0.0 ) return NULL;
    
    if (dot < 0.0) {
      dot1 = dot;
      tuv1[0] = tuv[0]; tuv1[1] = tuv[1];
    }else{
      dot0 = dot;
      tuv0[0] = tuv[0]; tuv0[1] = tuv[1];
    }

    if (ABS(dot) <=10.e-10 || ABS(dot1-dot0) <= 1.0e-10){
      ratio = 0.5;
      keep_going = FALSE;
    }else{
      ratio = dot0 / (dot0-dot1);
    }

  }
  
  gridBarycentricCoordinateTri(xyz0,xyz1,xyz2,curve,bary);

  //printf("bary %23.15e%23.15e%23.15e\n",bary[0],bary[1],bary[2]);

  if ( bary[0] < 0.0 || bary[1] < 0.0 || bary[2] < 0.0 ) {
    return NULL;
  }else{
    return grid;
  }
}

Grid *gridLineSegmentIntersectsFace(Grid *grid, int node0, int node1,
				    int faceId, double *uv_guess,
				    double *uv, double *xyz, double *bary )
{
  double xyz0[3], xyz1[3];
  double edge_direction[3], length;
  double nearest_point_on_edge[3];
  double disp[3];
  double normal_displacement[3], radius;
  double last_bary;
  int iteration;
  GridBool keep_going;

  if (grid != gridNodeXYZ(grid, node0, xyz0) ) return NULL;
  if (grid != gridNodeXYZ(grid, node1, xyz1) ) return NULL;

  gridSubtractVector(xyz1, xyz0, edge_direction);
  length = gridVectorLength(edge_direction);
  gridVectorNormalize(edge_direction);
  
  (*bary) = 0.5;
  last_bary = DBL_MAX;

  nearest_point_on_edge[0] = (*bary)*xyz1[0] + (1.0-(*bary))*xyz0[0];
  nearest_point_on_edge[1] = (*bary)*xyz1[1] + (1.0-(*bary))*xyz0[1];
  nearest_point_on_edge[2] = (*bary)*xyz1[2] + (1.0-(*bary))*xyz0[2];

  iteration = 0;
  keep_going = TRUE;
  while (keep_going) {
    iteration++; 
    if (iteration > 10000) {
      printf("iterations %d exhasted bary%12.8f diff %e\n",
	     iteration,(*bary),radius);
      return NULL;
    }
    gridProjectToFace(grid, faceId, nearest_point_on_edge, uv, xyz);

    gridSubtractVector(xyz, xyz0, disp);
    (*bary) = gridDotProduct( edge_direction, disp )/length;

    nearest_point_on_edge[0] = (*bary)*xyz1[0] + (1.0-(*bary))*xyz0[0];
    nearest_point_on_edge[1] = (*bary)*xyz1[1] + (1.0-(*bary))*xyz0[1];
    nearest_point_on_edge[2] = (*bary)*xyz1[2] + (1.0-(*bary))*xyz0[2];

    gridSubtractVector(xyz, nearest_point_on_edge, normal_displacement);
    radius = sqrt( gridDotProduct( normal_displacement, normal_displacement ) );

    if ( (*bary) > 1.1 || (*bary) < -0.1 ) keep_going = FALSE;
    if ( ABS((*bary)-last_bary) < 1.0e-8 ) keep_going = FALSE;
    last_bary = (*bary);

  }
  printf("bary%12.8f diff %e iter %d\n",(*bary),radius,iteration);
  if (radius < 1.0e-6 && (*bary) <= 1.0 && (*bary) >= 0.0) {
    return grid;
  }else{
    return NULL;
  }

}

Grid *gridSmoothNearNode1(Grid *grid, int node )
{
#define SMOOTHDEG (500)
  int i, nodes[4];
  int nlist, smooth, look, nodelist[SMOOTHDEG];
  GridBool looking;
  AdjIterator it;

  if (!gridValidNode(grid,node)) return NULL;
  
  nlist =0;
  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it); 
	it = adjNext(it) ){
    gridCell(grid, adjItem(it), nodes);
    for (i=0;i<4;i++) {
      if (!gridNodeFrozen( grid, nodes[i])) {
	looking = (nlist<=SMOOTHDEG);
	look = 0;
	for (look=0;look<nlist && looking ; look++){
	  looking = (nodelist[look] != nodes[i]);
	}
	if (looking && nlist<=SMOOTHDEG){
	  nodelist[nlist] = nodes[i];
	  nlist++;
	}
      }
    }
  }      

  for (smooth=0;smooth<5;smooth++)
    for (i=0;i<nlist;i++) 
      gridSmoothNode( grid, nodelist[i], TRUE);

  return grid;
}

Grid *gridSmoothNearNode(Grid *grid, int node )
{
  int i, i0, nodes[4], nodes0[4];
  int nlist, smooth, look, nodelist[SMOOTHDEG];
  GridBool looking;
  AdjIterator it, it0;

  if (!gridValidNode(grid,node)) return NULL;

  nlist =0;
  for ( it0 = adjFirst(gridCellAdj(grid),node); 
	adjValid(it0); 
	it0 = adjNext(it0) ){
    gridCell(grid, adjItem(it0), nodes0);
    for (i0=0;i0<4;i0++) {
      for ( it = adjFirst(gridCellAdj(grid),nodes0[i0]); 
	    adjValid(it); 
	    it = adjNext(it) ){
	gridCell(grid, adjItem(it), nodes);
	for (i=0;i<4;i++) {
	  if (!gridNodeFrozen( grid, nodes[i])) {
	    looking = (nlist<=SMOOTHDEG);
	    look = 0;
	    for (look=0;look<nlist && looking ; look++){
	      looking = (nodelist[look] != nodes[i]);
	    }
	    if (looking && nlist<=SMOOTHDEG){
	      nodelist[nlist] = nodes[i];
	      nlist++;
	    }
	  }
	}
      }
    }
  }      

  for (smooth=0;smooth<2;smooth++)
    for (i=0;i<nlist;i++) 
      gridSmoothNode( grid, nodelist[i], TRUE);

  return grid;
}

Grid *gridSmoothNode(Grid *grid, int node, GridBool smoothOnSurface )
{
  double xyzProj[3], uv[2];
  double ar, dARdx[3];
  double mr, dMRdx[3];
  double du[3], dv[3];
  double dMRdu[2];
  int vol =1;
  int nodes[3];
  int face, faceId;

  int maxsmooth;
  GridBool callAgain;

  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridGeometryBetweenFace( grid, node ) &&
       !gridGeometryEdge( grid, node ) ) return grid;

  /* skip boundary nodes if we have been asked not to smooth on surface */ 
  if ( gridGeometryFace( grid, node ) && !smoothOnSurface ) return grid;

  /* edge smooth */
  if ( gridGeometryEdge( grid, node ) ) {
    return gridLineSearchT(grid, node, gridMinSurfaceSmoothCost(grid) );
  }

  /* face smooth */
  if ( gridGeometryFace( grid, node ) ) {
    if (FALSE) {
      for (maxsmooth=0;maxsmooth<3;maxsmooth++) {
	face = adjItem(adjFirst(gridFaceAdj(grid), node));
	gridFace(grid,face,nodes,&faceId);
	gridNodeFaceMRDerivative ( grid, node, &mr, dMRdx);
	gridNodeUV( grid, node, faceId, uv);
	if ( !CADGeom_PointOnFace( vol, faceId,   
				   uv, xyzProj, 1, du, dv, NULL, NULL, NULL) )
	  printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
       
	dMRdu[0] = dMRdx[0]*du[0] + dMRdx[1]*du[1] + dMRdx[2]*du[2] ; 
	dMRdu[1] = dMRdx[0]*dv[0] + dMRdx[1]*dv[1] + dMRdx[2]*dv[2] ; 
	if (grid != gridLineSearchUV( grid, node, dMRdu, 
				      gridMinSurfaceSmoothCost(grid) ) ) 
	  return NULL;
      }
      return grid;
    }else{
      maxsmooth = 40;
      callAgain = TRUE;
      while ( callAgain && maxsmooth > 0 ) {
	maxsmooth--;
	if (grid != gridLinearProgramUV( grid, node, &callAgain ) ) {
	  return NULL;
	}
      }
      return grid;
    }
  }

  /* volume node smooth */
  if (FALSE) {
    gridNodeARDerivative ( grid, node, &ar, dARdx);
    return gridOptimizeXYZ( grid, node, dARdx );
  }else{
    maxsmooth = 40;
    callAgain = TRUE;
    while ( callAgain && maxsmooth > 0 ) {
      maxsmooth--;
      if (grid != gridLinearProgramXYZ( grid, node, &callAgain ) ) {
	return NULL;
      }
    }  
    return grid;
  }
}


Grid *gridSmoothNodeFaceMR(Grid *grid, int node )
{
  double xyzProj[3], uv[2];
  double mr, dMRdx[3];
  double du[3], dv[3];
  double dMRdu[2];
  int vol =1;
  int nodes[3];
  int face, faceId;

  if ( gridGeometryEdge( grid, node ) ) return grid;
  if ( !gridGeometryFace( grid, node ) ) return grid;


  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  gridFace(grid,face,nodes,&faceId);
  gridNodeFaceMRDerivative ( grid, node, &mr, dMRdx);
  gridNodeUV( grid, node, faceId, uv);
  if ( !CADGeom_PointOnFace( vol, faceId,   
			     uv, xyzProj, 1, du, dv, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );

  dMRdu[0] = dMRdx[0]*du[0] + dMRdx[1]*du[1] + dMRdx[2]*du[2] ; 
  dMRdu[1] = dMRdx[0]*dv[0] + dMRdx[1]*dv[1] + dMRdx[2]*dv[2] ; 

  return gridOptimizeFaceUV( grid, node, dMRdu );
}

Grid *gridLineSearchT(Grid *grid, int node, double optimized_cost_limit )
{
  double dt, t;
  double gold;
  double alpha[2], ar[2], equality[2];
  int iter;

  int nodes1[2], nodes2[2];
  int edge1, edge2, edgeId1, edgeId2;
  int node1, node2, nodeTemp;
  double ratio1, ratio2;
  double tStart, tEnd;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  /* do not allow geometry or ghost nodes to move */
  if ( gridGeometryNode( grid, node ) ) return grid;
  if ( gridNodeGhost(grid, node ) ) return NULL;

  /* get the two edges involved (incident to node) */
  edge1 = adjItem(adjFirst(gridEdgeAdj(grid), node));
  if (EMPTY == edge1) return NULL;
  edge2 = adjItem(adjNext(adjFirst(gridEdgeAdj(grid), node)));
  if (EMPTY == edge2) return NULL;

  /* make sure the we are on a single edge and get the adjacent nodes */
  if ( grid != gridEdge(grid,edge1,nodes1,&edgeId1) ) return NULL;
  if ( grid != gridEdge(grid,edge2,nodes2,&edgeId2) ) return NULL;
  if ( edgeId1 != edgeId2 ) return NULL;
  node1 = nodes1[0] + nodes1[1] - node; 
  node2 = nodes2[0] + nodes2[1] - node; 

  /* get lengths in mapped space */
  ratio1 = gridEdgeRatio(grid, node, node1);
  ratio2 = gridEdgeRatio(grid, node, node2);

  /* if necessary, switch node1 and node2 so we are moving toward node1 */
  if ( ratio2 > ratio1 ) {
    nodeTemp = node1; node1 = node2; node2 = nodeTemp;
  }

  /* calculate dt to point in the parameter direction of `longer' edge */
  if ( grid != gridNodeT( grid, node, edgeId1, &tStart ) ) return NULL;
  if ( grid != gridNodeT( grid, node1, edgeId1, &tEnd ) ) return NULL;
  dt = tEnd - tStart;

  /* initialize the alpha, ar, and equality arrays */
  alpha[0] = 0.0;
  t = tStart + alpha[0]*dt;
  if (grid != gridEvaluateEdgeAtT(grid, node, t ) ) return NULL;
  gridNodeAR( grid, node, &ar[0] );
  ratio1 = gridEdgeRatio(grid, node, node1);
  ratio2 = gridEdgeRatio(grid, node, node2);
  equality[0] = ratio2/ratio1;

  alpha[1] = 1.0e-10;
  t = tStart + alpha[1]*dt;
  if (grid != gridEvaluateEdgeAtT(grid, node, t ) ) return NULL;
  gridNodeAR( grid, node, &ar[1] );
  ratio1 = gridEdgeRatio(grid, node, node1);
  ratio2 = gridEdgeRatio(grid, node, node2);
  equality[1] = ratio2/ratio1;

  /* try larger alphas if equality is improving and ar is valid */
  iter = 0;
  while ( equality[1] > equality[0] && 
	  equality[1] < 1.0 &&
	  ar[1] > optimized_cost_limit && 
	  iter < 100 ) {
    iter++;
    alpha[0] = alpha[1]; ar[0] = ar[1]; equality[0] =  equality[1];
    alpha[1] = alpha[0] * gold;
    t = tStart + alpha[1]*dt;
    if (grid != gridEvaluateEdgeAtT(grid, node, t ) ) return NULL;
    gridNodeAR( grid, node, &ar[1] );
    ratio1 = gridEdgeRatio(grid, node, node1);
    ratio2 = gridEdgeRatio(grid, node, node2);
    equality[1] = ratio2/ratio1;
  }

  /* use the `best' t and update node xyz, edge t, and face uv */
  t = tStart + alpha[0]*dt;
  if (grid != gridEvaluateEdgeAtT(grid, node, t ) ) return NULL;

  return grid;
}

Grid *gridLineSearchUV(Grid *grid, int node, double *dudv,
		       double optimized_cost_limit )
{
  double uvOrig[2], uv[2];
  int nodes[3];
  int face, faceId;
  double gold;
  double alpha[2], ar[2], mr[2];
  int iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  gridNodeUV( grid, node, faceId, uvOrig);

  alpha[0] = 0.0;
  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
  gridNodeAR( grid, node, &ar[0] );
  gridNodeFaceMR( grid, node, &mr[0] );

  alpha[1] = 1.0e-10;
  uv[0] = uvOrig[0] + alpha[1]*dudv[0];
  uv[1] = uvOrig[1] + alpha[1]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
  gridNodeAR( grid, node, &ar[1] );
  gridNodeFaceMR( grid, node, &mr[1] );

  iter = 0;
  while ( mr[1] > mr[0] && ar[1] > optimized_cost_limit && iter < 100){
    iter++;
    alpha[0] = alpha[1]; ar[0] = ar[1]; mr[0] = mr[1];
    alpha[1] = alpha[0] * gold;
    uv[0] = uvOrig[0] + alpha[1]*dudv[0];
    uv[1] = uvOrig[1] + alpha[1]*dudv[1];
    if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
    gridNodeAR( grid, node, &ar[1] );
    gridNodeFaceMR( grid, node, &mr[1] );
  }

  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;

//printf("node %d alpha %e ar %f uv %e %e\n",node,alpha[0],ar[0],uv[0],uv[1]);
  
  return grid;
}

Grid *gridOptimizeFaceUV(Grid *grid, int node, double *dudv )
{
  double uvOrig[2], uv[2];
  int nodes[3];
  int face, faceId;
  double gold;
  double alpha[2], mr[2], ar;
  int iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  gridNodeUV( grid, node, faceId, uvOrig);

  alpha[0] = 0.0;
  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  
  gridNodeFaceMR( grid, node, &mr[0] );

  alpha[1] = 1.0e-10;
  uv[0] = uvOrig[0] + alpha[1]*dudv[0];
  uv[1] = uvOrig[1] + alpha[1]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
  gridNodeFaceMR( grid, node, &mr[1] );

  iter = 0;
  gridNodeAR( grid, node, &ar ); 
  while ( mr[1] > mr[0] && mr[1] > 0.0 && iter < 100 && ar > 0.1){
    iter++;
    alpha[0] = alpha[1]; mr[0] = mr[1];
    alpha[1] = alpha[0] * gold;
    uv[0] = uvOrig[0] + alpha[1]*dudv[0];
    uv[1] = uvOrig[1] + alpha[1]*dudv[1];
    if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
    gridNodeFaceMR( grid, node, &mr[1] );
    gridNodeAR( grid, node, &ar );
  }

  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;

//printf("node %d alpha %e mr %f uv %e %e\n",node,alpha[0],mr[0],uv[0],uv[1]);
  
  return grid;
}

Grid *gridLinearProgramUV(Grid *grid, int node, GridBool *callAgain )
{
  int i, minFace, nearestFace;
  double minCost, nearestCost, nearestDifference, newCost, searchDirection[3];
  double g00, g01, g11, minRatio, nearestRatio;
  double length, projection;
  double deltaCost, currentAlpha, alpha, lastAlpha;
  double predictedImprovement, actualImprovement, lastImprovement;
  double minDirection[3], nearestDirection[3], dCostdX[3];
  int face, faceId, nodes[3];
  double origUV[2], uv[2];
  double denom;
  GridBool searchFlag, goodStep;
  int iteration;
  double constraint, parameterArea;

  *callAgain = FALSE;

  if ( !gridValidNode(grid, node) ) return NULL; 
  if ( !gridGeometryFace(grid, node) || 
       gridGeometryEdge(grid, node)) return grid;

  face = adjItem(adjFirst(gridFaceAdj(grid),node));
  if ( grid != gridFace(grid, face, nodes, &faceId)) return NULL;
  if ( grid != gridNodeUV(grid, node, faceId, origUV)) return NULL;    

  if ( grid != gridStoreFaceCostParameterDerivatives(grid, node ) )return NULL;

  minCost =2.1;
  minFace = EMPTY;
  for (i=0;i<gridStoredCostDegree(grid);i++){
    if (gridStoredCost(grid,i)<minCost){
      minCost = gridStoredCost(grid,i);
      minFace = i;
    }
  }

  searchFlag = FALSE;
  if (searchFlag) {
    gridStoredCostDerivative(grid, minFace, searchDirection);
  }else{
    nearestFace=EMPTY;
    nearestCost = 2.1;
    for (i=0;i<gridStoredCostDegree(grid);i++){
      if ( i != minFace){
	nearestDifference = ABS(gridStoredCost(grid,i)-minCost);
	if (nearestDifference<nearestCost) {
	  nearestFace=i;
	  nearestCost = nearestDifference;
	}
      }
    }
    if (nearestFace == EMPTY || nearestCost > 0.001 ){
      gridStoredCostDerivative(grid, minFace, searchDirection);
      gridStoredCostDerivative(grid, minFace, minDirection);
    }else{
      gridStoredCostDerivative(grid, minFace, minDirection);
      gridStoredCostDerivative(grid, nearestFace, nearestDirection);
      g00 = gridDotProduct(minDirection,minDirection);
      g11 = gridDotProduct(nearestDirection,nearestDirection);
      g01 = gridDotProduct(minDirection,nearestDirection);
      /*
       * Note: If two incedent cells have the same Cost (more specifically
       *       CostDerivative), then nearestDirection == minDirection
       *       which will result in 0/0 for nearestRatio (g00 == g11).
       *       Could have check nearestDifference != 0.0 in above loop
       *       before setting nearestFace.  Would then have same result
       *       (e.g. searchDirection == minDirection) since nearestFace
       *       would be EMPTY and previous block would execute.
       */
      denom = g00 + g11 - 2*g01;
      if( abs(denom) < 1.0e-12 ) {
        nearestRatio = 0.0;
      } else {
        nearestRatio = (g00-g01)/denom;
      }
      if (nearestRatio < 1.0 && nearestRatio > 0.0 ) {
	minRatio = 1.0 - nearestRatio;
	for (i=0;i<3;i++) searchDirection[i] 
			    = minRatio*minDirection[i]
			    + nearestRatio*nearestDirection[i];
	/* reset length to the projection of min cell to search dir*/
	length = sqrt(gridDotProduct(searchDirection,searchDirection));
	if (ABS(length) > 1.0e-12) {
	  for (i=0;i<3;i++) searchDirection[i] = searchDirection[i]/length;
	  projection = gridDotProduct(searchDirection,minDirection);
	  for (i=0;i<3;i++) searchDirection[i] = projection*searchDirection[i];
	}else{
	  gridStoredCostDerivative(grid, minFace, searchDirection);
	  gridStoredCostDerivative(grid, minFace, minDirection);
	}
      }else{
	gridStoredCostDerivative(grid, minFace, searchDirection);
	gridStoredCostDerivative(grid, minFace, minDirection);
      }
    }
  }

  length = sqrt(gridDotProduct(searchDirection,searchDirection));
  if (ABS(length) < 1.0e-12) return NULL;
  for (i=0;i<3;i++) searchDirection[i] = searchDirection[i]/length;

  alpha = 1.0;
  for (i=0;i<gridStoredCostDegree(grid);i++){
    if (i != minFace ) {
      gridStoredCostDerivative(grid,i,dCostdX);
      projection = gridDotProduct(searchDirection,dCostdX);
      deltaCost = gridStoredCost(grid,i) - minCost;
      if (ABS(length-projection) < 1.0e-8){
	currentAlpha=1.0; /* no intersection */
      }else{
	currentAlpha = deltaCost / ( length + projection);
      }
      if (currentAlpha > 0 && currentAlpha < alpha ) alpha = currentAlpha;
    }
  }

  /* printf( "node %5d deg %3d active %3d old %12.9f\n",
     node, gridStoredCostDegree(grid), minFace, minCost ); */
  
  goodStep = FALSE;
  actualImprovement = 0.0;
  lastImprovement = -10.0;
  constraint = 1.0;
  lastAlpha = alpha;
  iteration = 0;
  while (alpha > 0.1e-9 && !goodStep && iteration < 30 ) {
    iteration++;

    predictedImprovement = length*alpha;
  
    for (i=0;i<2;i++) uv[i] = origUV[i] + alpha*searchDirection[i];
    gridSetNodeUV(grid,node,faceId,uv[0],uv[1]);
    gridMinFaceAreaUV(grid,node,&parameterArea);
    if ( parameterArea < 1.0e-14 ) {
      alpha = alpha*0.6;      
      continue;
    }
    gridEvaluateFaceAtUV(grid, node, uv );
    gridNodeFaceMR(grid,node,&newCost);
    gridNodeAR(grid,node,&constraint);
    actualImprovement = newCost-minCost;
    /* printf(" alpha %12.5e predicted %12.9f actual %12.9f new %12.9f\n",
       alpha, predictedImprovement, actualImprovement, newCost); */

    if ( actualImprovement < lastImprovement &&
	 constraint > gridMinSurfaceSmoothCost(grid) ) {
      for (i=0;i<2;i++) uv[i] = origUV[i] + lastAlpha*searchDirection[i];
      gridSetNodeUV(grid,node,faceId,uv[0],uv[1]);
      gridMinFaceAreaUV(grid,node,&parameterArea);
      if ( parameterArea < 1.0e-14 ) {
	alpha = alpha*0.6;      
	continue;
      }
      gridEvaluateFaceAtUV(grid,node,uv);
      gridNodeFaceMR(grid,node,&newCost);
      gridNodeAR(grid,node,&constraint);
      actualImprovement = newCost-minCost;
      break;
    }
    
    if ( actualImprovement > 0.9*predictedImprovement &&
	 constraint > gridMinSurfaceSmoothCost(grid) ) {
      goodStep = TRUE;
    }else{
      lastImprovement = actualImprovement;
      lastAlpha = alpha;
      alpha =alpha*0.6;
    }
  }

  /* printf( "node %5d deg %3d active %3d old %8.5f new %8.5f\n",
     node, gridStoredCostDegree(grid), minFace, minCost, newCost ); */

  if ( actualImprovement <= 0.0 ||
       constraint < gridMinSurfaceSmoothCost(grid) ||
       parameterArea < 1.0e-14 ) {
    gridEvaluateFaceAtUV(grid,node,origUV);
    return NULL;
  }

  if ( actualImprovement > 1.0e-10 || goodStep ) *callAgain = TRUE;

  return grid;
}

Grid *gridOptimizeUVForVolume(Grid *grid, int node, double *dudv )
{
  double uvOrig[2], uv[2];
  int nodes[3];
  int face, faceId;
  double gold;
  double alpha[2], volume[2], area[2];
  int iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  gridNodeUV( grid, node, faceId, uvOrig);

  alpha[0] = 0.0;
  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  if (grid != gridEvaluateFaceAtUV(grid, node, uv ) ) return NULL;
  gridNodeVolume( grid, node, &volume[0] );
  gridMinFaceAreaUV(grid,node,&area[0]);

  alpha[1] = 1.0e-10;
  uv[0] = uvOrig[0] + alpha[1]*dudv[0];
  uv[1] = uvOrig[1] + alpha[1]*dudv[1];
  gridEvaluateFaceAtUV(grid, node, uv );
  gridNodeVolume( grid, node, &volume[1] );
  gridMinFaceAreaUV(grid,node,&area[1]);

  iter = 0;
  while ( ( volume[1] > volume[0] ) && 
	  ( area[1] > 1.0e-12 || area[1] > area[0] ) &&
	  ( iter < 200 ) ) {
    iter++;
    alpha[0] = alpha[1]; volume[0] = volume[1]; area[0] = area[1];
    alpha[1] = alpha[0] * gold;
    uv[0] = uvOrig[0] + alpha[1]*dudv[0];
    uv[1] = uvOrig[1] + alpha[1]*dudv[1];
    gridEvaluateFaceAtUV(grid, node, uv );
    gridNodeVolume( grid, node, &volume[1] );
    gridMinFaceAreaUV(grid,node,&area[1]);
  }

  uv[0] = uvOrig[0] + alpha[0]*dudv[0];
  uv[1] = uvOrig[1] + alpha[0]*dudv[1];
  gridEvaluateFaceAtUV(grid, node, uv );

  //printf("node %d alpha %e vol %f uv %e %e\n",node,alpha[0],volume[0],uv[0],uv[1]);
  
  return grid;
}

Grid *gridOptimizeXYZ(Grid *grid, int node, double *dxdydz )
{
  double xyzOrig[3];
  double xyz[3];
  double gold;
  double alpha[2], ar[2];
  int i, iter;

  gold = ( 1.0 + sqrt(5.0) ) / 2.0;

  if (grid != gridNodeXYZ(grid,node,xyzOrig)) return NULL;

  alpha[0] = 0.0;
  for(i=0;i<3;i++) xyz[i] = xyzOrig[i] + alpha[0]*dxdydz[i];
  gridSetNodeXYZ(grid,node,xyz);
  gridNodeAR( grid, node, &ar[0] );

  alpha[1] = 1.0e-10;
  for(i=0;i<3;i++) xyz[i] = xyzOrig[i] + alpha[1]*dxdydz[i];
  gridSetNodeXYZ(grid,node,xyz);
  gridNodeAR( grid, node, &ar[1] );

  iter = 0;
  while ( ar[1] > ar[0] && ar[1] > 0.0 && iter < 100){
    iter++;
    alpha[0] = alpha[1]; ar[0] = ar[1];
    alpha[1] = alpha[0] * gold;
    for(i=0;i<3;i++) xyz[i] = xyzOrig[i] + alpha[1]*dxdydz[i];
    gridSetNodeXYZ(grid,node,xyz);
    gridNodeAR( grid, node, &ar[1] );
  }

  for(i=0;i<3;i++) xyz[i] = xyzOrig[i] + alpha[0]*dxdydz[i];
  gridSetNodeXYZ(grid,node,xyz);
  
  return grid;
}

static GridBool gridThisNodeCanBeModifiedInThisPhase( Grid *grid, int node )
{
  if ( gridGeometryNode(grid, node) ) return FALSE;

  switch( gridPhase(grid) ) {
  case (gridALL_PHASE) :
    return TRUE;
    break;
  case (gridEDGE_PHASE) :
    return ( gridGeometryEdge(grid, node) );
    break;
  case (gridFACE_PHASE) :
    if ( gridGeometryEdge(grid, node) ) return FALSE;
    return ( gridGeometryFace(grid, node ) );
    break;
  case (gridVOL_PHASE) :
    return ( !gridGeometryFace(grid, node ) );
    break;
  default : 
    printf("%s: %d: gridThisEdgeCanBeModifiedInThisPhase: phase %d unknown?\n",
	   __FILE__,__LINE__,gridPhase(grid));
    return FALSE;
  }
}

Grid *gridSmooth( Grid *grid, double optimizationLimit, double laplacianLimit )
{
  int cell, nodes[4];
  int i, node, ranking;
  double ar;
  double *cost;
  Plan *plan;
  
  if ( optimizationLimit < 0.0 ) optimizationLimit = 0.40;
  if ( laplacianLimit    < 0.0 ) laplacianLimit    = 0.60;

  cost = (double *)malloc(gridMaxNode(grid)*sizeof(double));
  for (node=0;node<gridMaxNode(grid);node++) cost[node]=2.0;

  plan = planCreate( gridNNode(grid)/2, MAX(gridNNode(grid)/10,1000) );
  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid,cell,nodes)) {
      ar = gridAR(grid, nodes);
      for(i=0;i<4;i++) {
	cost[nodes[i]] = MIN(cost[nodes[i]],ar);
      }
    }
  }
  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode(grid,node) && !gridNodeFrozen( grid, node ) &&
	 gridThisNodeCanBeModifiedInThisPhase( grid, node ) ) {
      if ( cost[node] < laplacianLimit ) {
	planAddItemWithPriority( plan, node, 1.0 - cost[node] );
      }
    }
  }
  free(cost);
  planDeriveRankingsFromPriorities(plan);
  for ( ranking=planSize(plan)-1; ranking>=0; ranking-- ) { 
    node = planItemWithThisRanking(plan,ranking);
    gridNodeAR(grid,node,&ar);
    if (ar < optimizationLimit) {
      gridSmoothNode( grid, node, TRUE );
    }else{
      if (ar < laplacianLimit && !gridGeometryFace( grid, node )) {
	gridSmartLaplacian( grid, node ); 
      }
    }
  }
  planFree(plan);
  return grid;
}

Grid *gridSmoothFaceMR( Grid *grid, double optimizationLimit )
{
  int node;
  double mr;

  if ( optimizationLimit < 0.0 ) optimizationLimit = 0.40;

  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) && gridGeometryFace( grid, node )) {
      gridNodeFaceMR(grid,node,&mr);
      if (mr < optimizationLimit) {
	gridSmoothNodeFaceMR( grid, node );
	gridSmoothNodeFaceMR( grid, node );
      }
    }
  }
  return grid;
}

Grid *gridSmoothVolume( Grid *grid )
{
  int node;
  double ar, optimizationLimit, laplacianLimit;
  optimizationLimit =0.30;
  laplacianLimit =0.60;
  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) && 
	 !gridGeometryFace( grid, node ) &&
	 gridNodeLocal(grid,node) ) {
      gridNodeAR(grid,node,&ar);
      if (ar < laplacianLimit && !gridGeometryFace( grid, node )) {
	gridSmartLaplacian( grid, node ); 
	gridNodeAR(grid,node,&ar);
      }
      if (ar < optimizationLimit) {
	gridSmoothNode( grid, node, TRUE );
      }
    }
  }
  return grid;
}

Grid *gridSmartLaplacian(Grid *grid, int node )
{
  double origAR, newAR;
  double origXYZ[3], xyz[3], nodeXYZ[3];
  double oneOverNCell;
  AdjIterator it;
  int nodes[4], ncell, inode, ixyz;
  
  gridNodeAR(grid, node, &origAR);
  if ( NULL == gridNodeXYZ(grid, node, origXYZ)) return NULL;

  xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
  ncell =0;

  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it) ; 
	it = adjNext(it) ){
    ncell++;
    gridCell(grid,adjItem(it),nodes);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      gridNodeXYZ(grid,nodes[inode],nodeXYZ);
      for (ixyz = 0 ; ixyz < 3 ; ixyz++ ) xyz[ixyz] += nodeXYZ[ixyz];
    }
  }
  oneOverNCell = 1.0/(double)(ncell*3);
  for (ixyz = 0 ; ixyz < 3 ; ixyz++ ){  
    xyz[ixyz] -= origXYZ[ixyz] * (double)ncell ;
    xyz[ixyz] = xyz[ixyz] * oneOverNCell;
  }
  gridSetNodeXYZ(grid,node,xyz);
  gridNodeAR(grid, node, &newAR);
  
  if ( origAR > newAR ) {
    gridSetNodeXYZ(grid,node,origXYZ);
    return NULL;
  }

  return grid;
}

Grid *gridSmartVolumeLaplacian(Grid *grid, int node )
{
  double origVol, newVol;
  double origXYZ[3];
  
  if ( NULL == gridNodeXYZ(grid, node, origXYZ)) return NULL;
  if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
    gridNodeMinCellJacDet2(grid,node,&origVol);
  } else {
    gridNodeVolume(grid,node,&origVol);
  }

  if (grid != gridNodeLaplacian( grid, node )) {
    gridSetNodeXYZ(grid,node,origXYZ);
    return NULL;
  }

  if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
    gridNodeMinCellJacDet2(grid,node,&newVol);
  } else {
    gridNodeVolume(grid,node,&newVol);
  }
  if ( origVol > newVol ) {
    gridSetNodeXYZ(grid,node,origXYZ);
    return NULL;
  }

  return grid;
}

Grid *gridNodeLaplacian(Grid *grid, int node )
{
  double origXYZ[3], xyz[3], nodeXYZ[3];
  double oneOverNCell;
  AdjIterator it;
  int nodes[4], ncell, inode, ixyz;
  
  if ( NULL == gridNodeXYZ(grid, node, origXYZ)) return NULL;

  xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
  ncell =0;

  for ( it = adjFirst(gridCellAdj(grid),node); 
	adjValid(it) ; 
	it = adjNext(it) ){
    ncell++;
    gridCell(grid,adjItem(it),nodes);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      gridNodeXYZ(grid,nodes[inode],nodeXYZ);
      for (ixyz = 0 ; ixyz < 3 ; ixyz++ ) xyz[ixyz] += nodeXYZ[ixyz];
    }
  }
  oneOverNCell = 1.0/(double)(ncell*3);
  for (ixyz = 0 ; ixyz < 3 ; ixyz++ ){  
    xyz[ixyz] -= origXYZ[ixyz] * (double)ncell ;
    xyz[ixyz] = xyz[ixyz] * oneOverNCell;
  }
  gridSetNodeXYZ(grid,node,xyz);

  return grid;
}

Grid *gridSmartAreaUVLaplacian(Grid *grid, int node )
{

  double origUV[2], avgUV[2], nodeUV[2];

  double origArea, newArea;

  double oneOverNFace2;
  AdjIterator it;
  int face, nodes[3], faceId;
  int nface, inode;
  
  if (!gridGeometryFace(grid,node)) return NULL;
  if (gridGeometryBetweenFace(grid,node)) return grid;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  if ( NULL == gridNodeUV(grid, node, faceId, origUV)) return NULL;

  gridMinFaceAreaUV(grid,node,&origArea);

  avgUV[0] = 0.0; avgUV[1] = 0.0;
  nface =0;

  for ( it = adjFirst(gridFaceAdj(grid),node); 
	adjValid(it) ; 
	it = adjNext(it) ){
    nface++;
    gridFace(grid,adjItem(it),nodes,&faceId);
    for ( inode = 0 ; inode < 3 ; inode++ ){
      if (nodes[inode] != node) { /* skip central node */
	gridNodeUV(grid,nodes[inode],faceId,nodeUV);
	avgUV[0] += nodeUV[0];
	avgUV[1] += nodeUV[1];
      }
    }
  }

  /* each surrounding node is added to avgUV twice, so divide by 2*nface */
  oneOverNFace2 = 1.0/((double)(nface*2));
  avgUV[0] *= oneOverNFace2;
  avgUV[1] *= oneOverNFace2;

  gridSetNodeUV(grid,node,faceId,avgUV[0],avgUV[1]);
  gridMinFaceAreaUV(grid,node,&newArea);
  
  if ( origArea > newArea ) {
    gridSetNodeUV(grid,node,faceId,origUV[0],origUV[1]);
    return NULL;
  }

  return grid;
}

Grid *gridStoreVolumeCostDerivatives (Grid *grid, int node )
{
  AdjIterator it;
  int nodes[4], orientedNodes[4];
  double AR, dARdX[3];

  if ( !gridValidNode( grid, node) ) return NULL;

  gridClearStoredCost( grid );
  for ( it = adjFirst(gridCellAdj(grid),node); adjValid(it); it = adjNext(it) ){
    gridCell(grid,adjItem(it),nodes);
    orientedNodes[0] = node;
    if (node == nodes[0]){
      orientedNodes[1] = nodes[1];
    }else{
      orientedNodes[1] = nodes[0];
    }
    gridOrient( grid, nodes, orientedNodes);
    if (grid != gridCellARDerivative(grid, orientedNodes, &AR, dARdX ) ) {
      gridClearStoredCost( grid );
      return NULL;
    }
    if (grid != gridStoreCost(grid, AR, dARdX ) ) {
      gridClearStoredCost( grid );
      return NULL;
    }
  }

  return grid;
}

Grid *gridStoreFaceCostParameterDerivatives (Grid *grid, int node )
{
  AdjIterator it;
  int face, faceId, nodes[3];
  int swapnode;
  double cost, dMRdx[3], costDerivative[3];
  double uv[2], xyzProj[3], du[3], dv[3];
  int vol=1;

  if ( !gridValidNode( grid, node) ) return NULL;

  gridClearStoredCost( grid );
  for ( it = adjFirst(gridFaceAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ){
    face = adjItem(it);
    if ( grid != gridFace(grid,face,nodes,&faceId) ) {
      gridClearStoredCost( grid );
      return NULL;
    }
    /* orient face so that nodes[0] is node for differentiation */
    if (node == nodes[1]) {
      swapnode = nodes[0];
      nodes[0] = nodes[1];
      nodes[1] = nodes[2];
      nodes[2] = swapnode;
    }
    if (node == nodes[2]) {
      swapnode = nodes[2];
      nodes[2] = nodes[1];
      nodes[1] = nodes[0];
      nodes[0] = swapnode;
    }
    if (grid != gridFaceMRDerivative(grid, nodes, &cost, dMRdx ) ) {
      gridClearStoredCost( grid );
      return NULL;
    }
    gridNodeUV( grid, node, faceId, uv);
    if ( !CADGeom_PointOnFace( vol, faceId,   
			       uv, xyzProj, 1, du, dv, NULL, NULL, NULL) ) {
      printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
      gridClearStoredCost( grid );
      return NULL;      
    }
    costDerivative[0] = dMRdx[0]*du[0] + dMRdx[1]*du[1] + dMRdx[2]*du[2] ; 
    costDerivative[1] = dMRdx[0]*dv[0] + dMRdx[1]*dv[1] + dMRdx[2]*dv[2] ; 
    costDerivative[2] = 0.0;
    if (grid != gridStoreCost(grid, cost, costDerivative ) ) {
      gridClearStoredCost( grid );
      return NULL;
    }
  }
  return grid;
}

Grid *gridLinearProgramXYZ(Grid *grid, int node, GridBool *callAgain )
{
  int i, minCell, nearestCell;
  double minAR, nearestAR, nearestDifference, newAR, searchDirection[3];
  double g00, g01, g11, minRatio, nearestRatio;
  double length, projection;
  double deltaAR, currentAlpha, alpha, lastAlpha;
  double predictedImprovement, actualImprovement, lastImprovement;
  double minDirection[3], nearestDirection[3], dARdX[3];
  double origXYZ[3], xyz[3];
  double denom;
  GridBool searchFlag, goodStep;
  int iteration;

  *callAgain = FALSE;

  if ( grid != gridNodeXYZ(grid, node, origXYZ)) return NULL;
  if ( grid != gridStoreVolumeCostDerivatives(grid, node ) ) return NULL;

  minAR =2.1;
  minCell = EMPTY;
  for (i=0;i<gridStoredCostDegree(grid);i++){
    if (gridStoredCost(grid,i)<minAR){
      minAR = gridStoredCost(grid,i);
      minCell = i;
    }
  }

  searchFlag = FALSE;
  if (searchFlag) {
    gridStoredCostDerivative(grid, minCell, searchDirection);
  }else{
    nearestCell=EMPTY;
    nearestAR = 2.1;
    for (i=0;i<gridStoredCostDegree(grid);i++){
      if ( i != minCell){
	nearestDifference = ABS(gridStoredCost(grid,i)-minAR);
	if (nearestDifference<nearestAR) {
	  nearestCell=i;
	  nearestAR = nearestDifference;
	}
      }
    }
    if (nearestCell == EMPTY || nearestAR > 0.001 ){
      gridStoredCostDerivative(grid, minCell, searchDirection);
      gridStoredCostDerivative(grid, minCell, minDirection);
    }else{
      gridStoredCostDerivative(grid, minCell, minDirection);
      gridStoredCostDerivative(grid, nearestCell, nearestDirection);
      g00 = gridDotProduct(minDirection,minDirection);
      g11 = gridDotProduct(nearestDirection,nearestDirection);
      g01 = gridDotProduct(minDirection,nearestDirection);
      /*
       * Note: If two incedent cells have the same AR (more specifically
       *       ARDerivative), then nearestDirection == minDirection
       *       which will result in 0/0 for nearestRatio (g00 == g11).
       *       Could have check nearestDifference != 0.0 in above loop
       *       before setting nearestCell.  Would then have same result
       *       (e.g. searchDirection == minDirection) since nearestCell
       *       would be EMPTY and previous block would execute.
       */
      denom = g00 + g11 - 2*g01;
      if( ABS(denom) < 1.0e-12 ) {
        nearestRatio = 0.0;
      } else {
        nearestRatio = (g00-g01)/denom;
      }
      if (nearestRatio < 1.0 && nearestRatio > 0.0 ) {
	minRatio = 1.0 - nearestRatio;
	for (i=0;i<3;i++) searchDirection[i] 
			    = minRatio*minDirection[i]
			    + nearestRatio*nearestDirection[i];
	/* reset length to the projection of min cell to search dir*/
	length = sqrt(gridDotProduct(searchDirection,searchDirection));
	for (i=0;i<3;i++) searchDirection[i] = searchDirection[i]/length;
	projection = gridDotProduct(searchDirection,minDirection);
	for (i=0;i<3;i++) searchDirection[i] = projection*searchDirection[i];
	//printf("node %5d min %10.7f near %10.7f\n",node,minRatio,nearestRatio);
      }else{
	gridStoredCostDerivative(grid, minCell, searchDirection);
	gridStoredCostDerivative(grid, minCell, minDirection);
      }
    }
  }

  length = sqrt(gridDotProduct(searchDirection,searchDirection));
  if (ABS(length) < 1.0e-12) return NULL;
  for (i=0;i<3;i++) searchDirection[i] = searchDirection[i]/length;

  alpha = 1.0;
  for (i=0;i<gridStoredCostDegree(grid);i++){
    if (i != minCell ) {
      gridStoredCostDerivative(grid,i,dARdX);
      projection = gridDotProduct(searchDirection,dARdX);
      deltaAR = gridStoredCost(grid,i) - minAR;
      if (ABS(length-projection) < 1.0e-12){
	currentAlpha=1.0; /* no intersection */
      }else{
	currentAlpha = deltaAR / ( length + projection);
      }
      if (currentAlpha > 0 && currentAlpha < alpha ) alpha = currentAlpha;
    }
  }

  //printf( "node %5d deg %3d active %3d old %12.9f\n",
  //	  node, gridStoredCostDegree(grid), minCell, minAR );
  
  goodStep = FALSE;
  actualImprovement = 0.0;
  lastImprovement = -10.0;
  lastAlpha = alpha;
  iteration = 0;
  while (alpha > 1.0e-12 && !goodStep && iteration < 30 ) {
    iteration++;

    predictedImprovement = length*alpha;
  
    for (i=0;i<3;i++) xyz[i] = origXYZ[i] + alpha*searchDirection[i];
    gridSetNodeXYZ(grid,node,xyz);
    gridNodeAR(grid,node,&newAR);
    actualImprovement = newAR-minAR;
    //printf(" alpha %12.5e predicted %12.9f actual %12.9f new %12.9f\n",
    //	   alpha, predictedImprovement, actualImprovement, newAR);

    if ( actualImprovement > 0.0 && actualImprovement < lastImprovement) {
      for (i=0;i<3;i++) xyz[i] = origXYZ[i] + lastAlpha*searchDirection[i];
      gridSetNodeXYZ(grid,node,xyz);
      gridNodeAR(grid,node,&newAR);
      actualImprovement = newAR-minAR;
      goodStep = TRUE;
    }
    
    if ( actualImprovement > 0.9*predictedImprovement  ){
      goodStep = TRUE;
    }else{
      lastImprovement = actualImprovement;
      lastAlpha = alpha;
      alpha =alpha*0.5;
    }
  }

  //printf( "node %5d deg %3d active %3d old %8.5f new %8.5f\n",
  //  node, gridStoredCostDegree(grid), minCell, minAR, newAR );

  if ( actualImprovement <= 0.0  ){
    gridSetNodeXYZ(grid,node,origXYZ);
    return NULL;
  }

  if ( (actualImprovement > 1.0e-12 || goodStep) && newAR < 0.999) 
    *callAgain = TRUE;

  return grid;
}

Grid *gridSmoothInvalidCellNodes(Grid *grid)
{
  int cell, nodes[4], i, node;

  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, nodes)) {
      if ( -0.5 > gridAR(grid,nodes) ) {
	for (i=0;i<4;i++) {
	  node = nodes[i];
	  if (!gridGeometryFace(grid,node)) {
	    gridSmartVolumeLaplacian( grid, node );	    
	    gridSmoothNodeVolumeSimplex(grid, node);
	  }
	}
      }
    }
  }
  return grid;
}

Grid *gridSmoothNodeVolume( Grid *grid, int node )
{
  if ( !gridValidNode(grid, node)   ||
       gridNodeFrozen(grid, node)   ||
       gridGeometryBetweenFace(grid, node) ||
       gridNodeGhost(grid, node)    ) return NULL;
  if (gridGeometryFace(grid, node)) {
    gridSmoothNodeVolumeWithSurf( grid, node );
  }else{
    gridSmartVolumeLaplacian( grid, node );
    gridSmoothNodeVolumeSimplex( grid, node );
  }
  return grid;
}

Grid *gridSmoothNodeVolumeWithSurf( Grid *grid, int node )
{
  int face, nodes[3], faceId;
  double volume, dVoldx[3];
  double uv[2], xyzProj[3], du[3], dv[3];
  double dVoldu[2];
  int vol=1;
  
  if ( !gridValidNode(grid, node)   ||
       gridNodeFrozen(grid, node)   ||
       gridGeometryEdge(grid, node) ||
       !gridGeometryFace(grid, node) ||
       gridNodeGhost(grid, node)    ) return NULL;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  gridFace(grid,face,nodes,&faceId);
  gridNodeVolumeDerivative ( grid, node, &volume, dVoldx);
  gridNodeUV( grid, node, faceId, uv);
  if ( !CADGeom_PointOnFace( vol, faceId,   
			     uv, xyzProj, 1, du, dv, NULL, NULL, NULL) )
    printf ( "ERROR: CADGeom_PointOnFace, %d: %s\n",__LINE__,__FILE__ );
  
  dVoldu[0] = dVoldx[0]*du[0] + dVoldx[1]*du[1] + dVoldx[2]*du[2] ; 
  dVoldu[1] = dVoldx[0]*dv[0] + dVoldx[1]*dv[1] + dVoldx[2]*dv[2] ; 
  if (grid != gridOptimizeUVForVolume( grid, node, dVoldu ) ) return NULL;

  return grid;
}

static double reflect( Grid *grid,
		       double simplex[4][3], double volume[4], double avgXYZ[3],
		       int node, int worst, double factor)
{
  int i;
  double factor1, factor2;
  double reflectedXYZ[3];
  double reflectedVolume;

  factor1 = (1.0-factor) / 3.0;
  factor2 = factor1 - factor;

  for(i=0;i<3;i++) 
    reflectedXYZ[i] = factor1*avgXYZ[i] - factor2*simplex[worst][i];

  gridSetNodeXYZ(grid,node,reflectedXYZ );
  if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
    gridNodeMinCellJacDet2(grid,node,&reflectedVolume);
  } else {
    gridNodeVolume(grid,node,&reflectedVolume);
  }

  if ( reflectedVolume > volume[worst] ) {
    volume[worst] = reflectedVolume;
    for(i=0;i<3;i++) avgXYZ[i] += ( reflectedXYZ[i] - simplex[worst][i] );
    for(i=0;i<3;i++) simplex[worst][i] = reflectedXYZ[i];
  }

  return reflectedVolume;
}

static Grid *gridMakeFacesFromSimplex(Grid *grid, 
				      double simplex[4][3], int faceId)
{
  int node;
  int nodes[4];
  for(node=0;node<4;node++)
    nodes[node]=gridAddNode(grid,
			    simplex[node][0],
			    simplex[node][1],
			    simplex[node][2]);
  gridAddFace(grid,nodes[0],nodes[1],nodes[2],faceId);
  gridAddFace(grid,nodes[0],nodes[1],nodes[3],faceId);
  gridAddFace(grid,nodes[1],nodes[2],nodes[3],faceId);
  gridAddFace(grid,nodes[0],nodes[3],nodes[2],faceId);
  return grid;
}

Grid *gridSmoothNodeVolumeSimplex( Grid *grid, int node )
{
  int evaluations;
  int s, i;
  double origXYZ[3], avgXYZ[3];
  double simplex[4][3];
  double volume[4];
  double lengthScale;
  int best, worst, secondworst; 
  double newVolume, savedVolume;
  GridBool makefaces = FALSE;
  int faceId = 1;

  if ( NULL == gridNodeXYZ(grid, node, origXYZ)) return NULL;

  lengthScale = 0.1*gridAverageEdgeLength(grid, node );

  for(s=0;s<4;s++)
    for(i=0;i<3;i++)
      simplex[s][i] = origXYZ[i];

  simplex[1][0] += lengthScale;
  simplex[2][1] += lengthScale;
  simplex[3][2] += lengthScale;

  for(s=0;s<4;s++) {
    gridSetNodeXYZ(grid, node, simplex[s]);
    if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
      gridNodeMinCellJacDet2(grid,node,&volume[s]);
    } else {
      gridNodeVolume(grid,node,&volume[s]);
    }
  }

  for(i=0;i<3;i++) avgXYZ[i] = 0.0;
  for(s=0;s<4;s++)
    for(i=0;i<3;i++) avgXYZ[i] += simplex[s][i];

  evaluations = 4;
  while (evaluations < 1000 ) {

    best = 0;
    if ( volume[0] > volume[1] ) {
      secondworst = 0;
      worst = 1;
    }else{
      secondworst = 1;
      worst = 0;
    }
    
    for(s=0;s<4;s++) {
      if (volume[s]>=volume[best]) best = s;
      if (volume[s]<volume[worst]) {
	secondworst = worst;
	worst = s;
      }else{
	if ( s!=worst && volume[s]<volume[secondworst]) secondworst = s;
      }
    }

    /* printf( "evaluations%6d best%20.15f worst%20.15f\n", 
               evaluations, volume[best], volume[worst]); */
    if (makefaces) gridMakeFacesFromSimplex(grid, simplex, ++faceId);

    if (volume[best]-volume[worst] < ABS(1.0e-10*volume[best])) break;

    evaluations++;
    newVolume = reflect( grid, simplex, volume, avgXYZ, node, worst, -1.0 );
    if ( newVolume >= volume[best] ) {
      evaluations++;
      newVolume = reflect( grid, simplex, volume, avgXYZ, node, worst, 2.0 );
    } else {
      if (newVolume <= volume[secondworst]) {
	savedVolume = volume[worst];
	evaluations++;
	newVolume = reflect( grid, simplex, volume, avgXYZ, node, worst, 0.5 );
	if (newVolume <= savedVolume) {
	  for(s=0;s<4;s++) {
	    if (s != best) {
	      for(i=0;i<3;i++) 
		simplex[s][i]=0.5*(simplex[s][i]+simplex[best][i]);
	      gridSetNodeXYZ(grid, node, simplex[s]);
	      if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
		gridNodeMinCellJacDet2(grid,node,&volume[s]);
	      } else {
		gridNodeVolume(grid,node,&volume[s]);
	      }
	    }
	  }
	}      
      }
    }
  }    

  best = 0;
  for(s=1;s<4;s++) if (volume[s]>=volume[best]) best = s;

  gridSetNodeXYZ(grid, node, simplex[best]);

  return grid;
}

Grid *gridRelaxNegativeCells(Grid *grid, GridBool dumpTecplot )
{
  int cell, nodes[4], i, node;
  double volume;
  char filename[256];

  if (dumpTecplot) {
    sprintf(filename,"gridNegativeCell%04d.t",gridPartId(grid));
    gridWriteTecplotSurfaceGeom(grid, filename);
  }

  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, nodes)) {
      volume = gridVolume(grid,nodes);
      if (0.0>=volume){
	if (dumpTecplot) {
	  double costs[4];
	  costs[0] = costs[1] = costs[2] = costs[3] = volume;
	  gridWriteTecplotCellGeom(grid,nodes,costs,filename);
	}
	for (i=0;i<4;i++) {
	  node = nodes[i];
	  gridSmoothVolumeNearNode(grid, node, FALSE);
	}
      }
    }
  }
  return grid;
}

Grid *gridSmoothVolumeNearNode(Grid *grid, int node, GridBool smoothOnSurface )
{
  int i, i0, nodes[4], nodes0[4];
  int nlist, smooth, look, nodelist[SMOOTHDEG];
  GridBool looking;
  AdjIterator it, it0;

  if (!gridValidNode(grid,node)) return NULL;

  nlist =0;
  for ( it0 = adjFirst(gridCellAdj(grid),node); 
	adjValid(it0); 
	it0 = adjNext(it0) ){
    gridCell(grid, adjItem(it0), nodes0);
    for (i0=0;i0<4;i0++) {
      for ( it = adjFirst(gridCellAdj(grid),nodes0[i0]); 
	    adjValid(it); 
	    it = adjNext(it) ){
	gridCell(grid, adjItem(it), nodes);
	for (i=0;i<4;i++) {
	  if ( !smoothOnSurface && gridGeometryFace(grid, nodes[i]) ) continue;
	  looking = (nlist<=SMOOTHDEG);
	  look = 0;
	  for (look=0;look<nlist && looking ; look++){
	    looking = (nodelist[look] != nodes[i]);
	  }
	  if (looking && nlist<=SMOOTHDEG){
	    nodelist[nlist] = nodes[i];
	    nlist++;
	  }
	}
      }
    }
  }      

  for (smooth=0;smooth<10;smooth++)
    for (i=0;i<nlist;i++) 
      gridSmoothNodeVolume( grid, nodelist[i]);

  return grid;
}

GridBool nearestOnEdge(int vol, int edgeId, double *xyz, double *t,
                       double *xyznew)
{
  double ptLocal[3], pt[3];

  /* Local coordinate system */
  if( !CADGeom_DisplacementIsIdentity(vol) ) {
    if( !CADGeom_UnMapPoint(vol,xyz,ptLocal) ) {
      printf("%s: %d: CADGeom_UnMapPoint failed.\n",__FILE__,__LINE__);
      return FALSE;
    }
  } else {
    ptLocal[0] = xyz[0]; ptLocal[1] = xyz[1]; ptLocal[2] = xyz[2];
  }

  /* Project Point */
  if (!CADGeom_NearestOnEdge( vol, edgeId, ptLocal, t, pt) ) {
    printf("%s: %d: CADGeom_NearestOnEdge failed.\n",__FILE__,__LINE__);
    return FALSE;  
  }

  /* Global coordinate system */
  if( !CADGeom_DisplacementIsIdentity(vol) ) {
    if( !CADGeom_MapPoint(vol,pt,xyznew) ) {
      printf("%s: %d: CADGeom_MapPoint failed.\n",__FILE__,__LINE__);
      return FALSE;  
    }
  } else {
    xyznew[0] = pt[0]; xyznew[1] = pt[1]; xyznew[2] = pt[2];
  }

  return TRUE;
}

GridBool nearestOnFace(int vol, int faceId, double *xyz, double *uv,
                       double *xyznew)
{
  double ptLocal[3], pt[3];

  /* Local coordinate system */
  if( !CADGeom_DisplacementIsIdentity(vol) ) {
    if( !CADGeom_UnMapPoint(vol,xyz,ptLocal) ) {
      printf("%s: %d: CADGeom_UnMapPoint failed.\n",__FILE__,__LINE__);
      return FALSE;
    }
  } else {
    ptLocal[0] = xyz[0]; ptLocal[1] = xyz[1]; ptLocal[2] = xyz[2];
  }

  /* Project Point */
  if (!CADGeom_NearestOnFace( vol, faceId, ptLocal, uv, pt) ) {
    printf("%s: %d: CADGeom_NearestOnFace failed.\n",__FILE__,__LINE__);
    return FALSE;  
  }

  /* Global coordinate system */
  if( !CADGeom_DisplacementIsIdentity(vol) ) {
    if( !CADGeom_MapPoint(vol,pt,xyznew) ) {
      printf("%s: %d: CADGeom_MapPoint failed.\n",__FILE__,__LINE__);
      return FALSE;  
    }
  } else {
    xyznew[0] = pt[0]; xyznew[1] = pt[1]; xyznew[2] = pt[2];
  }

  return TRUE;
}

Grid *gridUntangleBadFaceParameters(Grid *grid)
{
  int face, nodes[3], faceId;
  int node, *hits;
  
  hits = (int *)malloc( gridMaxNode(grid) * sizeof(int) );
  for (node=0;node<gridMaxNode(grid); node++) hits[node]=0;

  for (face=0;face<gridMaxFace(grid);face++) {
    if (grid == gridFace(grid,face,nodes,&faceId) ) {
      if ( 1.0e-14 > gridFaceAreaUV(grid, face) ) {
	hits[nodes[0]]++; hits[nodes[1]]++; hits[nodes[2]]++;
      }
    }
  }

  for (node=0;node<gridMaxNode(grid); node++) {
    if (hits[node] > 2) {
      printf("untangling node%10d with hits%3d\n",node,hits[node]);
      gridSmoothNodeFaceAreaUV(grid, node );
    }
  }  
  return grid;
}

Grid *gridSmoothNodeFaceAreaUV(Grid *grid, int node )
{
  if (!gridGeometryFace(grid,node)) return NULL;
  if (gridGeometryBetweenFace(grid,node)) return grid;
  gridSmartAreaUVLaplacian(grid, node);
  if (grid != gridSmoothNodeFaceAreaUVSimplex(grid, node )) return NULL;
  if (grid != gridSmoothNodeFaceAreaUVSimplex(grid, node )) return NULL;
  if (grid != gridSmoothNodeFaceAreaUVSimplex(grid, node )) return NULL;
  return grid;
}

static double reflectFaceAreaUV( Grid *grid,
		       double simplex[3][2], double area[3], double avgUV[2],
		       int node, int faceId, int worst, double factor)
{
  int i;
  double factor1, factor2;
  double reflectedUV[2];
  double reflectedArea;

  factor1 = (1.0-factor) / 2.0;
  factor2 = factor1 - factor;

  for(i=0;i<2;i++) 
    reflectedUV[i] = factor1*avgUV[i] - factor2*simplex[worst][i];

  gridSetNodeUV(grid, node, faceId, reflectedUV[0],  reflectedUV[1] );
  gridMinFaceAreaUV(grid,node,&reflectedArea);

  if ( reflectedArea > area[worst] ) {
    area[worst] = reflectedArea;
    for(i=0;i<2;i++) avgUV[i] += ( reflectedUV[i] - simplex[worst][i] );
    for(i=0;i<2;i++) simplex[worst][i] = reflectedUV[i];
  }

  return reflectedArea;
}

Grid *gridSmoothNodeFaceAreaUVSimplex( Grid *grid, int node )
{
  int evaluations;
  int s, i;
  double origUV[2], avgUV[2];
  double simplex[3][2];
  double area[3];
  double lengthScale;
  int best, worst, middle; 
  double newArea, savedArea;
  int face, faceId, nodes[3];
  /* GridBool makefaces = FALSE; int debugFaceId = 1; */

  if (!gridGeometryFace(grid,node)) return NULL;
  if (gridGeometryBetweenFace(grid,node)) return NULL;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  if ( NULL == gridNodeUV(grid, node, faceId, origUV)) return NULL;

  lengthScale = 0.1*gridAverageEdgeLength(grid, node );/* FAILED??? */

  for(s=0;s<3;s++)
    for(i=0;i<2;i++)
      simplex[s][i] = origUV[i];

  simplex[1][0] += lengthScale;
  simplex[2][1] += lengthScale;

  for(s=0;s<3;s++) {
    gridSetNodeUV(grid, node, faceId, simplex[s][0], simplex[s][1]);
    gridMinFaceAreaUV(grid,node,&area[s]);
  }

  for(i=0;i<2;i++) avgUV[i] = 0.0;
  for(s=0;s<3;s++)
    for(i=0;i<2;i++) avgUV[i] += simplex[s][i];

  evaluations = 3;
  while (evaluations < 1000 ) {

    best = 0;
    if ( area[0] > area[1] ) {
      middle = 0;
      worst = 1;
    }else{
      middle = 1;
      worst = 0;
    }
    
    for(s=0;s<3;s++) {
      if (area[s]>=area[best]) best = s;
      if (area[s]<area[worst]) {
	middle = worst;
	worst = s;
      }else{
	if ( s!=worst && area[s]<area[middle]) middle = s;
      }
    }

    /* printf( "evaluations%6d best%20.15f mid%20.15f  worst%20.15f\n", 
       evaluations, area[best], area[middle], area[worst]);*/
    /* if (makefaces) gridMakeFacesFromSimplex(grid, simplex, ++faceId); */

    if (area[best]-area[worst] < ABS(1.0e-8*area[best])) break;

    evaluations++;
    newArea = reflectFaceAreaUV( grid, simplex, area, avgUV, node,
				 faceId, worst, -1.0 );
    if ( newArea >= area[best] ) {
      evaluations++;
      newArea = reflectFaceAreaUV( grid, simplex, area, avgUV, node,
				   faceId, worst, 2.0 );
    } else {
      if (newArea <= area[middle]) {
	savedArea = area[worst];
	evaluations++;
	newArea = reflectFaceAreaUV( grid, simplex, area, avgUV, node,
				     faceId, worst, 0.5 );
	if (newArea <= savedArea) {
	  for(s=0;s<3;s++) {
	    if (s != best) {
	      for(i=0;i<2;i++) 
		simplex[s][i]=0.5*(simplex[s][i]+simplex[best][i]);
	      gridSetNodeUV(grid, node, faceId, simplex[s][0], simplex[s][1]);
	      gridMinFaceAreaUV(grid,node,&area[s]);
	    }
	  }
	}      
      }
    }
  }    

  best = 0;
  for(s=1;s<3;s++) if (area[s]>=area[best]) best = s;

  gridEvaluateFaceAtUV(grid, node, simplex[best]);
  
  return grid;
}

static double reflectVolumeUV( Grid *grid,
		       double simplex[3][2], double area[3], double avgUV[2],
		       int node, int faceId, int worst, double factor)
{
  int i;
  double factor1, factor2;
  double reflectedUV[2];
  double reflectedArea;

  factor1 = (1.0-factor) / 2.0;
  factor2 = factor1 - factor;

  for(i=0;i<2;i++) 
    reflectedUV[i] = factor1*avgUV[i] - factor2*simplex[worst][i];

  gridEvaluateFaceAtUV(grid, node, reflectedUV );
  if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
    gridNodeMinCellJacDet2(grid,node,&reflectedArea);
  } else {
    gridNodeVolume(grid,node,&reflectedArea);
  }

  if ( reflectedArea > area[worst] ) {
    area[worst] = reflectedArea;
    for(i=0;i<2;i++) avgUV[i] += ( reflectedUV[i] - simplex[worst][i] );
    for(i=0;i<2;i++) simplex[worst][i] = reflectedUV[i];
  }

  return reflectedArea;
}

Grid *gridSmoothNodeVolumeUVSimplex( Grid *grid, int node )
{
  int evaluations;
  int s, i;
  double origUV[2], avgUV[2];
  double simplex[3][2];
  double area[3];
  double lengthScale;
  int best, worst, middle; 
  double newArea, savedArea;
  int face, faceId, nodes[3];
  /* GridBool makefaces = FALSE; int debugFaceId = 1; */

  if (!gridGeometryFace(grid,node)) return NULL;
  if (gridGeometryBetweenFace(grid,node)) return NULL;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  if ( NULL == gridNodeUV(grid, node, faceId, origUV)) return NULL;

  lengthScale = 0.1*gridAverageEdgeLength(grid, node );/* FAILED??? */

  for(s=0;s<3;s++)
    for(i=0;i<2;i++)
      simplex[s][i] = origUV[i];

  simplex[1][0] += lengthScale;
  simplex[2][1] += lengthScale;

  for(s=0;s<3;s++) {
    gridEvaluateFaceAtUV(grid, node, simplex[s] );
    if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
      gridNodeMinCellJacDet2(grid,node,&area[s]);
    } else {
      gridNodeVolume(grid,node,&area[s]);
    }
  }

  for(i=0;i<2;i++) avgUV[i] = 0.0;
  for(s=0;s<3;s++)
    for(i=0;i<2;i++) avgUV[i] += simplex[s][i];

  evaluations = 3;
  while (evaluations < 1000 ) {

    best = 0;
    if ( area[0] > area[1] ) {
      middle = 0;
      worst = 1;
    }else{
      middle = 1;
      worst = 0;
    }
    
    for(s=0;s<3;s++) {
      if (area[s]>=area[best]) best = s;
      if (area[s]<area[worst]) {
	middle = worst;
	worst = s;
      }else{
	if ( s!=worst && area[s]<area[middle]) middle = s;
      }
    }

    /* printf( "evaluations%6d best%20.15f mid%20.15f  worst%20.15f\n", 
       evaluations, area[best], area[middle], area[worst]);*/
    /* if (makefaces) gridMakeFacesFromSimplex(grid, simplex, ++faceId); */

    if (area[best]-area[worst] < ABS(1.0e-8*area[best])) break;

    evaluations++;
    newArea = reflectVolumeUV( grid, simplex, area, avgUV, node,
				 faceId, worst, -1.0 );
    if ( newArea >= area[best] ) {
      evaluations++;
      newArea = reflectVolumeUV( grid, simplex, area, avgUV, node,
				   faceId, worst, 2.0 );
    } else {
      if (newArea <= area[middle]) {
	savedArea = area[worst];
	evaluations++;
	newArea = reflectVolumeUV( grid, simplex, area, avgUV, node,
				     faceId, worst, 0.5 );
	if (newArea <= savedArea) {
	  for(s=0;s<3;s++) {
	    if (s != best) {
	      for(i=0;i<2;i++) 
		simplex[s][i]=0.5*(simplex[s][i]+simplex[best][i]);
	      gridEvaluateFaceAtUV(grid, node, simplex[s] );
	      if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
		gridNodeMinCellJacDet2(grid,node,&area[s]);
	      } else {
		gridNodeVolume(grid,node,&area[s]);
	      }
	    }
	  }
	}      
      }
    }
  }    

  best = 0;
  for(s=1;s<3;s++) if (area[s]>=area[best]) best = s;

  gridEvaluateFaceAtUV(grid, node, simplex[best]);
  
  return grid;
}

#define FREE_A_C_TABLEAU(a,c,f,tableau) free(a);free(c);free(f);tableauFree(tableau);

Grid *gridUntangleAreaUV( Grid *grid, int node, int recursive_depth )
{
  int face, nodes[3], faceId, temp;
  double orig_uv[2], new_uv[2];
  double original_area, new_area;
  double uv1[2], uv2[2];
  double length2;
  int degree;
  int m, n, i, j;
  double b[3]= {0.0, 0.0, 1.0};
  double *a, *c;
  int *f;
  int basis[3];
  double at[12];
  AdjIterator it;
  Tableau *tableau;
  int vol = 1;

  if (!gridGeometryFace(grid,node)) return NULL;
  if (gridGeometryBetweenFace(grid,node)) return NULL;

  face = adjItem(adjFirst(gridFaceAdj(grid), node));
  if ( grid != gridFace(grid,face,nodes,&faceId)) return NULL;

  if ( NULL == gridNodeUV(grid, node, faceId, orig_uv)) return NULL;
  gridMinFaceAreaUV(grid, node, &original_area);

  degree = adjDegree( gridFaceAdj(grid), node );
  m = 3;
  n = degree;
  a = (double *)malloc(m*n*sizeof(double));
  c = (double *)malloc(n*sizeof(double));
  f = (int *)malloc(m*n*sizeof(double));
  tableau = tableauCreate( m, n );
  j = -1;
  for ( it = adjFirst(gridFaceAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ){
    j++;
    face = adjItem(it);
    if ( grid != gridFace(grid, face, nodes, &faceId) ) {
      FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
    }
    /* orient nodes so that the central node is in position 0 */
    if (node == nodes[1]) {
      temp = nodes[0];
      nodes[0] = nodes[1];
      nodes[1] = nodes[2];
      nodes[2] = temp;
    }
    if (node == nodes[2]) {
      temp = nodes[0];
      nodes[0] = nodes[2];
      nodes[2] = nodes[1];
      nodes[1] = temp;
    }
    if ( grid != gridNodeUV(grid, nodes[1], faceId, uv1 ) ||
	 grid != gridNodeUV(grid, nodes[2], faceId, uv2 ) ) {
      FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
    }
    if ( CADGeom_ReversedSurfaceNormal(vol, faceId) ) {
      a[0+m*j]=-0.5*(uv1[1]-uv2[1]);
      a[1+m*j]=-0.5*(uv2[0]-uv1[0]);
      c[j] = 0.5*(uv1[0]*uv2[1] - uv2[0]*uv1[1]);
    }else{
      a[0+m*j]= 0.5*(uv1[1]-uv2[1]);
      a[1+m*j]= 0.5*(uv2[0]-uv1[0]);
      c[j] = 0.5*(uv2[0]*uv1[1] - uv1[0]*uv2[1]);
    }
    a[2+m*j]=1.0;

    f[0+m*j] = nodes[1];
    f[1+m*j] = nodes[2];
    f[2+m*j] = face;

    /* remove triangles that have a very small opposite side */
    length2 = a[0+m*j]*a[0+m*j] + a[1+m*j]*a[1+m*j];
    if (length2<1.0e-20) {
      n--;j--;
      /* resize tableau */
      tableauFree( tableau );
      tableau = tableauCreate( m, n );
    }
  }
  
  /* do not contiune if there are less than m usable elements */
  if (n<m) {
    FREE_A_C_TABLEAU(a,c,f,tableau); return grid;
  }

  /* solve primal linear program with tableau method */
  tableauConstraintMatrix( tableau, a );
  tableauConstraint( tableau, b );
  tableauCost( tableau, c );
  if ( tableau != tableauSolve( tableau ) ) {
    printf( "%s: %d: %s: tableauSolve NULL\n",
	    __FILE__, __LINE__, "gridUntangleAreaUV");
    tableauFree( tableau );
    FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
  }
  tableauBasis( tableau, basis );

  /* form dual linear program and invert basis */
  for (j = 0; j<m ; j++) {
    for (i = 0; i<m ; i++) {
      at[j+m*i] = a[i+m*basis[j]];
    }
  }
  for (i = 0; i<m ; i++) {
    at[i+m*m] = c[basis[i]];
  }
  if ( !gridGaussianElimination( m, m+1, at ) ) {
    printf( "%s: %d: %s: gridGaussianElimination FALSE\n",
	    __FILE__, __LINE__, "gridUntangleAreaUV");
    FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
  }
  if ( !gridGaussianBacksolve( m, m+1, at ) ) {
    printf( "%s: %d: %s: gridGaussianBacksolve FALSE\n",
	    __FILE__, __LINE__, "gridUntangleAreaUV");
    FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
  }

  new_uv[0] = at[0+m*m];
  new_uv[1] = at[1+m*m];

  gridEvaluateFaceAtUV(grid, node, new_uv);

  gridMinFaceAreaUV(grid, node, &new_area);
  if ( new_area < original_area ) {
    //gridEvaluateFaceAtUV(grid, node, orig_uv);
  }

  if ( recursive_depth > 0 ) {
    for (j = 0; j<m ; j++) {
      for (i = 0; i<m-1 ; i++) {
	gridUntangleAreaUV(grid, f[i+m*basis[j]], recursive_depth-1 );
      }
    }
    gridUntangleAreaUV(grid, node, 0 );
  }

  FREE_A_C_TABLEAU(a,c,f,tableau); return grid;
}

Grid *gridUntangleVolume( Grid *grid, int node, int recursive_depth )
{
  int cell, unsorted_nodes[4], nodes[4], temp_node;
  double orig_xyz[3], new_xyz[3];
  double original_volume, new_volume;
  double xyz1[3], xyz2[3], xyz3[3];
  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double area2;
  int degree;
  int m, n, i, j;
  double b[4]= {0.0, 0.0, 0.0, 1.0};
  double *a, *c;
  int *f;
  int basis[4];
  double at[20];
  AdjIterator it;
  Tableau *tableau;

  double d0[9], d1[9], d2[9], d3[9];

  if (gridGeometryFace(grid,node)) return NULL;
  if ( NULL == gridNodeXYZ(grid, node, orig_xyz)) return NULL;

  gridNodeVolume(grid, node, &original_volume );

  degree = gridCellDegree( grid, node );
  m = 4;
  n = degree;
  a = (double *)malloc(m*n*sizeof(double));
  c = (double *)malloc(n*sizeof(double));
  f = (int *)malloc(m*n*sizeof(double));
  tableau = tableauCreate( m, n );
  j = -1;
  for ( it = adjFirst(gridCellAdj(grid),node);
	adjValid(it);
	it = adjNext(it) ){
    j++;
    cell = adjItem(it);
    if ( grid != gridCell(grid, cell, unsorted_nodes ) ) {
      FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
    }

    /* orient nodes so that the central node is in position 0 */    
    nodes[0] = node;
    if (unsorted_nodes[0] == node) {
      nodes[1] = unsorted_nodes[1];
    }else{
      nodes[1] = unsorted_nodes[0];
    }
    gridOrient( grid, unsorted_nodes, nodes);

    /* need an outward facing normal */
    temp_node = nodes[1];
    nodes[1] = nodes[2];
    nodes[2] = temp_node;
    if ( grid != gridNodeXYZ(grid, nodes[1], xyz1 ) ||
	 grid != gridNodeXYZ(grid, nodes[2], xyz2 ) ||
	 grid != gridNodeXYZ(grid, nodes[3], xyz3 ) ) {
      FREE_A_C_TABLEAU(a,c,f,tableau);  return NULL;
    }
    x1 = xyz1[0]; y1 = xyz1[1]; z1 = xyz1[2];
    x2 = xyz2[0]; y2 = xyz2[1]; z2 = xyz2[2];
    x3 = xyz3[0]; y3 = xyz3[1]; z3 = xyz3[2];

    d0[0] = d1[0] = d2[0] = d3[0] = x1;
    d0[1] = d1[1] = d2[1] = d3[1] = x2;
    d0[2] = d1[2] = d2[2] = d3[2] = x3;

    d0[3] = d1[3] = d2[3] = d3[3] = y1;
    d0[4] = d1[4] = d2[4] = d3[4] = y2;
    d0[5] = d1[5] = d2[5] = d3[5] = y3;

    d0[6] = d1[6] = d2[6] = d3[6] = z1;
    d0[7] = d1[7] = d2[7] = d3[7] = z2;
    d0[8] = d1[8] = d2[8] = d3[8] = z3;

    d0[0] = d0[1] = d0[2] = 1.0;
    d1[3] = d1[4] = d1[5] = 1.0;
    d2[6] = d2[7] = d2[8] = 1.0;

    a[0+m*j]=-(1.0/6.0)*gridMatrixDeterminate(d0);
    a[1+m*j]=-(1.0/6.0)*gridMatrixDeterminate(d1);
    a[2+m*j]=-(1.0/6.0)*gridMatrixDeterminate(d2);
    a[3+m*j]=1.0;
    c[j] = -(1.0/6.0)*gridMatrixDeterminate(d3);

    f[0+m*j] = nodes[1];
    f[1+m*j] = nodes[2];
    f[2+m*j] = nodes[3];
    f[3+m*j] = cell;

    /* remove tets that have a very small opposite face*/
    area2 = a[0+m*j]*a[0+m*j] + a[1+m*j]*a[1+m*j] + a[2+m*j]*a[2+m*j];
    if (area2<1.0e-20) {
      n--;j--;
      /* resize tableau */
      tableauFree( tableau );
      tableau = tableauCreate( m, n );
    }

  }
  
  /* do not contiune if there are less than m usable elements */
  if (n<m) {
    printf("node %6d too small %25.15e\n",node,original_volume);
    FREE_A_C_TABLEAU(a,c,f,tableau); return grid;
  }

  /* solve primal linear program with tableau method */
  tableauConstraintMatrix( tableau, a );
  tableauConstraint( tableau, b );
  tableauCost( tableau, c );
  if ( tableau != tableauSolve( tableau ) ) {
    printf( "%s: %d: %s: tableauSolve NULL\n",
	    __FILE__, __LINE__, "gridUntangleAreaUV");
    tableauShowTransposed( tableau );
    FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
  }
  tableauBasis( tableau, basis );

  /* form dual linear program and invert basis */
  for (j = 0; j<m ; j++) {
    for (i = 0; i<m ; i++) {
      at[j+m*i] = a[i+m*basis[j]];
    }
  }
  for (i = 0; i<m ; i++) {
    at[i+m*m] = c[basis[i]];
  }

  if ( !gridGaussianElimination( m, m+1, at ) ) {
    printf( "%s: %d: %s: gridGaussianElimination FALSE\n",
	    __FILE__, __LINE__, "gridUntangleAreaUV");
    FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
  }
  if ( !gridGaussianBacksolve( m, m+1, at ) ) {
    printf( "%s: %d: %s: gridGaussianBacksolve FALSE\n",
	    __FILE__, __LINE__, "gridUntangleAreaUV");
    FREE_A_C_TABLEAU(a,c,f,tableau); return NULL;
  }

  new_xyz[0] = at[0+m*m];
  new_xyz[1] = at[1+m*m];
  new_xyz[2] = at[2+m*m];

  gridSetNodeXYZ(grid, node,  new_xyz);

  gridNodeVolume(grid, node, &new_volume );
  if ( new_volume < original_volume ) {
    /* printf("node %6d not improved %25.15e %25.15e\n",
       node,original_volume,new_volume); */
  }

  if ( recursive_depth > 0 ) {
    for (j = 0; j<m ; j++) {
      for (i = 0; i<m-1 ; i++) {
	gridUntangleVolume(grid, f[i+m*basis[j]], recursive_depth-1 );
      }
    }
    gridUntangleVolume(grid, node, 0 );
  }

  FREE_A_C_TABLEAU(a,c,f,tableau); return grid;
}

