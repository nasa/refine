
/* Compute the intersection of geometry elements
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include "intersect.h"
#include "gridmath.h"

/* TRUE if 0.0 < segment0,segment1 < 1.0 unless colinear */
GridBool intersectSegmentSegment(double *segment0_node0,double *segment0_node1,
				 double *segment1_node0,double *segment1_node1,
				 double *segment0, double *segment1 )
{
  double u[3], v[3], w[3];
  double uu, uv, vv, uw, vw;
  double denom;

  gridSubtractVector(segment0_node1,segment0_node0,u);
  gridSubtractVector(segment1_node1,segment1_node0,v);
  gridSubtractVector(segment0_node0,segment1_node0,w);

  uu = gridDotProduct(u,u);
  uv = gridDotProduct(u,v);
  vv = gridDotProduct(v,v);
  uw = gridDotProduct(u,w);
  vw = gridDotProduct(v,w);

  denom = uu*vv - uv*uv;

  if (denom < 1.0e-12) { /* colinear */
    *segment0 = 0.0;
    if ( vv > abs(uv) ) {
      *segment1 = vw/vv;
    } else {
      *segment1 = uw/uv;
    }
    return FALSE;
  }else{
    *segment0 = (uv*vw - vv*uw) / denom;
    *segment1 = (uu*vw - uv*uw) / denom;
  }
  
  return ( ( 0.0 <= *segment0 && *segment0 <= 1.0 ) &&
	   ( 0.0 <= *segment1 && *segment1 <= 1.0 ) );
}

GridBool intersectAbove( double *vertex0, double *vertex1, double *vertex2,
		     double *node )
{
  double edge1[3], edge2[3], normal[3], direction[3];
  double height;
  gridSubtractVector(vertex1, vertex0, edge1);
  gridSubtractVector(vertex2, vertex0, edge2);
  gridCrossProduct(edge1, edge2, normal);
  gridSubtractVector(node, vertex0, direction);
  height = gridDotProduct(normal,direction);
  if ( height > 0 ) return TRUE;
  return FALSE;
}

GridBool intersectTriangleNode( double *vertex0, double *vertex1, double *vertex2,
			    double *node)
{
  double side[3], otherSide[3], target[3];
  double triangleNormal[3], targetNormal[3];

  gridSubtractVector(vertex1, vertex0, side);

  gridSubtractVector(vertex2, vertex0, otherSide);
  gridCrossProduct(side, otherSide, triangleNormal);

  gridSubtractVector(node, vertex0, target);
  gridCrossProduct(side, target, targetNormal);
  if ( 0 >= gridDotProduct( triangleNormal, targetNormal) ) return FALSE;

  gridSubtractVector(vertex2, vertex1, side);
  gridSubtractVector(node, vertex1, target);
  gridCrossProduct(side, target, targetNormal);
  if ( 0 >= gridDotProduct( triangleNormal, targetNormal) ) return FALSE;

  gridSubtractVector(vertex0, vertex2, side);
  gridSubtractVector(node, vertex2, target);
  gridCrossProduct(side, target, targetNormal);
  if ( 0 >= gridDotProduct( triangleNormal, targetNormal) ) return FALSE;

  return TRUE;
}

GridBool intersectTriangleSegment(double *vertex0, double *vertex1, double *vertex2,
			      double *node0, double *node1, double *ratio )
{
  int i;
  double edge1[3], edge2[3], normal[3];
  double dir0[3], dir1[3];
  double h0, h1;
  GridBool coplanar;
  double denom, intersection[3];

  *ratio = -1.0;

  gridSubtractVector(vertex1, vertex0, edge1);
  gridSubtractVector(vertex2, vertex0, edge2);
  gridCrossProduct(edge1, edge2, normal);
  gridVectorNormalize( normal );
  gridSubtractVector(node0, vertex0, dir0);  
  gridSubtractVector(node1, vertex0, dir1);  

  h0 = gridDotProduct(normal,dir0);
  h1 = gridDotProduct(normal,dir1);

  if (h0>0 && h1>0) return FALSE;
  if (h0<0 && h1<0) return FALSE;

  coplanar = (ABS(h0)<1.0e-14 && ABS(h1)<1.0e-14);

  if (!coplanar) {
    denom = 1/(h0-h1);
    for(i=0;i<3;i++) intersection[i] = ( h0*node1[i] - h1*node0[i] )*denom;
    if ( intersectTriangleNode(vertex0, vertex1, vertex2, intersection) ) {
      *ratio = h0*denom;
      return TRUE;
    } else {
      return FALSE;
    }
  }else{
    if (intersectTriangleNode(vertex0, vertex1, vertex2, node0) ) {
      *ratio = 0.0;
      return TRUE;
    }
    if (intersectTriangleNode(vertex0, vertex1, vertex2, node1) ) {
      *ratio = 1.0;
      return TRUE;
    }
    return FALSE;
  }
}

GridBool intersectInsideTet(double *vertex0, double *vertex1, 
			double *vertex2, double *vertex3,
			double *node)
{
  if ( !intersectAbove( vertex1, vertex3, vertex2, node )) return FALSE;
  if ( !intersectAbove( vertex2, vertex3, vertex0, node )) return FALSE;
  if ( !intersectAbove( vertex3, vertex1, vertex0, node )) return FALSE;
  if ( !intersectAbove( vertex0, vertex1, vertex2, node )) return FALSE;

  return TRUE;
}

GridBool intersectTetSegment(double *vertex0, double *vertex1, 
			 double *vertex2, double *vertex3,
			 double *node0, double *node1)
{
  double ratio;
  if (intersectInsideTet( vertex0, vertex1, vertex2, vertex3, node0 )) 
    return TRUE;

  if (intersectInsideTet( vertex0, vertex1, vertex2, vertex3, node1 )) 
    return TRUE;

  if (intersectTriangleSegment( vertex1, vertex3, vertex2, node0, node1, &ratio)) 
    return TRUE;
  if (intersectTriangleSegment( vertex2, vertex3, vertex0, node0, node1, &ratio)) 
    return TRUE;
  if (intersectTriangleSegment( vertex3, vertex1, vertex0, node0, node1, &ratio)) 
    return TRUE;
  if (intersectTriangleSegment( vertex0, vertex1, vertex2, node0, node1, &ratio)) 
    return TRUE;

  return FALSE;
}

GridBool intersectTetTet(double *vertex0, double *vertex1, 
		     double *vertex2, double *vertex3,
		     double *node0, double *node1,
		     double *node2, double *node3 )
{
  if (intersectInsideTet( vertex0, vertex1, vertex2, vertex3, node0 )) 
    return TRUE;
  if (intersectInsideTet( vertex0, vertex1, vertex2, vertex3, node1 )) 
    return TRUE;
  if (intersectInsideTet( vertex0, vertex1, vertex2, vertex3, node2 )) 
    return TRUE;
  if (intersectInsideTet( vertex0, vertex1, vertex2, vertex3, node3 )) 
    return TRUE;

  if (intersectInsideTet( node0, node1, node2, node3, vertex0 )) 
    return TRUE;
  if (intersectInsideTet( node0, node1, node2, node3, vertex1 )) 
    return TRUE;
  if (intersectInsideTet( node0, node1, node2, node3, vertex2 )) 
    return TRUE;
  if (intersectInsideTet( node0, node1, node2, node3, vertex3 )) 
    return TRUE;

  if (intersectTetSegment( vertex0, vertex1, vertex2, vertex3, node0, node1 )) 
    return TRUE;
  if (intersectTetSegment( vertex0, vertex1, vertex2, vertex3, node0, node2 )) 
    return TRUE;
  if (intersectTetSegment( vertex0, vertex1, vertex2, vertex3, node0, node3 )) 
    return TRUE;
  if (intersectTetSegment( vertex0, vertex1, vertex2, vertex3, node1, node2 )) 
    return TRUE;
  if (intersectTetSegment( vertex0, vertex1, vertex2, vertex3, node1, node3 )) 
    return TRUE;
  if (intersectTetSegment( vertex0, vertex1, vertex2, vertex3, node2, node3 )) 
    return TRUE;

  if (intersectTetSegment( node0, node1, node2, node3, vertex0, vertex1 )) 
    return TRUE;
  if (intersectTetSegment( node0, node1, node2, node3, vertex0, vertex2 )) 
    return TRUE;
  if (intersectTetSegment( node0, node1, node2, node3, vertex0, vertex3 )) 
    return TRUE;
  if (intersectTetSegment( node0, node1, node2, node3, vertex1, vertex2 )) 
    return TRUE;
  if (intersectTetSegment( node0, node1, node2, node3, vertex1, vertex3 )) 
    return TRUE;
  if (intersectTetSegment( node0, node1, node2, node3, vertex2, vertex3 )) 
    return TRUE;


  return FALSE;
}

/* debug chaff
    printf("\nnode0 %10.5f %10.5f %10.5f\n",node0[0],node0[1],node0[2]);
    printf("node1 %10.5f %10.5f %10.5f\n",node1[0],node1[1],node1[2]);
    printf("norm  %10.5f %10.5f %10.5f\n",normal[0],normal[1],normal[2]);
    printf("dir0  %10.5f %10.5f %10.5f\n",dir0[0],dir0[1],dir0[2]);
    printf("dir1  %10.5f %10.5f %10.5f\n",dir1[0],dir1[1],dir1[2]);
    printf("h0 %10.5f h1 %10.5f\n",h0,h1);
    printf("isect %10.5f %10.5f %10.5f\n",
	   intersection[0],intersection[1],intersection[2]);
*/
