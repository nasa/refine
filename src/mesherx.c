
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>         /* Needed in some systems for DBL_MAX definition */
#include <float.h>
#include "gridfiller.h"
#include "grid.h"
#include "layer.h"
#include "CADGeom/CADGeom.h"


/******************** Private Functions ******************************/

static CADCurvePtr *makePhantomEdges(int vol, int nGeomEdge, Layer *layer)
{
  Grid *grid;
  int normal, edgeId, globalNode;
  double direction[3],edgexyz[3],dumxyz[3];
  double t;
  double tangent[3],curv;
  double trange[2];
  int    nodes[2];
  CADCurvePtr *phantomEdge;
  
  grid = layerGrid(layer);

  if( (phantomEdge=(CADCurvePtr *)calloc(nGeomEdge,sizeof(CADCurvePtr))) == NULL
 ) {
    printf("ERROR: Allocation of Phantom Edges\n");
    return NULL;
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    edgeId = -layerConstrained(layer, normal );
    if (edgeId > 0){	/* Constrained to Edge */
      layerNormalDirection(layer, normal, direction);
      globalNode = layerNormalRoot(layer, normal );
      gridNodeXYZ(grid,globalNode,edgexyz);
      gridNodeT(grid,globalNode,edgeId, &t);

/* Remove when gridNodeT() returns good t value */
CADGeom_ResolveOnEdge(vol,edgeId,edgexyz,&t,dumxyz);

      if( !CADGeom_CurvOfEdge(vol,edgeId,t,tangent,&curv) ) {
        printf("ERROR: Tangent Broke\n");
        return NULL;
      }

      if( phantomEdge[edgeId-1] == NULL ) {
        if( (phantomEdge[edgeId-1]=CADCurve_New(2)) == NULL ) {
          printf("ERROR: Allocation of Phantom Edge CADCurve\n");
          return NULL;
        }

        CADGeom_GetEdge(vol,edgeId,trange,nodes);
        CADGeom_GetNode(vol,nodes[0],CADCurve_Point(phantomEdge[edgeId-1],0));
        CADCurve_Param(phantomEdge[edgeId-1],0) = trange[0];
        CADGeom_GetNode(vol,nodes[1],CADCurve_Point(phantomEdge[edgeId-1],1));
        CADCurve_Param(phantomEdge[edgeId-1],1) = trange[1];
      }

      if( gridDotProduct(direction,tangent) > 0.0 ) {
        CADCurve_Coord(phantomEdge[edgeId-1],0,X) = edgexyz[0];
        CADCurve_Coord(phantomEdge[edgeId-1],0,Y) = edgexyz[1];
        CADCurve_Coord(phantomEdge[edgeId-1],0,Z) = edgexyz[2];
        CADCurve_Param(phantomEdge[edgeId-1],0) = t;
      } else {
        CADCurve_Coord(phantomEdge[edgeId-1],1,X) = edgexyz[0];
        CADCurve_Coord(phantomEdge[edgeId-1],1,Y) = edgexyz[1];
        CADCurve_Coord(phantomEdge[edgeId-1],1,Z) = edgexyz[2];
        CADCurve_Param(phantomEdge[edgeId-1],1) = t;
      }

    }
  }

  return phantomEdge;
}


/******************** Public Functions ******************************/


int
MesherX_DiscretizeVolume( int npts, double *points, int ntri_b, int *tri_b,
                          int ntri, int *tri, int nsid, int *sid, int *npo,
                          int *nel, int **iel, double **xyz)
{
  int vol=1;
  Grid *grid;
  Layer *layer;
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;
  CADCurvePtr *phantomEdge;
  UGPatchPtr  *phantomFace;
  int normal;

  grid = gridFillFromPart( vol, npts*10 );

  layer = formAdvancingFront( grid, "box" );

  layerAdvance(layer,0.01);

/* Allocate Phantom Edge and Face Arrays */

  if( !CADGeom_GetVolume(1,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) ) {
    printf("ERROR: CADGeom_GetVolume, line %d of %s\n.",__LINE__, __FILE__);
    return 1;
  }

  if( (phantomEdge=makePhantomEdges(vol,nGeomEdge,layer)) == NULL ) {
    printf("ERROR: Could NOT create Phantom Edges line %d of %s\n.",__LINE__, __FILE__);
    return 1;
  }

  if( (phantomFace=(UGPatchPtr *)calloc(nGeomFace,sizeof(UGPatchPtr))) == NULL ) {
    printf("ERROR: Allocation of Phantom Faces\n");
    return 1;
  }

  return 1;
}
