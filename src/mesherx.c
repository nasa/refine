
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

Layer *layerRebuildEdges(Layer *layer, int vol){

  int i, edgeId, edgeEndPoints[2];

  double edgexyz[6],tRange[2];
  int nedgenode;
  double *newxyz;
  double *newt;
  int *newnodes;
  int i0, i1;

  Grid *grid;

  grid = layerGrid(layer);

  for (edgeId=1;edgeId<=gridNGeomEdge(grid);edgeId++) {
    if ( layerConstrainingGeometry(layer,-edgeId) ){
      edgeEndPoints[0]=gridGeomEdgeStart(grid,edgeId);
      edgeEndPoints[1]=gridGeomEdgeEnd(grid,edgeId);
      printf("rebuild edge %4d:  %10d <-> %10d\n",
	     edgeId,edgeEndPoints[0],edgeEndPoints[1]);
      edgeEndPoints[0] = gridFrozenEdgeEndPoint(grid,edgeId,edgeEndPoints[0]);
      edgeEndPoints[1] = gridFrozenEdgeEndPoint(grid,edgeId,edgeEndPoints[1]);
      printf("rebuild endpoints:  %10d <-> %10d\n",
	     edgeEndPoints[0],edgeEndPoints[1]);
      gridNodeXYZ(grid, edgeEndPoints[0], &edgexyz[0]);
      gridNodeXYZ(grid, edgeEndPoints[1], &edgexyz[3]);
      gridNodeT(grid, edgeEndPoints[0], edgeId, &tRange[0]);
      gridNodeT(grid, edgeEndPoints[1], edgeId, &tRange[1] );

      /* WTJ will change this to generic with DSO */
      if( !MeshMgr_MeshEdge(vol, edgeId, edgexyz, tRange, &nedgenode, &newxyz, &newt) ) {
	printf("Could NOT mesh Edge %d\n",edgeId);
	return NULL;
      }

      printf("number of rebuild edge points:  %d\n",nedgenode);

      newnodes = malloc( nedgenode * sizeof(int));
      newnodes[0] = edgeEndPoints[0];
      newnodes[nedgenode-1] = edgeEndPoints[1];
      for(i=1;i<(nedgenode-1);i++){
	newnodes[i]=gridAddNode(grid,newxyz[0+3*i],newxyz[1+3*i],newxyz[2+3*i]);
      }

      gridDeleteThawedEdgeSegments(grid,edgeId);
      for(i=1;i<(nedgenode-1);i++){
	i0 = i-1; i1 = i;
	gridAddEdge(grid,newnodes[i0],newnodes[i1],edgeId,newt[i0],newt[i1]);
      }
      i0 = nedgenode-2; i1 = nedgenode-1;
      gridAddEdge(grid,newnodes[i0],newnodes[i1],edgeId,newt[i0],newt[i1]);

      free(newxyz);
      free(newt);
      free(newnodes);
      
    }
  }

  return layer;
}

int
MesherX_DiscretizeVolume( int npts, double *points, int ntri_b, int *tri_b,
                          int ntri, int *tri, int nsid, int *sid, int *npo,
                          int *nel, int **iel, double **xyz)
{
  char *outputProject;
  int vol=1;
  Grid *grid;
  Layer *layer;
  int i;

  int edgeId;
  int faceId;
  double uv[4];
  int loop, nloop;
  int *loopLength;
  int *loopEdge;
  int edge;
  int nedge;
  int orient;

  grid = gridFillFromPart( vol, npts*10 );

  layer = formAdvancingFront( grid, "box" );

  /* only needed for formAdvancingFront freeze distant volume nodes */
  gridThawAll(grid); 
  layerFindParentEdges(layer);

  for (i=0;i<5;i++) layerAdvance(layer,0.01);

  printf(" -- REBUILD EDGES\n");
  layerRebuildEdges(layer,vol);

  printf(" -- REBUILD FACES\n");
  for (faceId=1;FALSE && faceId<=gridNGeomFace(grid);faceId++){
    CADGeom_GetFace(vol, faceId, uv, &nloop, &loopLength, &loopEdge);
    printf("face %d has %d loops\n",faceId,nloop);
    for (loop=0;loop<nloop;loop++){
      printf("  loop %d has %d edges\n",loop,loopLength[loop]);
      printf("   ");
      for(edge=0;edge<loopLength[loop];edge++){
	printf(" edge %d dir %d",loopEdge[0+2*edge],loopEdge[1+2*edge]);
      }
      printf("\n");
    }
  }

  for (faceId=1;faceId<=gridNGeomFace(grid);faceId++){
    if (layerConstrainingGeometry(layer,faceId)) {
      printf("faceId %d is a rebuild face.\n",faceId);    
      CADGeom_GetFace(vol, faceId, uv, &nloop, &loopLength, &loopEdge);
      nedge = 0;
      for (loop=0;loop<nloop;loop++) nedge += loopLength[loop];
      for(edge=0;edge<nedge;edge++){
	edgeId = loopEdge[0+2*edge];
	orient = loopEdge[1+2*edge];
	if (layerConstrainingGeometry(layer,-edgeId)) {
	  printf(" edge %4d, edgeId %4d %2d is rebuild.\n",edge,edgeId,orient);
	} else if ( layerNParentEdgeSegments(layer,edgeId)>0 ) {
	  printf(" edge %4d, edgeId %4d %2d is phantom.\n",edge,edgeId,orient);
	} else {
	  printf(" edge %4d, edgeId %4d %2d is original.\n",edge,edgeId,orient);
	}
      }
    }
  }

  printf(" -- DUMP PART\n");
  outputProject = "../test/MesherX";
  printf("writing DEBUG output project %s\n",outputProject);
  gridSavePart( grid, outputProject );

  return 1;
}
