
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
#include "Goolache/MeshMgr.h"
#include "MeatLib/ErrMgr.h"

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

      /* start hack
      {
	double d[3],dist;
	d[0] = newxyz[0+3*(nedgenode-1)] - newxyz[0+3*(nedgenode-2)];
	d[1] = newxyz[1+3*(nedgenode-1)] - newxyz[1+3*(nedgenode-2)];
	d[2] = newxyz[2+3*(nedgenode-1)] - newxyz[2+3*(nedgenode-2)];
	dist = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
	printf("last delta %e\n",dist);
	if (dist<1.0e-8){
	  printf("WARNING: duplicate node detected and deleted");
	  nedgenode--;
	}
      }
      end hack */

      printf("number of rebuild edge points:  %d\n",nedgenode);

      newnodes = malloc( nedgenode * sizeof(int));
      newnodes[0] = edgeEndPoints[0];
      newnodes[nedgenode-1] = edgeEndPoints[1];
      for(i=1;i<(nedgenode-1);i++){
	newnodes[i]=gridAddNode(grid,newxyz[0+3*i],newxyz[1+3*i],newxyz[2+3*i]);
	printf("node added %8d x %8.5f y %8.5f z %8.5f \n",
	       newnodes[i],newxyz[0+3*i],newxyz[1+3*i],newxyz[2+3*i]);
      }

      gridDeleteThawedEdgeSegments(grid,edgeId);
      for(i=1;i<nedgenode;i++){
	i0 = i-1; i1 = i;
	gridAddEdge(grid,newnodes[i0],newnodes[i1],edgeId,newt[i0],newt[i1]);
	printf("edge added %8d <-> %8d \n",newnodes[i0],newnodes[i1]);
      }

      free(newxyz);
      free(newt);
      free(newnodes);
      
    }
  }

  return layer;
}

Layer *layerRebuildFaces(Layer *layer, int vol){

  int maxnode, nnode;
  int edgeId;
  int faceId;
  double uv[4];
  int loop, nloop;
  int *loopLength;
  int *loopEdge;
  int edge;
  int nedge;
  int nshell;
  int nparent;
  int nthaw;
  int orient;
  int *shell;
  int nodes[3], front, side;
  int ncurve, *curve;
  int i,j;
  int n0,n1;
  int *l2g, *g2l;
  int node;
  double *shellxyz, *shelluv, resolved[3];
  int nfacenode, nfacetri, *newface;
  double *newxyz, *newuv;
  

  Grid *grid;
  grid = layerGrid(layer);

  maxnode = gridMaxNode(grid);

  l2g = malloc(maxnode*sizeof(int));
  g2l = malloc(maxnode*sizeof(int));

  for (faceId=1;faceId<=gridNGeomFace(grid);faceId++){
    if (layerConstrainingGeometry(layer,faceId)) {
      printf("faceId %4d is a rebuild face.\n",faceId);    
      CADGeom_GetFace(vol, faceId, uv, &nloop, &loopLength, &loopEdge);
      nedge = 0;
      for (loop=0;loop<nloop;loop++) nedge += loopLength[loop];
      nshell =0;
      for(edge=0;edge<nedge;edge++){
	edgeId = loopEdge[0+2*edge];
	orient = loopEdge[1+2*edge];
	nparent = layerNParentEdgeSegments(layer,edgeId);
	if (nparent > 0 ) {
	  nshell += nparent;
	  printf(" edge %4d, edgeId %4d %2d has %4d phantom.\n",
		 edge,edgeId,orient,nparent);
	} else if ( layerConstrainingGeometry(layer,-edgeId) ) {
	  nthaw = gridNThawedEdgeSegments(grid,edgeId);
	  nshell += nthaw;
	  printf(" edge %4d, edgeId %4d %2d has %4d rebuild.\n",
		 edge,edgeId,orient,nthaw);
	} else {
	  nthaw = gridGeomEdgeSize(grid,edgeId)-1;
	  nshell += nthaw;
	  printf(" edge %4d, edgeId %4d %2d has %4d original.\n",
		 edge,edgeId,orient,nthaw);
	}
      }
      printf("faceId %4d has %4d segments\n",faceId,nshell);
      shell = malloc(2*nshell*sizeof(int));
      nshell =0;
      for(edge=0;edge<nedge;edge++){
	edgeId = loopEdge[0+2*edge];
	orient = loopEdge[1+2*edge];
	nparent = layerNParentEdgeSegments(layer,edgeId);
	if (nparent > 0 ) {
	  for(front=0;front<layerNFront(layer);front++){
	    for(side=0;side<3;side++){
	      if (edgeId == layerParentEdge(layer,front,side)){
		n0 = side;
		n1 = side+1; if (n1>2) n1 = 0;
		layerFront(layer,front,nodes);
		n0 = nodes[n0];
		n1 = nodes[n1];
		shell[0+2*nshell] = n0;
		shell[1+2*nshell] = n1;
		nshell++;
 	      }
	    }
	  }
	  printf("phantom %4d %4d edge added. nshell %8d\n",edge,edgeId,nshell);
	} else if ( layerConstrainingGeometry(layer,-edgeId) ) {
	  ncurve = gridGeomEdgeSize(grid,edgeId);
	  curve = malloc( ncurve * sizeof(int) );
	  gridGeomEdge( grid, edgeId, curve );
	  for(i=1;i<ncurve;i++){
	    if ( !gridNodeFrozen(grid,curve[i-1]) || 
		 !gridNodeFrozen(grid,curve[i])   ){
	      if (orient>0){
		shell[0+2*nshell] = curve[i-1];
		shell[1+2*nshell] = curve[i];
	      }else{
		shell[0+2*nshell] = curve[i];
		shell[1+2*nshell] = curve[i-1];
	      }
	      nshell++;
	    }	      
	  }
	  free(curve);
	  printf("rebuild %4d %4d edge added. nshell %8d\n",edge,edgeId,nshell);
	} else {
	  ncurve = gridGeomEdgeSize(grid,edgeId);
	  curve = malloc( ncurve * sizeof(int) );
	  gridGeomEdge( grid, edgeId, curve );
	  for(i=1;i<ncurve;i++){
	    if (orient>0){
	      shell[0+2*nshell] = curve[i-1];
	      shell[1+2*nshell] = curve[i];
	    }else{
	      shell[0+2*nshell] = curve[i];
	      shell[1+2*nshell] = curve[i-1];
	    }
	    nshell++;	      
	  }
	  free(curve);
	  printf("original%4d %4d edge added. nshell %8d\n",edge,edgeId,nshell);
	}
      }
      printf("faceId %4d has %4d segments\n",faceId,nshell);

      for(i=0;i<maxnode;i++) l2g[i]=EMPTY;
      for(i=0;i<maxnode;i++) g2l[i]=EMPTY;
      nnode = 0;
      for(i=0;i<nshell;i++){
	if (EMPTY == g2l[shell[0+2*i]] ) { 
	  g2l[shell[0+2*i]] = nnode; 
	  nnode++;
	}
	if (EMPTY == g2l[shell[1+2*i]] ) { 
	  g2l[shell[1+2*i]] = nnode; 
	  nnode++;
	}
	shell[0+2*i] = g2l[shell[0+2*i]];
	shell[1+2*i] = g2l[shell[1+2*i]];
      }
      printf("the shell has %8d nodes.\n",nnode);
      for(i=0;i<maxnode;i++)if (EMPTY != g2l[i]) l2g[g2l[i]]=i;
      shellxyz = malloc(3*nnode*sizeof(double));
      shelluv  = malloc(2*nnode*sizeof(double));
      for(i=0;i<nnode;i++) gridNodeXYZ(grid,l2g[i],&shellxyz[3*i]);
      for(i=0;i<nnode;i++) {
	gridNodeUV(grid,l2g[i],faceId, &shelluv[2*i]);
	/* project all nodes for safety. rebuilt edge has no uv */
	CADGeom_ResolveOnFace(vol,faceId,&shellxyz[3*i],&shelluv[2*i],resolved);
      }
      for(i=0;i<nnode;i++) 
	printf("node %4d %8d x %8.5f y %8.5f z %8.5f u %11.3e v %11.3e\n",
	       i,l2g[i],
	       shellxyz[0+3*i],shellxyz[1+3*i],shellxyz[2+3*i],
	       shelluv[0+2*i],shelluv[1+2*i]
	       );
      for(i=0;i<nshell;i++) 
	printf("shell %4d: %8d <-> %8d or %8d <-> %8d\n",
	       i,shell[0+2*i],shell[1+2*i],l2g[shell[0+2*i]],l2g[shell[1+2*i]]);

/*
      for(j=0;j<nshell;j++) {
        printf("2\n");
        i = shell[0+2*j];
	printf("%8.5f %8.5f %8.5f\n",shellxyz[0+3*i],shellxyz[1+3*i],shellxyz[2+3*i]);
        i = shell[1+2*j];
	printf("%8.5f %8.5f %8.5f\n",shellxyz[0+3*i],shellxyz[1+3*i],shellxyz[2+3*i]);
      }

      for(j=0;j<nshell;j++) {
        printf("2\n");
        i = shell[0+2*j];
	printf("%8.5f %8.5f 0.0\n",shelluv[0+2*i],shelluv[1+2*i]);
        i = shell[1+2*j];
	printf("%8.5f %8.5f 0.0\n",shelluv[0+2*i],shelluv[1+2*i]);
      }
 */

      nfacenode = EMPTY;
      nfacetri  = EMPTY;
      if( !MeshMgr_MeshTriFace(vol, faceId, 
			       nnode, shellxyz, shelluv,
			       nshell, shell,
			       0, NULL,
			       &nfacenode, &nfacetri, 
			       &newface, &newxyz, &newuv) ) {
	printf("%s\nCould NOT mesh Face %d\n",ErrMgr_GetErrStr(),faceId);
	//return NULL;
      }
      printf("rebuild face has %d nodes %d faces\n",nfacenode,nfacetri);

      free(shell);
    }
  }

  free(g2l);
  free(l2g);

  return layer;
}

int
MesherX_DiscretizeVolume( int npts, double *points, int ntri_b, int *tri_b,
                          int ntri, int *tri, int nsid, int *sid, int *npo,
                          int *nel, int **iel, double **xyz, char *pj)
{
  char *outputProject;
  int vol=1;
  Grid *grid;
  Layer *layer;
  int i;

  grid = gridFillFromPart( vol, npts*10 );

  layer = formAdvancingFront( grid, pj );

  /* only needed for formAdvancingFront freeze distant volume nodes */
  gridThawAll(grid); 
  layerFindParentEdges(layer);

  for (i=0;i<5;i++) layerAdvance(layer,0.01);

  printf(" -- REBUILD EDGES\n");
  layerRebuildEdges(layer,vol);

  printf(" -- REBUILD FACES\n");
  layerRebuildFaces(layer,vol);

  printf(" -- DUMP PART\n");
  outputProject = "../test/MesherX";
  printf("writing DEBUG output project %s\n",outputProject);
  gridSavePart( grid, outputProject );

  return 1;
}
