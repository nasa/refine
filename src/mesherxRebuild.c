
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>         /* Needed in some systems for DBL_MAX definition */
#include <float.h>
#include "CADGeom/CADGeom.h"
#include "CADGeom/CADTopo.h"
#include "Goolache/MeshMgr.h"
#include "Goolache/UGPatch.h"
#include "MeatLib/ErrMgr.h"
#include "grid.h"
#include "mesherxRebuild.h"

Layer *layerRebuildInterior(Layer *layer, int vol)
{

  printf(" -- REBUILD EDGES\n");
  if ( layer != layerRebuildEdges(layer,vol) ) return NULL;

  printf(" -- REBUILD FACES\n");
  if ( layer != layerRebuildFaces(layer,vol) ) return NULL;

  printf(" -- REBUILD VOLUME\n");
  if ( layer != layerRebuildVolume(layer,vol) ) return NULL;

  return layer;
}

Layer *layerRebuildEdges(Layer *layer, int vol){

  int i, edgeId, edgeEndPoints[2];

  double edgexyz[6],tRange[2];
  int nedgenode;
  double *newxyz;
  double *newt;
  int *newnodes;
  int i0, i1;
  int edge, nodes[2], id;

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
      if( !MeshMgr_MeshEdge(vol, edgeId, edgexyz, tRange, 
			    &nedgenode, &newxyz, &newt) ) {
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

      for(edge=0;edge<gridMaxEdge(grid);edge++){
	if( (grid==gridEdge(grid,edge,nodes,&id)) && 
	    (id == edgeId) && 
	    !layerEdgeInLayer(layer,edge) ){
	  gridRemoveEdge(grid,edge);
	}
      }

      for(i=1;i<nedgenode;i++){
	i0 = i-1; i1 = i;
	gridAddEdge(grid,newnodes[i0],newnodes[i1],edgeId,newt[i0],newt[i1]);
      }

      free(newxyz);
      free(newt);
      free(newnodes);
      
    }
  }

  return layer;
}

void layerDumpFaceWire(char *fileroot, int faceId, int nshell, 
		       int *shell, double *shellxyz){
  FILE *mfile;
  char filename[256];
  int i;

  printf("Dumping face wire shell.\n");
  sprintf(filename,"%s%d.plt",fileroot,faceId);
  mfile = fopen(filename,"w");
  fprintf(mfile,"Variables = \"x\",\"y\",\"z\"\n");
  for(i=0;i<nshell;i++){
    fprintf(mfile,"Zone\n%20.10f %20.10f %20.10f\n%20.10f %20.10f %20.10f\n",
	    shellxyz[0+3*shell[0+2*i]],
	    shellxyz[1+3*shell[0+2*i]],
	    shellxyz[2+3*shell[0+2*i]],
	    shellxyz[0+3*shell[1+2*i]],
	    shellxyz[1+3*shell[1+2*i]],
	    shellxyz[2+3*shell[1+2*i]]);
  }
  fclose(mfile);
}

void layerDumpTecplotShell( char *filename, int nnode, int nshell, 
			    double *xyz, int *shell )
{
  int i;
  FILE *tecplotFile;
  
  printf("writing shell to tecplot file %s\n",filename);

  tecplotFile = fopen(filename,"w");
  fprintf(tecplotFile, "title=\"tecplot refine geometry file\"\n");
  fprintf(tecplotFile, "variables=\"X\",\"Y\",\"Z\"\n");
  fprintf(tecplotFile, "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	  nnode, nshell);
  for ( i=0; i<nnode ; i++ ){
    fprintf(tecplotFile, "%23.15e%23.15e%23.15e\n",
	    xyz[0+3*i],xyz[1+3*i],xyz[2+3*i]);
  }
  fprintf(tecplotFile, "\n");
  for ( i=0; i<nshell ; i++ ){
    fprintf(tecplotFile, " %9d %9d %9d\n",
	    shell[0+3*i]+1,shell[1+3*i]+1,shell[2+3*i]+1);
  }
  fclose(tecplotFile);
}

Layer *layerRebuildFaces(Layer *layer, int vol){

  int maxnode, nnode;
  int edgeId;
  int faceId, id, face;
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
  int nodes[3], triangle, side;
  int ncurve, *curve;
  int i,j;
  int n0,n1;
  int *l2g, *g2l;
  int node;
  double *shellxyz, *shelluv, resolved[3];
  int nfacenode, nfacetri, *newface;
  double *newxyz, *newuv;
  int l0,l1,l2,g0,g1,g2;

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
	nparent = layerNParentGeomEdgeSegments(layer,edgeId);
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
	nparent = layerNParentGeomEdgeSegments(layer,edgeId);
	if (nparent > 0 ) {
	  for(triangle=0;triangle<layerNTriangle(layer);triangle++){
	    for(side=0;side<3;side++){
	      if (edgeId == layerParentGeomEdge(layer,triangle,side)){
		n0 = side;
		n1 = side+1; if (n1>2) n1 = 0;
		layerTriangle(layer,triangle,nodes);
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
	shelluv[2*i]=DBL_MAX;
	shelluv[2*i+1]=DBL_MAX;
	gridNodeUV(grid,l2g[i],faceId, &shelluv[2*i]);
	/* project all nodes for safety. rebuilt edge has no uv */
	CADGeom_ResolveOnFace(vol,faceId,&shellxyz[3*i],&shelluv[2*i],resolved);
      }

      layerDumpFaceWire("faceWireOrig", faceId, nshell, shell, shellxyz);

      nfacenode = EMPTY;
      nfacetri  = EMPTY;
      if( !MeshMgr_MeshTriFace(vol, faceId, 
			       nnode, shellxyz, shelluv,
			       nshell, shell,
			       0, NULL,
			       &nfacenode, &nfacetri, 
			       &newface, &newxyz, &newuv) ) {
	printf("rebuild face %d has FAILED.\n", faceId );
	return NULL;
      }
      printf("rebuild face has %d nodes %d faces\n",nfacenode,nfacetri);
      
      if (TRUE) {
	char filename[256];
	sprintf(filename,"rebuildFace%d.plt",faceId);
	layerDumpTecplotShell( filename, nfacenode, nfacetri, 
			       newxyz,  newface);
      }
      for(face=0;face<gridMaxFace(grid);face++){
	if( (grid==gridFace(grid,face,nodes,&id)) && 
	    (id == faceId) && 
	    !layerFaceInLayer(layer,face) ){
	  gridRemoveFace(grid,face);
	}
      }

      for(i=nnode;i<nfacenode;i++){
	l2g[i]=gridAddNode(grid,newxyz[0+3*i],newxyz[1+3*i],newxyz[2+3*i]);
      }
      for(i=0;i<nfacetri;i++){
	l0 = newface[0+3*i]; 
	l1 = newface[1+3*i]; 
	l2 = newface[2+3*i];
	g0 = l2g[l0];
	g1 = l2g[l1];
	g2 = l2g[l2];
	gridAddFaceUV(grid,
		      g0, newuv[0+2*l0], newuv[1+2*l0],
		      g1, newuv[0+2*l1], newuv[1+2*l1],
		      g2, newuv[0+2*l2], newuv[1+2*l2],
		      faceId );
      }

      free(shell);
      free(newface);
      free(newxyz);
      free(newuv);
    }
  }

  free(g2l);
  free(l2g);

  return layer;
}

Layer *layerRebuildVolume(Layer *layer, int vol){

  int faceId;
  int nshell, *shell;
  int *l2g, *g2l;
  int maxnode, nnode;
  int maxface, face;
  int nodes[3];
  int triangle;
  int i, cell;
  double *shellxyz;

  int nvolnode, nvolcell;
  int *newcell;
  double *newxyz;

  Grid *grid;
  grid = layerGrid(layer);

  maxnode = gridMaxNode(grid);
  maxface = gridMaxFace(grid);

  l2g = malloc(maxnode*sizeof(int));
  g2l = malloc(maxnode*sizeof(int));

  nshell =0;
  for (face=0;face<maxface;face++){
    if (grid==gridFace(grid,face,nodes,&faceId)){
      if (!layerFaceInLayer(layer,face) ){
	nshell +=1;
      }
    }
  }
  nshell += layerNTriangle(layer);
  printf("allocating  %10d faces for shell.\n",nshell);

  shell = malloc(3*nshell*sizeof(int));

  nshell =0;
  for (face=0;face<maxface;face++){
    if (grid==gridFace(grid,face,nodes,&faceId)){
      if (!layerFaceInLayer(layer,face) ){
	shell[0+3*nshell] = nodes[0];
	shell[1+3*nshell] = nodes[1];
	shell[2+3*nshell] = nodes[2];
	nshell++;
      }
    }
  }
  for(triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangle(layer, triangle, nodes);
    shell[0+3*nshell] = nodes[0];
    shell[1+3*nshell] = nodes[1];
    shell[2+3*nshell] = nodes[2];
    nshell++;
  }
  printf("inserted  %12d faces into shell.\n",nshell);

  for(i=0;i<maxnode;i++) l2g[i]=EMPTY;
  for(i=0;i<maxnode;i++) g2l[i]=EMPTY;
  nnode = 0;
  for(i=0;i<nshell;i++){
    if (EMPTY == g2l[shell[0+3*i]] ) { 
      g2l[shell[0+3*i]] = nnode; 
      nnode++;
    }
    if (EMPTY == g2l[shell[1+3*i]] ) { 
      g2l[shell[1+3*i]] = nnode; 
      nnode++;
    }
    if (EMPTY == g2l[shell[2+3*i]] ) { 
      g2l[shell[2+3*i]] = nnode; 
      nnode++;
    }
    shell[0+3*i] = g2l[shell[0+3*i]];
    shell[1+3*i] = g2l[shell[1+3*i]];
    shell[2+3*i] = g2l[shell[2+3*i]];
  }
  printf("the shell has %8d nodes.\n",nnode);
  for(i=0;i<maxnode;i++)if (EMPTY != g2l[i]) l2g[g2l[i]]=i;
  shellxyz = malloc(3*nnode*sizeof(double));
  for(i=0;i<nnode;i++) gridNodeXYZ(grid,l2g[i],&shellxyz[3*i]);

  layerDumpTecplotShell("debugVolumeShell.t", nnode, nshell, shellxyz, shell);

  if( !MeshMgr_MeshTetVolume(vol, 
			     nnode, shellxyz,
			     &nshell, &shell,
			     0, NULL, 0, NULL,
			     &nvolnode, &nvolcell, 
			     &newcell, &newxyz) ) {

    printf("%s\nCould NOT mesh Volume %d\n",ErrMgr_GetErrStr(),vol);

    layerDumpTecplotShell("failedVolumeShell.t",nnode,nshell,shellxyz,shell);

    return NULL;
  }
  printf("rebuild volume has %d nodes %d cells\n",nvolnode,nvolcell);

  for(cell=0;cell<gridMaxCell(grid);cell++) 
    if ( !layerCellInLayer(layer,cell) ) gridRemoveCell(grid,cell);
  gridDeleteNodesNotUsed(grid);

  for(i=nnode;i<nvolnode;i++){
    l2g[i]=gridAddNode(grid,newxyz[0+3*i],newxyz[1+3*i],newxyz[2+3*i]);
  }
  for(i=0;i<nvolcell;i++){
    gridAddCell(grid,
		l2g[newcell[0+4*i]],
		l2g[newcell[1+4*i]],
		l2g[newcell[2+4*i]],
		l2g[newcell[3+4*i]]);
  }

  free(l2g);
  free(g2l);
  free(newcell);
  free(newxyz);

  return layer;
}
