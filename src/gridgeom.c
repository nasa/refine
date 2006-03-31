
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "gridgeom.h"
#ifdef HAVE_SDK
#include "CADGeom/CADGeom.h"
#include "CADGeom/CADTopo.h"
#else
#include "FAKEGeom.h"
#endif

#define UG_REFINE (6) /* Lifted from SDK/UG_API/UG_API.h to remove dependency */

Grid *gridParallelGeomLoad( Grid *grid, char *url, char *modeler,
                            char *project )
{
  int vol=1;
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;
  int i, iedge, inode;
  int nedgenode;
  double trange[2];
  int edgeEndPoint[2];
  CADCurvePtr edge;
  int face, localNode, globalNode, partitionNode;
  int patchDimensions[3];
  int volumeEdgeNode, patchEdgeNode;
  int *patch2global;
  UGPatchPtr localPatch;
  Iterator patchIterator;

  if ( ! CADGeom_Start( ) ){
    printf("ERROR: CADGeom_Start broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  

#ifdef HAVE_CAPRI2
  if ( ! CADGeom_LoadModel( url, modeler, project, &(grid->model) ) ){
#else
  if ( ! CADGeom_LoadPart( project ) ){
#endif
    printf("ERROR: CADGeom_LoadPart broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  if( !CADGeom_GetVolume(vol,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) ) {
    printf("ERROR: CADGeom_GetVolume. \n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  gridSetNGeomNode( grid, nGeomNode );
  gridSetNGeomEdge( grid, nGeomEdge );
  gridSetNGeomFace( grid, nGeomFace );

  inode = nGeomNode;

  for( iedge=1; iedge<=nGeomEdge; iedge++ ) {
    if( (edge=CADGeom_EdgeGrid(vol,iedge)) == NULL ) {
      printf("ERROR: CADGeom_EdgeGrid(%d).\n%s\n",iedge,ErrMgr_GetErrStr());
      return NULL;
    }
 
    nedgenode = CADCURVE_NUMPTS(edge);

    CADGeom_GetEdge( vol, iedge, trange, edgeEndPoint );

    edgeEndPoint[0]--; /* convert from fortran to c numbers */
    edgeEndPoint[1]--;

    gridAddGeomEdge( grid, iedge, edgeEndPoint[0], edgeEndPoint[1]);

    if (nedgenode == 2) {
      gridAddEdgeInGlobal(grid, edgeEndPoint[0], edgeEndPoint[1], 
			  iedge, trange[0], trange[1]);
    }else{
      gridAddEdgeInGlobal(grid, edgeEndPoint[0], inode, iedge,
			  edge->param[0], edge->param[1]);
      for( i=1 ; i < (nedgenode-2) ; i++ ) { // skip end segments  
	gridAddEdgeInGlobal(grid, inode, inode+1, iedge,
			    edge->param[i], edge->param[i+1]);
	inode++;
      }
      gridAddEdgeInGlobal(grid, inode, edgeEndPoint[1], iedge,
			  edge->param[nedgenode-2], 
			  edge->param[nedgenode-1]);
      inode++;
    }
  }

  if (!CADTopo_VolEdgePts( vol, &volumeEdgeNode )){
    printf("%s: %d: CADTopo_VolEdgePts failied.\n",__FILE__, __LINE__);
    return NULL;
  }

  if (inode!=volumeEdgeNode) {
    printf("%s: %d: Edge node count error %d %d.\n",
	   __FILE__, __LINE__, inode, volumeEdgeNode);
    inode = volumeEdgeNode;
    return NULL;
  }

  for( face=1; face<=nGeomFace; face++ ) {
    localPatch = CADGeom_FaceGrid(vol,face);
    UGPatch_GetDims(localPatch,patchDimensions);
    patch2global = malloc(patchDimensions[0]*sizeof(int));
    if (!CADTopo_VolFacePts(vol, face, patch2global, &patchEdgeNode)) {
      printf("%s: %d: CADTopo_VolFacePts failied.\n",__FILE__, __LINE__);
      free(patch2global);
      return NULL;
    }
    for( localNode=patchEdgeNode;localNode<patchDimensions[0]; localNode++ ) {
      patch2global[localNode] = inode; inode++;
    }
    for( localNode=0; localNode<patchDimensions[0]; localNode++ ) {
      globalNode = patch2global[localNode];
      partitionNode = gridGlobal2Local(grid, globalNode);
      if (partitionNode>EMPTY) 
	gridSetNodeUV( grid, partitionNode, face,
		       UGPatch_Parameter(localPatch,localNode,0), 
		       UGPatch_Parameter(localPatch,localNode,1));
    }
    free(patch2global);
  }

  return grid;
}

Grid *gridParallelGeomSave( Grid *grid, char *project )
{
  int vol=1;

  CADGeom_UseDefaultIOCallbacks();
#ifdef HAVE_CAPRI2
  if( !CADGeom_SaveModel(grid->model,project) ) {
    printf("%s: %d: Could not save CAPRI model.\n",
	   __FILE__, __LINE__);    
#else
  if( !CADGeom_SavePart(vol,project) ) {
    printf("%s: %d: Could not save CAPRI part.\n",
	   __FILE__, __LINE__);    
#endif
    return NULL;
  }

  return grid;
}

Grid *gridUpdateEdgeGrid( Grid *grid, int edgeId, int nCurveNode, 
			  double *xyz, double *t )
{
  int vol=1;
  int i;
  double *xyzCopy, *tCopy;

  xyzCopy = malloc(3*nCurveNode*sizeof(double));
  for(i=0;i<3*nCurveNode;i++) xyzCopy[i] = xyz[i];

  tCopy = malloc(nCurveNode*sizeof(double));
  for(i=0;i<nCurveNode;i++) tCopy[i] = t[i];

  if (!CADGeom_UpdateEdgeGrid( vol, edgeId, 
			       nCurveNode, xyzCopy, tCopy ) ) {
    free(xyzCopy);
    free(tCopy);
    return NULL;
  }

  return grid;
}

int gridFaceEdgeCount( Grid *grid, int faceId )
{
  int vol=1;
  int count;
  if (!CADTopo_FaceNumEdgePts(vol, faceId, &count)) {
    printf("%s: %d: CADTopo_FaceNumEdgePts failied.\n",__FILE__, __LINE__);
    return EMPTY;
  }
  return count;
}

Grid *gridFaceEdgeLocal2Global( Grid *grid, int faceId, 
				int faceEdgeCount, int *local2global )
{
  int vol=1;
  int count;

  if (!CADTopo_VolEdgePts( vol, &count )){
    printf("%s: %d: CADTopo_VolEdgePts failied.\n",__FILE__, __LINE__);
    return NULL;
  }

  if (!CADTopo_VolFacePts(vol, faceId, local2global, &count)) {
    printf("%s: %d: CADTopo_VolFacePts failied.\n",__FILE__, __LINE__);
    return NULL;
  }
  if (count > faceEdgeCount) {
    printf("%s: %d: ran out of array length %d %d\n",
	   __FILE__, __LINE__, count, faceEdgeCount);
    return NULL;
  }
  return grid;
}

Grid *gridUpdateGeometryFace( Grid *grid, int faceId, 
			      int nnode, double *xyz, double *uv, 
			      int nface, int *f2n )
{
  int vol=1;
  int i;
  double *xyzCopy, *uvCopy;
  int *f2nCopy;
  
  xyzCopy = malloc(3*nnode*sizeof(double));
  for(i=0;i<3*nnode;i++) xyzCopy[i] = xyz[i];

  uvCopy = malloc(2*nnode*sizeof(double));
  for(i=0;i<2*nnode;i++) uvCopy[i] = uv[i];

  f2nCopy = malloc(3*nface*sizeof(int));
  for(i=0;i<3*nface;i++) f2nCopy[i] = f2n[i];

  if ( !CADGeom_UpdateFaceGrid( vol, faceId, nnode, xyzCopy, uvCopy, 
				nface, f2nCopy )) {
    printf("%s: %d: CADGeom_UpdateFaceGrid Failed.\n", __FILE__, __LINE__);
    return NULL;
  }

  UGrid_TIMESTAMP(UGPatch_Parent(CADGeom_FaceGrid(vol,faceId))) = time( NULL );
  UGrid_ALGORITHM(UGPatch_Parent(CADGeom_FaceGrid(vol,faceId))) = UG_REFINE;

  return grid;
}

Grid *gridCreateShellFromFaces( Grid *grid )
{
  int vol=1;
  int nc;               /* Total Edge points (Not used here) */
  int tPts;             /* Total UNIQUE Points in Volume Shell */
  int tTri;             /* Total Triangles in Volume Shell */
  int maxFace=0;        /* Number of Pts of largest Face */
  UGridPtr ugp;

  if( !CADTopo_ShellStats(vol,&nc,&tPts,&tTri,&maxFace) ) {
    printf("%s: %d: Could NOT total Shell Statistics\n",__FILE__,__LINE__);
    return NULL;
  }

  if( (ugp=CADTopo_AssembleTShell(vol,tPts,tTri,maxFace)) == NULL ) {
    printf("%s: %d: Failed to create Shell from faces\n",__FILE__,__LINE__);
    return NULL;
  }

  UGrid_TIMESTAMP(ugp) = time( NULL );
  UGrid_ALGORITHM(ugp) = UG_REFINE;

  return grid;
}
