
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
#include "UG_API/UG_API.h"
#else
#include "FAKEGeom.h"
#endif

Grid *gridParallelGeomLoad( Grid *grid, char *project )
{
  int vol=1;
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;
  UGridPtr ugrid;
  int i, iedge, inode;
  int nedgenode;
  double trange[2];
  int edgeEndPoint[2];
  CADCurvePtr edge;
  int face, localNode, globalNode, partitionNode;
  int patchDimensions[3];
  UGPatchPtr  localPatch, globalPatch;
  Iterator patchIterator;

  if ( ! MeshMgr_Initialize( ) ){
    printf("ERROR: MeshMgr_Initialize broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  

  if ( ! CADGeom_Start( ) ){
    printf("ERROR: CADGeom_Start broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  

  if ( ! GeoMesh_LoadPart( project ) ){
    printf("ERROR: GeoMesh_LoadPart broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  if( !CADGeom_GetVolume(vol,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) ) {
    printf("ERROR: CADGeom_GetVolume. \n%s\n",ErrMgr_GetErrStr());
  }

  gridSetNGeomNode( grid, nGeomNode );
  gridSetNGeomEdge( grid, nGeomEdge );
  gridSetNGeomFace( grid, nGeomFace );

  inode = nGeomNode;

  for( iedge=1; iedge<=nGeomEdge; iedge++ ) {
    if( (edge=CADGeom_EdgeGrid(vol,iedge)) == NULL ) 
      printf("ERROR: CADGeom_EdgeGrid(%d).\n%s\n",iedge,ErrMgr_GetErrStr());
 
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

  /* get uv vals for surface(s) */
  /* we use globalPatch to track with the localPatch so that we can get global
   * node numbering relative the volume grid and NOT the face grid as would
   * be the case of global index of upp
   */

  if (NULL == (ugrid = CADGeom_VolumeGrid(vol)) ) {
    printf("ERROR: Can not find grid in restart. \n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  globalPatch = DList_SetIteratorToHead(UGrid_PatchList(ugrid),&patchIterator);

  for( face=1; face<=nGeomFace; face++ ) {
    localPatch = CADGeom_FaceGrid(vol,face);
    UGPatch_GetDims(localPatch,patchDimensions);
    for( localNode=0; localNode<patchDimensions[0]; localNode++ ) {
      globalNode = UGPatch_GlobalIndex(globalPatch,localNode);
      partitionNode = gridGlobal2Local(grid, globalNode);
      if (partitionNode>EMPTY) 
	gridSetNodeUV( grid, partitionNode, face,
		       UGPatch_Parameter(localPatch,localNode,0), 
		       UGPatch_Parameter(localPatch,localNode,1));
    }
    globalPatch = DList_GetNextItem(&patchIterator);
  }

  return grid;
}

Grid *gridParallelGeomSave( Grid *grid, char *project )
{
  int vol=1;

  GeoMesh_UseDefaultIOCallbacks();
  if( !CADGeom_SavePart(vol,project) ) {
    printf("%s: %d: Could not save CAPRI part.\n",
	   __FILE__, __LINE__);    
    return NULL;
  }

  return grid;
}

Grid *gridUpdateEdgeGrid( Grid *grid, int edgeId, int nCurveNode, 
			  double *xyz, double *t )
{
  int vol=1;
  if (!CADGeom_UpdateEdgeGrid( vol, edgeId, nCurveNode, xyz, t)) return NULL;
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
  UGridPtr ugrid;  
  int iface;
  
  if ( !UGrid_FromArrays( &ugrid, nnode, xyz, nface, f2n, 0, NULL  )) {
    printf("%s: %d: Could not make UGrid_FromArrays.\n", __FILE__, __LINE__);
    return NULL;
  }

  for (iface = 0 ; iface < nface ; iface++ ) {
    UGrid_FlagValue(ugrid,iface) = faceId;
  }

  if( !UGPatch_InitSurfacePatches(ugrid) ) {
    printf("%s: %d: Could not make surface patches for new UGridPtr.\n",
	   __FILE__, __LINE__);    
    return NULL;
  }

  printf("WARNING %s: %d: set patch UV param's.\n", __FILE__, __LINE__);

  UGrid_TIMESTAMP(ugrid) = time( NULL );	/* Updated time */
  UGrid_ALGORITHM(ugrid) = UG_UNKNOWN;

  if( !CADGeom_SetFaceGrid(vol,faceId,ugrid) ) {
    printf("%s: %d: Could not set CAPRI face UGridPtr.\n",
	   __FILE__, __LINE__);    
    return NULL;
  }

  return grid;
}
