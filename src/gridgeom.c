
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
#else
#include "FAKEGeomExtras.h"
#endif

Grid *gridParallelGeomLoad( Grid *grid, char *project )
{
  int vol=1;
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;
  UGridPtr ugrid;
  UGPatchPtr  localPatch, globalPatch;
  Iterator patchIterator;
  int face, localNode, globalNode, partitionNode;
  int patchDimensions[3];


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

  printf("Geometry: %d nodes %d edges %d faces %d boundaries\n",
	 nGeomNode,nGeomEdge,nGeomFace,nGeomGroups);

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
