
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:mike.park@nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "CADGeom/CADGeom.h"
#include "Goolache/CAPrI_IO.h"

#define SubtractVector(v1,v2,result) \
(result)[0] = (v1)[0] - (v2)[0]; \
(result)[1] = (v1)[1] - (v2)[1]; \
(result)[2] = (v1)[2] - (v2)[2];

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  int i;
  char project[256];
  UGridPtr ugrid;  
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;
  UGPatchPtr  localPatch, globalPatch;
  Iterator patchIterator;
  CADCurvePtr edge;
  int inode, iedge, nedgenode;
  double trange[2];
  int edgeEndPoint[2];

  int patchDimensions[3];
  int face, localNode, globalNode;
  double pxyz[3], gxyz[3], dxyz[3];
  int vol;
  vol=1;

  i = 1;
  while( i < argc ) {
    if( strcmp(argv[i],"-p") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-p argument %d: %s\n",i, project);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  printf("calling CADGeom_Start ... \n");
  if ( ! CADGeom_Start( ) ){
    printf("ERROR: CADGeom_Start broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  

  printf("calling CADGeom_Load for project <%s> ... \n",project);
  if ( ! GeoMesh_LoadPart( project ) ){
    printf("ERROR: GeoMesh_LoadPart broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  if (NULL == (ugrid = CADGeom_VolumeGrid(vol)) ) {
    printf("ERROR: Can not find grid in restart. \n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  if( !CADGeom_GetVolume(vol,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) ) {
    printf("ERROR: CADGeom_GetVolume. \n%s\n",ErrMgr_GetErrStr());
  }

  inode = nGeomNode;
  for( iedge=1; iedge<=nGeomEdge; iedge++ ) {
    if( (edge=CADGeom_EdgeGrid(vol,iedge)) == NULL ) 
      printf("ERROR: CADGeom_EdgeGrid(%d).\n%s\n",iedge,ErrMgr_GetErrStr());
 
    nedgenode = CADCURVE_NUMPTS(edge);

    CADGeom_GetEdge( vol, iedge, trange, edgeEndPoint );

    edgeEndPoint[0]--; /* convert from fortran to c numbers */
    edgeEndPoint[1]--;

    printf("e%3d local%6d global%6d\n",iedge,0,edgeEndPoint[0]);
    for( i=1 ; i < (nedgenode-1) ; i++ ) { // skip end segments  
      printf("e%3d local%6d global%6d\n",iedge,i,inode);
      inode++;
    }
    printf("e%3d local%6d global%6d\n",iedge,nedgenode-1,edgeEndPoint[1]);
  }

  globalPatch = DList_SetIteratorToHead(UGrid_PatchList(ugrid),&patchIterator);

  for( face=1; face<=nGeomFace; face++ ) {
    localPatch = CADGeom_FaceGrid(vol,face);
    UGPatch_GetDims(localPatch,patchDimensions);
    for( localNode=0; localNode<patchDimensions[0]; localNode++ ) {
      globalNode = UGPatch_GlobalIndex(globalPatch,localNode);
      gxyz[0] = UGrid_PtValue(ugrid,globalNode,0);
      gxyz[1] = UGrid_PtValue(ugrid,globalNode,1);
      gxyz[2] = UGrid_PtValue(ugrid,globalNode,2);
      pxyz[0] = UGPatch_PtValue(localPatch,localNode,0);
      pxyz[1] = UGPatch_PtValue(localPatch,localNode,1);
      pxyz[2] = UGPatch_PtValue(localPatch,localNode,2);
      SubtractVector( pxyz, gxyz, dxyz);
      printf("f%3d p%6d%10.5f%10.5f%10.5f g%6d%10.5f%10.5f%10.5f\n",
	     face, 
	     localNode, pxyz[0], pxyz[1], pxyz[2], 
	     globalNode, gxyz[0], gxyz[1], gxyz[2]);
    }
    
    globalPatch = DList_GetNextItem(&patchIterator);
  }

}
