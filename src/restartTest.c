
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
  int patchDimensions[3];
  int face, localNode, globalNode;

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

  globalPatch = DList_SetIteratorToHead(UGrid_PatchList(ugrid),&patchIterator);

  for( face=1; face<=nGeomFace; face++ ) {
	double xyz[3], pxyz[3], gxyz[3], dxyz[3], txyz[3];
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
	printf("%2d n%5d%5d%10.2e%10.2e%10.2e\n",
	       face, localNode, globalNode,
	       dxyz[0],dxyz[1],dxyz[2]);
    }

    globalPatch = DList_GetNextItem(&patchIterator);
  }

}
