
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "FAKEGeom.h"
  
UGridPtr CADGeom_VolumeGrid( int vol )
{
  return NULL;
}

GridBool CADGeom_GetEdge(int a, int b, double *c, int *d) {
  c = NULL;
  d = NULL;
  return FALSE;
}

CADCurvePtr CADGeom_EdgeGrid( int a, int b)
{
  return NULL;
}

CADCurvePtr CADGeom_EdgeGrid( int, int );

void *DList_SetIteratorToHead(DListPtr dlp,Iterator *dli)
{
  return NULL;
}

void *DList_GetNextItem(Iterator *dli)
{
  return NULL;
}

void UGPatch_GetDims(UGPatchPtr upp, int *dims)
{
  dims[0]=EMPTY;
  dims[1]=EMPTY;
  dims[2]=EMPTY;
}

int UGPatch_GlobalIndex(UGPatchPtr upp, int ndx)
{
  return EMPTY;
}

char *ErrMgr_GetErrStr(void)
{
  return "Non-functional FAKEGeomExtras\n";
}

GridBool CADGeom_Start( void )
{
  return FALSE;
}

GridBool GeoMesh_LoadPart( char *project )
{
  return FALSE;
}

GridBool CADGeom_GetVolume(int a, int *b, int *c, int *d, int *e)
{
  return FALSE;
}

UGPatchPtr CADGeom_FaceGrid( int vol, int face )
{
  return FALSE;
}

GridBool CADGeom_NormalToFace( int vol, int faceId, 
			       double *uv, double *xyz, double *normal)
{
  return FALSE;
}
