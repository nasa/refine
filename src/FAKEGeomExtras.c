
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "FAKEGeom.h"
  
GridBool UGrid_FromArrays(UGridPtr *ugp,
			  int a,double *b,int c,int *d,int e,int *f)
{
  return FALSE;
}

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

GridBool CADGeom_UpdateEdgeGrid(int vol, int iedge, int nCurveNode,
				double *xyz, double *t)
{
  return FALSE;
}

GridBool CADGeom_UpdateFaceGrid(int vol, int faceId, int nnode,
				double *xyz, double *t,
				int nface, int *f2n)
{
  return FALSE;
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

GridBool UGPatch_InitSurfacePatches(UGridPtr ugp)
{
  return FALSE;
}

char *ErrMgr_GetErrStr(void)
{
  return "Non-functional FAKEGeomExtras\n";
}

GridBool UGMgr_LoadLibs( void )
{
  return FALSE;
}

GridBool CADGeom_Start( void )
{
  return FALSE;
}

#ifdef HAVE_CAPRI2
GridBool GeoMesh_LoadPart( char *modeler, char *project )
#else
GridBool GeoMesh_LoadPart( char *project )
#endif
{
  return FALSE;
}

#ifdef HAVE_CAPRI2
GridBool CADGeom_LoadPart( char *modeler, char *project )
#else
GridBool CADGeom_LoadPart( char *project )
#endif
{
  return FALSE;
}

void GeoMesh_UseDefaultIOCallbacks( void )
{
}

void CADGeom_UseDefaultIOCallbacks( void )
{
}

GridBool CADGeom_SavePart( int vol, char *project )
{
  return FALSE;
}

GridBool CADGeom_GetVolume(int a, int *b, int *c, int *d, int *e)
{
  return FALSE;
}

UGPatchPtr CADGeom_FaceGrid( int vol, int face )
{
  return NULL;
}

GridBool CADGeom_SetFaceGrid( int vol, int faceId, UGPatchPtr ugrid)
{
  return FALSE;
}

GridBool CADTopo_FaceNumEdgePts(int vol, int faceId, int *count)
{
  return FALSE;
}

GridBool CADTopo_VolFacePts(int vol, int faceId, int *count, int *l2g)
{
  return FALSE;
}

GridBool CADTopo_VolEdgePts(int vol, int *count)
{
  return FALSE;
}

GridBool CADTopo_ShellStats(int vol, int *nc, int *tPts, int *tTri, int *maxFac)
{
  return FALSE;
}

UGridPtr CADTopo_AssembleTShell(int vol,int tPts, int tTri, int maxFace)
{
  return NULL;
}

GridBool CADGeom_DisplacementIsIdentity(int vol)
{
  return TRUE;
}

GridBool CADGeom_UnMapPoint(int vol, double *xyz, double *pt)
{
  return FALSE;
}

GridBool CADGeom_MapPoint(int vol, double *xyz, double *pt)
{
  return FALSE;
}

#define DOTP(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define CROSS(a, b, r)                                          \
{                                                               \
  ((r)[0]) = ((a)[1]) * ((b)[2]) - ((a)[2]) * ((b)[1]);         \
  ((r)[1]) = ((a)[2]) * ((b)[0]) - ((a)[0]) * ((b)[2]);         \
  ((r)[2]) = ((a)[0]) * ((b)[1]) - ((a)[1]) * ((b)[0]);         \
}
GeoBool CADGeom_ReversedSurfaceNormal(int vol, int face)
{
  double     pt[3];
  double     du[3],dv[3];
  double     sNorm[3],fNorm[3];
  double     uv[2]={0.0,0.0};	/* Irrelevant anything should work for planes */

  /* Surface normal */
  if( !CADGeom_PointOnFace(vol,face,uv,pt,1,du,dv,NULL,NULL,NULL) ) {
    ErrMgr_Append(__FILE__,__LINE__,"Could not get Volume %d, Face %d surface derivatives\n",vol,face);
    return( FALSE );
  }
  CROSS(du,dv,sNorm);

  /* Face normal */
  if( !CADGeom_NormalToFace(vol,face,uv,pt,fNorm) ) {
    ErrMgr_Append(__FILE__,__LINE__,"Could not get Volume %d, Face %d normal\n",vol,face);
    return( FALSE );
  }
  if( DOTP(sNorm,fNorm) < 0.0 ) return( TRUE );

  return( FALSE );
}
