
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  


#include <stdlib.h>
#include <stdio.h>
#include "FAKEGeom.h"
  
GridBool UGrid_FromArrays(UGridPtr *ugp,
			  int a,double *b,int c,int *d,int e,int *f)
{
  SUPRESS_UNUSED_COMPILER_WARNING(ugp);
  SUPRESS_UNUSED_COMPILER_WARNING(a);
  SUPRESS_UNUSED_COMPILER_WARNING(b);
  SUPRESS_UNUSED_COMPILER_WARNING(c);
  SUPRESS_UNUSED_COMPILER_WARNING(d);
  SUPRESS_UNUSED_COMPILER_WARNING(e);
  SUPRESS_UNUSED_COMPILER_WARNING(f);
  return FALSE;
}

UGridPtr CADGeom_VolumeGrid( int vol )
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  return NULL;
}

GridBool CADGeom_GetEdge(int a, int b, double *c, int *d) {
  SUPRESS_UNUSED_COMPILER_WARNING(a);
  SUPRESS_UNUSED_COMPILER_WARNING(b);
  SUPRESS_UNUSED_COMPILER_WARNING(c);
  SUPRESS_UNUSED_COMPILER_WARNING(d);
  return FALSE;
}

GridBool CADGeom_GetFace(int a, int b, double *c, int *d, int **e, int **f) {
  SUPRESS_UNUSED_COMPILER_WARNING(a);
  SUPRESS_UNUSED_COMPILER_WARNING(b);
  SUPRESS_UNUSED_COMPILER_WARNING(c);
  SUPRESS_UNUSED_COMPILER_WARNING(d);
  *e = NULL;
  *f = NULL;
  return FALSE;
}

CADCurvePtr CADGeom_EdgeGrid( int a, int b)
{
  SUPRESS_UNUSED_COMPILER_WARNING(a);
  SUPRESS_UNUSED_COMPILER_WARNING(b);
  return NULL;
}

GridBool CADGeom_UpdateEdgeGrid(int vol, int iedge, int nCurveNode,
				double *xyz, double *t)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(iedge);
  SUPRESS_UNUSED_COMPILER_WARNING(nCurveNode);
  SUPRESS_UNUSED_COMPILER_WARNING(xyz);
  SUPRESS_UNUSED_COMPILER_WARNING(t);
  return FALSE;
}

GridBool CADGeom_UpdateFaceGrid(int vol, int faceId, int nnode,
				double *xyz, double *t,
				int nface, int *f2n)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(faceId);
  SUPRESS_UNUSED_COMPILER_WARNING(nnode);
  SUPRESS_UNUSED_COMPILER_WARNING(xyz);
  SUPRESS_UNUSED_COMPILER_WARNING(t);
  SUPRESS_UNUSED_COMPILER_WARNING(nface);
  SUPRESS_UNUSED_COMPILER_WARNING(f2n);
  return FALSE;
}

CADCurvePtr CADGeom_EdgeGrid( int, int );

void *DList_SetIteratorToHead(DListPtr dlp,Iterator *dli)
{
  SUPRESS_UNUSED_COMPILER_WARNING(dlp);
  SUPRESS_UNUSED_COMPILER_WARNING(dli);
  return NULL;
}

void *DList_GetNextItem(Iterator *dli)
{
  SUPRESS_UNUSED_COMPILER_WARNING(dli);
  return NULL;
}

void UGPatch_GetDims(UGPatchPtr upp, int *dims)
{
  SUPRESS_UNUSED_COMPILER_WARNING(upp);
  dims[0]=EMPTY;
  dims[1]=EMPTY;
  dims[2]=EMPTY;
}

int UGPatch_GlobalIndex(UGPatchPtr upp, int ndx)
{
  SUPRESS_UNUSED_COMPILER_WARNING(upp);
  SUPRESS_UNUSED_COMPILER_WARNING(ndx);
  return EMPTY;
}

GridBool UGPatch_InitSurfacePatches(UGridPtr ugp)
{
  SUPRESS_UNUSED_COMPILER_WARNING(ugp);
  return FALSE;
}

char *ErrMgr_GetErrStr(void)
{
  return "Non-functional FAKEGeomExtras\n";
}
void ErrMgr_Append(char *mod,int line,char *fmt,...)
{
  SUPRESS_UNUSED_COMPILER_WARNING(mod);
  SUPRESS_UNUSED_COMPILER_WARNING(line);
  SUPRESS_UNUSED_COMPILER_WARNING(fmt);
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
GridBool GeoMesh_LoadModel( char *url, char *modeler, char *project, int *mod )
#else
GridBool GeoMesh_LoadPart( char *project )
#endif
{
  SUPRESS_UNUSED_COMPILER_WARNING(project);
  return FALSE;
}

#ifdef HAVE_CAPRI2
GridBool CADGeom_LoadModel( char *url, char *modeler, char *project, int *mod )
#else
GridBool CADGeom_LoadPart( char *project )
#endif
{
  SUPRESS_UNUSED_COMPILER_WARNING(project);
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
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(project);
  return FALSE;
}

GridBool CADGeom_GetVolume(int a, int *b, int *c, int *d, int *e)
{
  SUPRESS_UNUSED_COMPILER_WARNING(a);
  SUPRESS_UNUSED_COMPILER_WARNING(b);
  SUPRESS_UNUSED_COMPILER_WARNING(c);
  SUPRESS_UNUSED_COMPILER_WARNING(d);
  SUPRESS_UNUSED_COMPILER_WARNING(e);
  return FALSE;
}

UGPatchPtr CADGeom_FaceGrid( int vol, int face )
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(face);
  return NULL;
}

GridBool CADGeom_SetFaceGrid( int vol, int faceId, UGPatchPtr ugrid)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(faceId);
  SUPRESS_UNUSED_COMPILER_WARNING(ugrid);
  return FALSE;
}

GridBool CADTopo_FaceNumEdgePts(int vol, int faceId, int *count)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(faceId);
  SUPRESS_UNUSED_COMPILER_WARNING(count);
  return FALSE;
}

GridBool CADTopo_VolFacePts(int vol, int faceId, int *count, int *l2g)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(faceId);
  SUPRESS_UNUSED_COMPILER_WARNING(count);
  SUPRESS_UNUSED_COMPILER_WARNING(l2g);
  return FALSE;
}

GridBool CADTopo_VolEdgePts(int vol, int *count)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(count);
  return FALSE;
}

GridBool CADTopo_ShellStats(int vol, int *nc, int *tPts, int *tTri, int *maxFac)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(nc);
  SUPRESS_UNUSED_COMPILER_WARNING(tPts);
  SUPRESS_UNUSED_COMPILER_WARNING(tTri);
  SUPRESS_UNUSED_COMPILER_WARNING(maxFac);
  return FALSE;
}

UGridPtr CADTopo_AssembleTShell(int vol,int tPts, int tTri, int maxFace)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(tPts);
  SUPRESS_UNUSED_COMPILER_WARNING(tTri);
  SUPRESS_UNUSED_COMPILER_WARNING(maxFace);
  return NULL;
}

GridBool CADGeom_DisplacementIsIdentity(int vol)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  return TRUE;
}

GridBool CADGeom_UnMapPoint(int vol, double *xyz, double *pt)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(xyz);
  SUPRESS_UNUSED_COMPILER_WARNING(pt);
  return FALSE;
}

GridBool CADGeom_MapPoint(int vol, double *xyz, double *pt)
{
  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(xyz);
  SUPRESS_UNUSED_COMPILER_WARNING(pt);
  return FALSE;
}

#define DOTP(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])
#define CROSS(a, b, r)                                          \
{                                                               \
  ((r)[0]) = ((a)[1]) * ((b)[2]) - ((a)[2]) * ((b)[1]);         \
  ((r)[1]) = ((a)[2]) * ((b)[0]) - ((a)[0]) * ((b)[2]);         \
  ((r)[2]) = ((a)[0]) * ((b)[1]) - ((a)[1]) * ((b)[0]);         \
}
GridBool CADGeom_ReversedSurfaceNormal(int vol, int face)
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

GridBool CADGeom_ResolveOnEdgeWCS(int vol, int edge, double *coor, double *t,
                                  double *point)
{
  return CADGeom_NearestOnEdge(vol,edge,coor,t,point);
}

GridBool CADGeom_ResolveOnFaceWCS(int vol, int face, double *coor, double *uv,
                                  double *point)
{
  return CADGeom_NearestOnFace(vol,face,coor,uv,point);
}
