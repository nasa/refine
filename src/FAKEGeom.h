
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef CADGEOM_H
#define CADGEOM_H

#include <stdlib.h>
#include <time.h>
#include "refine_defs.h"

BEGIN_C_DECLORATION

GridBool CADGeom_NearestOnEdge(int vol, int edgeId, 
			   double *xyz, double *t, double *xyznew);
GridBool CADGeom_NearestOnFace(int vol, int faceId, 
			   double *xyz, double *uv, double *xyznew);

GridBool CADGeom_PointOnEdge(int vol, int edgeId, 
			 double t, double *xyz, 
			 int derivativeFlag, double *dt, double *dtdt );
GridBool CADGeom_PointOnFace(int vol, int faceId, 
			 double *uv, double *xyz, 
			 int derivativeFlag, double *du, double *dv,
			 double *dudu, double *dudv, double *dvdv );

GridBool CADGeom_NormalToFace( int vol, int faceId, 
			       double *uv, double *xyz, double *normal);

#ifdef HAVE_SDK
#else
typedef unsigned char magic_t;		/* Magic Number type */
typedef unsigned char algo_t;          /* Computed Algorithm type */
typedef void          *Iterator;
#endif /* HAVE_SDK */

typedef struct _UGridPtr {
   magic_t    magic;		/* Magic Number */
   int flags;
   time_t     updated;		/* Last Updated Timestamp */
   algo_t     algorithm;	/* Algorithm used to compute grid */
  double *param;
} UGrid, *UGridPtr;
#define UGrid_FlagValue(ugp,i)   ((ugp)->flags)
#define UGrid_PatchList(ugp)     (NULL)
#define UGrid_TIMESTAMP(ugp)     ((ugp)->updated)
#define UGrid_ALGORITHM(ugp)     ((ugp)->algorithm)

GridBool UGrid_FromArrays(UGridPtr *,int,double *,int,int *,int,int *);

UGridPtr CADGeom_VolumeGrid( int );

extern int UGrid_GetDims(UGridPtr ugp, int *dims);
extern int UGrid_BuildConnectivity(UGridPtr);

#define UGrid_PtValue(ugp,i,l)   (ugp)->param[l]
#define UGrid_VertValue(ugp,i,l)   (ugp)->param[l]
#define UGrid_TetValue(ugp,i,l)   (ugp)->param[l]

typedef struct {
  magic_t   magic;
  double *param;
} CADCurve, *CADCurvePtr;
 
GridBool CADGeom_GetEdge(int, int, double *, int *);
CADCurvePtr CADGeom_EdgeGrid( int, int );
/* to suppress unused warn */
#define CADCURVE_NUMPTS(edge) (0==(edge)->magic?-1:(edge)->magic)

GridBool CADGeom_UpdateEdgeGrid(int, int, int, double *, double *);
GridBool CADGeom_UpdateFaceGrid(int, int, int, double *, double *, int, int *);

typedef struct _DList {
   magic_t    magic;
} DList,*DListPtr;

void *DList_SetIteratorToHead(DListPtr dlp,Iterator *dli);
void *DList_GetNextItem(Iterator *dli);

typedef struct _UGPatchPtr {
   magic_t    magic;    /* Magic Number */
   double *param;
} UGPatch, *UGPatchPtr;

#define BC_NOSLIP     0	

void UGPatch_GetDims(UGPatchPtr upp, int *dims);
int UGPatch_GlobalIndex(UGPatchPtr upp, int ndx);
#define UGPatch_Parameter(upp,i,l) ((upp)->param[i*l])
#define UGPatch_Parent(upp) ((UGridPtr)malloc(sizeof(UGrid)))

GridBool UGPatch_InitSurfacePatches(UGridPtr ugp);
void ErrMgr_Append(char *mod,int line,char *fmt,...);

char *ErrMgr_GetErrStr(void);

GridBool UGMgr_LoadLibs( void );
GridBool CADGeom_Start( void );
#ifdef HAVE_CAPRI2
GridBool GeoMesh_LoadModel( char *url, char *modeler, char *project, int *mdl );
GridBool CADGeom_LoadModel( char *url, char *modeler, char *project, int *mdl );
#else
GridBool GeoMesh_LoadPart( char *project );
GridBool CADGeom_LoadPart( char *project );
#endif

GridBool CADGeom_SavePart(int vol, char *project);

extern GridBool CADGeom_SetVolumeGrid( int, UGridPtr );

GridBool CADGeom_GetVolume(int, int *, int *, int *, int *);
UGPatchPtr CADGeom_FaceGrid( int, int );
GridBool CADGeom_SetFaceGrid( int vol, int faceId, UGPatchPtr ugrid);
GridBool CADGeom_NormalToFace( int vol, int faceId, 
			       double *uv, double *xyz, double *normal);

extern GridBool CADTopo_FlushPatches(int, UGridPtr);
GridBool CADTopo_FaceNumEdgePts(int vol, int faceId, int *count);
GridBool CADTopo_VolFacePts(int vol, int faceId, int *count, int *l2g);
GridBool CADTopo_VolEdgePts(int vol, int *count);
GridBool CADTopo_ShellStats(int vol, int *nc, int *tPts, int *tTri, int *maxF);
UGridPtr CADTopo_AssembleTShell(int vol,int tPts, int tTri, int maxFace);
GridBool CADGeom_DisplacementIsIdentity(int vol);
GridBool CADGeom_UnMapPoint(int vol, double *xyz, double *pt);
GridBool CADGeom_MapPoint(int vol, double *xyz, double *pt);

GridBool CADGeom_ReversedSurfaceNormal(int vol, int face);
GridBool CADGeom_ResolveOnEdgeWCS(int vol, int edge, double *coor, double *t,
                                  double *point);
GridBool CADGeom_ResolveOnFaceWCS(int vol, int face, double *coor, double *uv,
                                  double *point);

void GeoMesh_UseDefaultIOCallbacks( void );
void CADGeom_UseDefaultIOCallbacks( void );

#define UG_SUCCESS           0


END_C_DECLORATION

#endif /* CADGEOM_H */
