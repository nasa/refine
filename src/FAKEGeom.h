
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
} UGrid, *UGridPtr;
#define UGrid_FlagValue(ugp,i)   ((ugp)->flags)
#define UGrid_PatchList(ugp)     (NULL)
#define UGrid_TIMESTAMP(ugp)     ((ugp)->updated)
#define UGrid_ALGORITHM(ugp)     ((ugp)->algorithm)

GridBool UGrid_FromArrays(UGridPtr *,int,double *,int,int *,int,int *);
int UGrid_BuildConnectivity(UGridPtr);

UGridPtr CADGeom_VolumeGrid( int );

typedef struct {
  magic_t   magic;
  double *param;
} CADCurve, *CADCurvePtr;
 
GridBool CADGeom_GetEdge(int, int, double *, int *);
CADCurvePtr CADGeom_EdgeGrid( int, int );
#define CADCURVE_NUMPTS(edge) (-1)

GridBool CADGeom_UpdateEdgeGrid(int, int, int, double *, double *);

typedef struct _DList {
   magic_t    magic;
} DList,*DListPtr;

void *DList_SetIteratorToHead(DListPtr dlp,Iterator *dli);
void *DList_GetNextItem(Iterator *dli);

typedef struct _UGPatchPtr {
   magic_t    magic;    /* Magic Number */
} UGPatch, *UGPatchPtr;

void UGPatch_GetDims(UGPatchPtr upp, int *dims);
int UGPatch_GlobalIndex(UGPatchPtr upp, int ndx);
#define UGPatch_Parameter(upp,i,l) (DBL_MAX)

GridBool UGPatch_InitSurfacePatches(UGridPtr ugp);

char *ErrMgr_GetErrStr(void);

GridBool CADGeom_Start( void );
GridBool GeoMesh_LoadPart( char *project );
GridBool CADGeom_GetVolume(int, int *, int *, int *, int *);
UGPatchPtr CADGeom_FaceGrid( int, int );
GridBool CADGeom_NormalToFace( int vol, int faceId, 
			       double *uv, double *xyz, double *normal);

GridBool CADTopo_FaceNumEdgePts(int vol, int faceId, int *count);
GridBool CADTopo_VolFacePts(int vol, int faceId, int *count, int *l2g);
GridBool CADTopo_VolEdgePts(int vol, int *count);

END_C_DECLORATION

#endif /* CADGEOM_H */
