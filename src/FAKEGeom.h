
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef CADGEOM_H
#define CADGEOM_H

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
typedef void          *Iterator;
#endif /* HAVE_SDK */

typedef struct _UGridPtr {
   magic_t    magic;		/* Magic Number */
} UGrid, *UGridPtr;
#define UGrid_PatchList(ugp)     (NULL)

UGridPtr CADGeom_VolumeGrid( int );

typedef struct {
  magic_t   magic;
  double *param;
} CADCurve, *CADCurvePtr;
 
GridBool CADGeom_GetEdge(int, int, double *, int *);
CADCurvePtr CADGeom_EdgeGrid( int, int );

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

char *ErrMgr_GetErrStr(void);

GridBool CADGeom_Start( void );
GridBool GeoMesh_LoadPart( char *project );
GridBool CADGeom_GetVolume(int, int *, int *, int *, int *);
UGPatchPtr CADGeom_FaceGrid( int, int );
GridBool CADGeom_NormalToFace( int vol, int faceId, 
			       double *uv, double *xyz, double *normal);

END_C_DECLORATION

#endif /* CADGEOM_H */
