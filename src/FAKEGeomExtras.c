
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include "FAKEGeomExtras.h"
  
char *ErrMgr_GetErrStr(void)
{
  return "Non-functional FAKEGeomExtras";
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

GridBool CADGeom_NormalToFace( int vol, int faceId, 
			       double *uv, double *xyz, double *normal)
{
  return FALSE;
}
