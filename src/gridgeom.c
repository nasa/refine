
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
#include <CADGeom/CADGeom.h>

Grid *gridParallelGeomLoad( Grid *grid, char *project )
{
  int vol=1;
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;

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

  return grid;

}
