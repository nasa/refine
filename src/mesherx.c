
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>         /* Needed in some systems for DBL_MAX definition */
#include <float.h>
#include "gridfiller.h"
#include "grid.h"
#include "layer.h"
#include "CADGeom/CADGeom.h"

/******************** EXTERNAL FUNCTIONS ******************************/


int
MesherX_DiscretizeVolume( int npts, double *points, int ntri_b, int *tri_b,
                          int ntri, int *tri, int nsid, int *sid, int *npo,
                          int *nel, int **iel, double **xyz)
{
  int vol=1;
  Grid *grid;
  Layer *layer;
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;
  CADCurvePtr *phantomEdge;
  UGPatchPtr  *phantomFace;

  grid = gridFillFromPart( vol, npts*10 );

  layer = formAdvancingFront( grid, "box" );

  layerAdvance(layer,0.01);

/* Allocate Phantom Edge and Face Arrays */

  if( !CADGeom_GetVolume(1,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) ) {
    printf("ERROR: CADGeom_GetVolume, line %d of %s\n.",__LINE__, __FILE__);
    return 1;
  }

  if( (phantomEdge=(CADCurvePtr *)calloc(nGeomEdge,sizeof(CADCurvePtr))) == NULL
 ) {
    printf("ERROR: Allocation of Phantom Edges\n");
    return 1;
  }

  if( (phantomFace=(UGPatchPtr *)calloc(nGeomFace,sizeof(UGPatchPtr))) == NULL ) {
    printf("ERROR: Allocation of Phantom Faces\n");
    return 1;
  }

  return 1;
}
