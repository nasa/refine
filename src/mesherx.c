
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
#include "grid.h"
#include "layer.h"

/******************** EXTERNAL FUNCTIONS ******************************/


int
MesherX_DiscretizeVolume( int npts, double *points, int ntri_b, int *tri_b,
                          int ntri, int *tri, int nsid, int *sid, int *npo,
                          int *nel, int **iel, double **xyz)
{
  int vol=1;
  Grid *grid;
  Layer *layer;

  grid = gridFillFromPart( vol, npts*10 );

  layer = formAdvancingFront( grid, "box" );

  layerAdvance(layer,0.01);

  /* Make edge and face lists */

  return 1;
}
