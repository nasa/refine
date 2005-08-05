
/* Ring, continuous loop of linear segements on a geometry face that
 *       is meshed by filling with triangles.
 *
 * Michael A. Park
 * Computational AeroSciences Branch
 * NASA Langley Research Center
 * Phone: (757) 864-6604
 * Email: Mike.Park@NASA.Gov
 */

/* $Id$ */

#include <stdlib.h>
#ifndef __APPLE__       /* Not needed on Mac OS X */
#include <malloc.h>
#endif
#include "ring.h"
#include "sort.h"

Ring* ringCreate( void )
{
  Ring *ring;

  ring = (Ring *)malloc( sizeof(Ring) );

  ring->segments  = 0;
  ring->triangles = 0;

  return ring;
}

void ringFree( Ring *ring )
{
  free( ring );
}

int ringSegments( Ring *ring )
{
  return ring->segments;
}

Ring *ringAddSegment( Ring *ring,
		      int node0, int node1,
		      double *uv0, double *uv1 )
{
  ring->segments++;
  return ring;
}

int ringTriangles( Ring *ring )
{
  return ring->triangles;
}
