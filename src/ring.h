
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

#ifndef RING_H
#define RING_H

#include "refine_defs.h"
#include <stdio.h>

BEGIN_C_DECLORATION

typedef struct Ring Ring;

struct Ring {

  int malloc_chunk_size;

  int segments;
  int malloced_segments;
  int *segment_nodes;
  double *segment_uvs;

  int triangles;
  int malloced_triangles;
  int *triangle_nodes;
  double *triangle_uvs;

  /* originals are for debugging */
  int originals;
  int *original_nodes;
  double *original_uvs;
};

Ring *ringCreate( void );
void ringFree( Ring * );

Ring *ringInspect( Ring * );
Ring *ringTecplot( Ring *, char *filename );

int ringSegments( Ring * );
Ring *ringAddSegment( Ring *, int node0, int node1, double *uv0, double *uv1 );
Ring *ringSegment( Ring *, int segment, int *node0, int *node1,
		   double *uv0, double *uv1 );
GridBool ringSegmentsContainNode( Ring *, int node, double *uv );
GridBool ringIntersectsSegment( Ring *, int node0, int node1,
				double *uv0, double *uv1 );

int ringTriangles( Ring * );
Ring *ringAddTriangle( Ring *, int node0, int node1, int node2, double *uv2 );
Ring *ringTriangle( Ring *, int triangle, int *node0, int *node1, int *node2,
		    double *uv0, double *uv1, double *uv2 );
GridBool ringSurroundsTriangle( Ring *,
				int node0, int node1, int node2, double *uv2 );

double ringArea( Ring * );

END_C_DECLORATION

#endif /* RING_H */
