
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

  ring->malloc_chunk_size = 100;

  ring->segments = 0;
  ring->malloced_segments = 100;
  ring->segment_nodes = (int *)malloc(ring->malloced_segments*2*sizeof(int));
  ring->segment_uvs = 
    (double *)malloc(ring->malloced_segments*4*sizeof(double));

  ring->triangles = 0;
  ring->malloced_triangles = 100;
  ring->triangle_nodes = (int *)malloc(ring->malloced_triangles*3*sizeof(int));
  ring->triangle_uvs = 
    (double *)malloc(ring->malloced_triangles*6*sizeof(double));

  return ring;
}

void ringFree( Ring *ring )
{
  if ( NULL != ring->segment_nodes  ) free( ring->segment_nodes );
  if ( NULL != ring->segment_uvs    ) free( ring->segment_uvs );
  if ( NULL != ring->triangle_nodes ) free( ring->triangle_nodes );
  if ( NULL != ring->triangle_uvs   ) free( ring->triangle_uvs );
  free( ring );
}

Ring *ringInspect( Ring *ring )
{
  int segment;
  printf("ring inspection,%4d segments\n",ringSegments( ring ));
  for ( segment=0 ; segment<ringSegments( ring ) ; segment++ ) {
    printf("segment%4d: nodes%10d%10d\n",
	   segment,
	   ring->segment_nodes[0+2*segment],ring->segment_nodes[1+2*segment]);
  }
  return ring;
}


int ringSegments( Ring *ring )
{
  return ring->segments;
}

Ring *ringAddSegment( Ring *ring,
		      int node0, int node1,
		      double *uv0, double *uv1 )
{
  int segment;
  int remove_segment;

  for ( segment = 0 ; segment < ringSegments(ring) ; segment++ ) {
    if ( ring->segment_nodes[0+2*segment] == node0 &&
	 ring->segment_nodes[1+2*segment] == node1 ) {
      return NULL;
    }
  }

  remove_segment = EMPTY;
  for ( segment = 0 ;
	EMPTY == remove_segment && segment < ringSegments(ring) ;
	segment++ ) {
    if ( ring->segment_nodes[0+2*segment] == node1 &&
	 ring->segment_nodes[1+2*segment] == node0 ) {
      remove_segment = segment;
    }
  }

  if ( EMPTY != remove_segment ) {
    for ( segment = remove_segment ; 
	  segment < ( ringSegments(ring) - 1 ) ; 
	  segment++ ) {
      ring->segment_nodes[0+2*segment] = ring->segment_nodes[0+2*(segment+1)];
      ring->segment_nodes[1+2*segment] = ring->segment_nodes[1+2*(segment+1)];
      ring->segment_uvs[0+4*segment] = ring->segment_uvs[0+4*(segment+1)]; 
      ring->segment_uvs[1+4*segment] = ring->segment_uvs[1+4*(segment+1)]; 
      ring->segment_uvs[2+4*segment] = ring->segment_uvs[2+4*(segment+1)]; 
      ring->segment_uvs[3+4*segment] = ring->segment_uvs[3+4*(segment+1)]; 
    }
    ring->segments--;
    return ring;
  }

  if ( ring->segments >= ring->malloced_segments ) {
    ring->malloced_segments += ring->malloc_chunk_size;
    
    ring->segment_nodes = (int *)realloc(ring->segment_nodes,
					 ring->malloced_segments*2*sizeof(int));
    ring->segment_uvs =
      (double *)realloc(ring->segment_uvs,
			ring->malloced_segments*4*sizeof(double));
  }

  ring->segment_nodes[0+2*ring->segments] = node0;
  ring->segment_nodes[1+2*ring->segments] = node1;
  ring->segment_uvs[0+4*ring->segments] = uv0[0];
  ring->segment_uvs[1+4*ring->segments] = uv0[1];
  ring->segment_uvs[2+4*ring->segments] = uv1[0];
  ring->segment_uvs[3+4*ring->segments] = uv1[1];

  ring->segments++;
  return ring;
}

Ring *ringSegment( Ring *ring, int segment, int *node0, int *node1,
		   double *uv0, double *uv1 )
{
  if ( segment < 0 || segment >= ringSegments( ring ) ) return NULL;

  (*node0) = ring->segment_nodes[0+2*segment];
  (*node1) = ring->segment_nodes[1+2*segment];

  uv0[0] = ring->segment_uvs[0+4*segment];
  uv0[1] = ring->segment_uvs[1+4*segment];

  uv1[0] = ring->segment_uvs[2+4*segment];
  uv1[1] = ring->segment_uvs[3+4*segment];

  return ring;
}

GridBool ringSegmentsContainNode( Ring *ring, int node, double *uv )
{
  int segment;

  for ( segment = 0 ; segment < ringSegments(ring) ; segment++ ) {

    if ( ring->segment_nodes[0+2*segment] == node ) {
      uv[0] = ring->segment_uvs[0+4*segment];
      uv[1] = ring->segment_uvs[1+4*segment];
      return TRUE;
    }

    if ( ring->segment_nodes[1+2*segment] == node ) {
      uv[0] = ring->segment_uvs[2+4*segment];
      uv[1] = ring->segment_uvs[3+4*segment];
      return TRUE;
    }

  }

  return FALSE;
}

int ringTriangles( Ring *ring )
{
  return ring->triangles;
}

Ring *ringAddTriangle( Ring *ring,
		       int node0, int node1, int node2,
		       double *uv2 )
{
  int segment;
  double uv0[2], uv1[2];
  for ( segment = 0 ; segment < ringSegments(ring) ; segment++ ) {
    if ( ring->segment_nodes[0+2*segment] == node0 &&
	 ring->segment_nodes[1+2*segment] == node1 ) {

      uv0[0] = ring->segment_uvs[0+4*segment];
      uv0[1] = ring->segment_uvs[1+4*segment];
      uv1[0] = ring->segment_uvs[2+4*segment];
      uv1[1] = ring->segment_uvs[3+4*segment];
      
      if ( ring != ringAddSegment(ring,node1,node0,uv1,uv0) ) {
	printf("%s: %d: ringAddTriangle: ringAddSegment 1 0 returned NULL.\n",
	       __FILE__,__LINE__);
	return NULL;
      }

      if ( ring != ringAddSegment(ring,node2,node1,uv2,uv1) ) {
	printf("%s: %d: ringAddTriangle: ringAddSegment 2 1 returned NULL.\n",
	       __FILE__,__LINE__);
	return NULL;
      }

      if ( ring != ringAddSegment(ring,node0,node2,uv0,uv2) ) {
	printf("%s: %d: ringAddTriangle: ringAddSegment 0 2 returned NULL.\n",
	       __FILE__,__LINE__);
	return NULL;
      }

      if ( ring->triangles >= ring->malloced_triangles ) {
	ring->malloced_triangles += ring->malloc_chunk_size;
    
	ring->triangle_nodes = (int *)realloc(ring->triangle_nodes,
					      ring->malloced_triangles*
					      3*sizeof(int));
	ring->triangle_uvs =
	  (double *)realloc(ring->triangle_uvs,
			    ring->malloced_triangles*6*sizeof(double));
      }

      ring->triangle_nodes[0+3*ring->triangles] = node0;
      ring->triangle_nodes[1+3*ring->triangles] = node1;
      ring->triangle_nodes[2+3*ring->triangles] = node2;
      ring->triangle_uvs[0+6*ring->triangles] = uv0[0];
      ring->triangle_uvs[1+6*ring->triangles] = uv0[1];
      ring->triangle_uvs[2+6*ring->triangles] = uv1[0];
      ring->triangle_uvs[3+6*ring->triangles] = uv1[1];
      ring->triangle_uvs[4+6*ring->triangles] = uv2[0];
      ring->triangle_uvs[5+6*ring->triangles] = uv2[1];
      
      ring->triangles++;
      return ring;
    }
  }
  return NULL;
}

Ring *ringTriangle( Ring *ring, int triangle,
		    int *node0, int *node1, int *node2,
		    double *uv0, double *uv1, double *uv2 )
{
  if ( triangle < 0 || triangle >= ringTriangles( ring ) ) return NULL;

  (*node0) = ring->triangle_nodes[0+3*triangle];
  (*node1) = ring->triangle_nodes[1+3*triangle];
  (*node2) = ring->triangle_nodes[2+3*triangle];

  uv0[0] = ring->triangle_uvs[0+6*triangle];
  uv0[1] = ring->triangle_uvs[1+6*triangle];

  uv1[0] = ring->triangle_uvs[2+6*triangle];
  uv1[1] = ring->triangle_uvs[3+6*triangle];

  uv2[0] = ring->triangle_uvs[4+6*triangle];
  uv2[1] = ring->triangle_uvs[5+6*triangle];

  return ring;
}

