
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
#include <stdio.h>
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

  /* originals are for debugging */
  ring->originals = 0;
  ring->original_nodes = NULL;
  ring->original_uvs   = NULL;

  return ring;
}

void ringFree( Ring *ring )
{
  if ( NULL != ring->segment_nodes  ) free( ring->segment_nodes );
  if ( NULL != ring->segment_uvs    ) free( ring->segment_uvs );
  if ( NULL != ring->triangle_nodes ) free( ring->triangle_nodes );
  if ( NULL != ring->triangle_uvs   ) free( ring->triangle_uvs );
  /* originals are for debugging */
  if ( NULL != ring->original_nodes ) free( ring->original_nodes );
  if ( NULL != ring->original_uvs   ) free( ring->original_uvs );
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

Ring *ringTecplot( Ring *ring, char *filename )
{
  FILE *file;
  int original, segment, triangle;
  if (NULL == filename) {
    file = fopen("ring_triangles_in_uv.t","w");
  }else{
    file = fopen(filename,"w");
  } 
  fprintf(file, "title=\"tecplot ring triangle uv geometry file\"\n");
  fprintf(file, "variables=\"U\",\"V\",\"Node\"\n");

  /* originals are for debugging */
  for ( original = 0 ; original < ring->originals ; original++ ) {
    fprintf(file, "zone t=orig, i=3, j=1, f=fepoint, et=triangle\n");
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->original_uvs[0+4*original], ring->original_uvs[1+4*original],
	    ring->original_nodes[0+2*original]);
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->original_uvs[2+4*original], ring->original_uvs[3+4*original],
	    ring->original_nodes[1+2*original]);
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->original_uvs[2+4*original], ring->original_uvs[3+4*original],
	    ring->original_nodes[1+2*original]);
    fprintf(file, "1 2 3\n");
  }

  for ( segment = 0 ; segment < ringSegments( ring ) ; segment++ ) {
    fprintf(file, "zone t=seg, i=3, j=1, f=fepoint, et=triangle\n");
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->segment_uvs[0+4*segment], ring->segment_uvs[1+4*segment],
	    ring->segment_nodes[0+2*segment]);
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->segment_uvs[2+4*segment], ring->segment_uvs[3+4*segment],
	    ring->segment_nodes[1+2*segment]);
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->segment_uvs[2+4*segment], ring->segment_uvs[3+4*segment],
	    ring->segment_nodes[1+2*segment]);
    fprintf(file, "1 2 3\n");
  }

  for ( triangle = 0 ; triangle < ringTriangles( ring ) ; triangle++ ) {
    fprintf(file, "zone t=tri, i=3, j=1, f=fepoint, et=triangle\n");
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->triangle_uvs[0+6*triangle], ring->triangle_uvs[1+6*triangle],
	    ring->triangle_nodes[0+3*triangle]);
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->triangle_uvs[2+6*triangle], ring->triangle_uvs[3+6*triangle],
	    ring->triangle_nodes[1+3*triangle]);
    fprintf(file, "%23.15e%23.15e%10d\n",
	    ring->triangle_uvs[4+6*triangle], ring->triangle_uvs[5+6*triangle],
	    ring->triangle_nodes[2+3*triangle]);
    fprintf(file, "1 2 3\n");
  }

  fclose(file);
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

GridBool ringIntersectsSegment( Ring *ring,
				int node0, int node1,
				double *uv0, double *uv1 )
{
  int segment;
  int node2, node3;
  double uv2[2], uv3[2];
  double denom, ratio_a, ratio_b;
  double tol;
  GridBool meet00, meet01, meet10, meet11;

  tol = 1.0e-14;

  for ( segment = 0 ; segment < ringSegments(ring) ; segment++ ) {
    if ( ring->segment_nodes[0+2*segment] == node0 &&
	 ring->segment_nodes[1+2*segment] == node1 ) {
      return TRUE;
    }
    if ( ring->segment_nodes[0+2*segment] == node1 &&
	 ring->segment_nodes[1+2*segment] == node0 ) {
      return FALSE;
    }
  }

  for ( segment = 0 ; segment < ringSegments(ring) ; segment++ ) {

    ringSegment( ring, segment, &node2, &node3, uv2, uv3 );

    denom = (uv3[1] - uv2[1])*(uv1[0] - uv0[0])
          - (uv3[0] - uv2[0])*(uv1[1] - uv0[1]);

    ratio_a = (uv3[0] - uv2[0])*(uv0[1] - uv2[1])
            - (uv3[1] - uv2[1])*(uv0[0] - uv2[0]);
    ratio_b = (uv1[0] - uv0[0])*(uv0[1] - uv2[1])
            - (uv1[1] - uv0[1])*(uv0[0] - uv2[0]);

    if ( ABS( 0.0 - denom ) > tol ) {
      ratio_a /= denom;
      ratio_b /= denom;
    }else{
      if ( ( ABS( 0.0 - ratio_a ) < tol ) &&
	   ( ABS( 0.0 - ratio_b ) < tol ) ) {
	return TRUE;
      }
    }

    if (FALSE) { 
      printf("\ndenom %12.9f a %12.9f b %12.9f\n",
	     denom,ratio_a,ratio_b);
      printf("node0 %d node1 %d node2 %d node3 %d\n",
	     node0, node1, node2, node3);
    }

    meet00 = ( ( ABS( 0.0 - ratio_a ) < tol ) &&
	       ( ABS( 0.0 - ratio_b ) < tol ) );

    meet01 = ( ( ABS( 0.0 - ratio_a ) < tol ) &&
	       ( ABS( 1.0 - ratio_b ) < tol ) );

    meet10 = ( ( ABS( 1.0 - ratio_a ) < tol ) &&
	       ( ABS( 0.0 - ratio_b ) < tol ) );

    meet11 = ( ( ABS( 1.0 - ratio_a ) < tol ) &&
	       ( ABS( 1.0 - ratio_b ) < tol ) );

    if ( ( meet00 && (node0 != node2) ) ||
	 ( meet01 && (node0 != node3) ) ||
	 ( meet10 && (node1 != node2) ) ||
	 ( meet11 && (node1 != node3) ) ) {
      return TRUE; /* segment ends meet without matching node ids */
    }
    if ( !(meet00 || meet01 || meet10 || meet11) ) {
      if ( ( 1.0 >= ratio_a && 0.0 <= ratio_a ) &&
	   ( 1.0 >= ratio_b && 0.0 <= ratio_b ) ) {
	return TRUE; /* intersection */
      }
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

      /* originals are for debugging */
      if ( 0 == ring->originals ) {
	int original;
	ring->originals = ringSegments(ring);
	ring->original_nodes = (int    *)malloc(ring->originals*2*sizeof(int));
	ring->original_uvs   = (double *)malloc(ring->originals*4*sizeof(double));
	for ( original = 0 ; original < (2*ring->originals) ; original++ ) {
	  ring->original_nodes[original] = ring->segment_nodes[original];
	}
	for ( original = 0 ; original < (4*ring->originals) ; original++ ) {
	  ring->original_uvs[original] = ring->segment_uvs[original];
	}
      }

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

static 
double ringTriangleArea(Ring *ring, double *uv0,  double *uv1,  double *uv2)
{
  double edge0[2], edge1[2];

  edge0[0] = uv1[0]-uv0[0];
  edge0[1] = uv1[1]-uv0[1];

  edge1[0] = uv2[0]-uv0[0];
  edge1[1] = uv2[1]-uv0[1];

  return 0.5 * ( edge0[0]*edge1[1] - edge0[1]*edge1[0] );
}

GridBool ringSurroundsTriangle( Ring *ring,
				int node0, int node1, int node2, double *uv2 )
{
  int segment;
  double uv0[2], uv1[2];
  double area;

  for ( segment = 0 ; segment < ringSegments(ring) ; segment++ ) {
    if ( ring->segment_nodes[0+2*segment] == node0 &&
	 ring->segment_nodes[1+2*segment] == node1 ) {
      uv0[0] = ring->segment_uvs[0+4*segment];
      uv0[1] = ring->segment_uvs[1+4*segment];
      uv1[0] = ring->segment_uvs[2+4*segment];
      uv1[1] = ring->segment_uvs[3+4*segment];
      area = ringTriangleArea(ring,uv0,uv1,uv2);
      if (area > 1.0e-14 ) {
	return TRUE;
      } else {
	return FALSE;
      }
    }
  }

  return FALSE;
}
