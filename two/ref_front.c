
#include <stdlib.h>
#include <stdio.h>

#include "ref_front.h"
#include "ref_malloc.h"

REF_STATUS ref_front_create( REF_FRONT *ref_front_ptr, REF_INT node_per )
{
  REF_FRONT ref_front;
  REF_INT face;

  ref_malloc( *ref_front_ptr, 1, REF_FRONT_STRUCT );
  ref_front = (*ref_front_ptr);

  ref_front_n(ref_front) = 0;
  ref_front_node_per(ref_front) = node_per;

  ref_front_max(ref_front) = 10;

  ref_malloc( ref_front->f2n, ref_front_max(ref_front) *
	      ref_front_node_per(ref_front), REF_INT);
  for ( face = 0 ; face < ref_front_max(ref_front) ; face++ ) 
    {
      ref_front_f2n(ref_front,0,face) = REF_EMPTY;
      ref_front_f2n(ref_front,1,face) = face+1;
    }
  ref_front_f2n(ref_front,1,ref_front_max(ref_front)-1) = REF_EMPTY;
  ref_front_blank(ref_front) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_front_free( REF_FRONT ref_front )
{
  if ( NULL == (void *)ref_front ) return REF_NULL;
  ref_free( ref_front->f2n );
  ref_free( ref_front );
  return REF_SUCCESS;
}

REF_STATUS ref_front_insert( REF_FRONT ref_front, REF_INT *nodes )
{
  REF_INT node, face;
  REF_INT orig, chunk;

  if ( REF_EMPTY == ref_front_blank(ref_front) ) 
    {
      orig = ref_front_max(ref_front);
      chunk = MAX(100,(REF_INT)(1.5*(REF_DBL)orig));
      ref_front_max(ref_front) = orig + chunk;

      ref_realloc( ref_front->f2n, ref_front_node_per(ref_front) *
		   ref_front_max(ref_front), REF_INT );

      for (face=orig;face < ref_front_max(ref_front); face++ ) 
	{
	  ref_front_f2n(ref_front,0,face)= REF_EMPTY; 
	  ref_front_f2n(ref_front,1,face) = face+1; 
	}
      ref_front_f2n(ref_front,1,(ref_front->max)-1) = REF_EMPTY; 
      ref_front_blank(ref_front) = orig;
    }

  face = ref_front_blank(ref_front);
  ref_front_blank(ref_front) = ref_front_f2n(ref_front,1,face);
  for ( node = 0 ; node < ref_front_node_per(ref_front) ; node++ )
    ref_front_f2n(ref_front,node,face) = nodes[node];

  ref_front_n(ref_front)++;

  return REF_SUCCESS;
}

