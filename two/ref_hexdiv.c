
#include <stdlib.h>
#include <stdio.h>

#include "ref_hexdiv.h"

REF_STATUS ref_hexdiv_create( REF_HEXDIV *ref_hexdiv_ptr, REF_GRID ref_grid )
{
  REF_HEXDIV ref_hexdiv;
  REF_INT face;

  (*ref_hexdiv_ptr) = NULL;
  (*ref_hexdiv_ptr) = (REF_HEXDIV)malloc( sizeof(REF_HEXDIV_STRUCT) );
  RNS(*ref_hexdiv_ptr,"malloc ref_hexdiv NULL");

  ref_hexdiv = *ref_hexdiv_ptr;

  ref_hexdiv_grid(ref_hexdiv) = ref_grid;

  RSS( ref_face_create( &(ref_hexdiv_face( ref_hexdiv )), 
			ref_hexdiv_grid(ref_hexdiv) ), "create face" );

  ref_hexdiv->mark = (REF_INT *)malloc( ref_face_n(ref_hexdiv_face(ref_hexdiv)) 
					* sizeof(REF_INT));
  RNS(ref_hexdiv->mark,"malloc mark NULL");

  for ( face=0 ; face < ref_face_n(ref_hexdiv_face(ref_hexdiv)) ; face++ )
    ref_hexdiv_mark( ref_hexdiv, face ) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_free( REF_HEXDIV ref_hexdiv )
{
  if ( NULL == (void *)ref_hexdiv ) return REF_NULL;

  free( ref_hexdiv->mark );
  RSS( ref_face_free( ref_hexdiv_face( ref_hexdiv ) ), "free face" );

  ref_cond_free( ref_hexdiv );

  return REF_SUCCESS;
}
