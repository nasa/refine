
#include <stdlib.h>
#include <stdio.h>

#include "ref_face.h"

REF_STATUS ref_face_create( REF_FACE *ref_face_ptr, REF_GRID ref_grid )
{
  REF_FACE ref_face;

  (*ref_face_ptr) = NULL;
  (*ref_face_ptr) = (REF_FACE)malloc( sizeof(REF_FACE_STRUCT) );
  RNS(*ref_face_ptr,"malloc ref_face NULL");

  ref_face = *ref_face_ptr;

  ref_face_n(ref_face) = 0;
  ref_face_max(ref_face) = 10;

  ref_face->f2n = (REF_INT *)malloc( 4 * ref_face_max(ref_face) * 
				     sizeof(REF_INT) );
  RNS(ref_face->f2n,"malloc ref_face->f2n NULL");

  RSS( ref_adj_create( &(ref_face_adj( ref_face )) ), "create adj" );

  return REF_SUCCESS;
}

REF_STATUS ref_face_free( REF_FACE ref_face )
{
  if ( NULL == (void *)ref_face ) return REF_NULL;

  RSS( ref_adj_free( ref_face_adj( ref_face ) ), "free adj" );
  free( ref_face->f2n );

  ref_cond_free( ref_face );

  return REF_SUCCESS;
}
