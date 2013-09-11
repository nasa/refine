
#include <stdlib.h>
#include <stdio.h>

#include "ref_recover.h"
#include "ref_malloc.h"

REF_STATUS ref_recover_create( REF_RECOVER *ref_recover_ptr )
{
  REF_RECOVER ref_recover;

  ref_malloc( *ref_recover_ptr, 1, REF_RECOVER_STRUCT );
  ref_recover = (*ref_recover_ptr);

  ref_recover_n(ref_recover) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_recover_free( REF_RECOVER ref_recover )
{
  if ( NULL == (void *)ref_recover ) return REF_NULL;
  ref_free( ref_recover );
  return REF_SUCCESS;
}

