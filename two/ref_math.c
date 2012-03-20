
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_math.h"

REF_STATUS ref_math_normalize( REF_DBL *normal )
{
  REF_DBL length;

  length = ref_math_dot(normal,normal);

  length = sqrt(length);

  if ( !ref_math_divisible(normal[0],length) ||
       !ref_math_divisible(normal[1],length) ||
       !ref_math_divisible(normal[2],length) ) return REF_DIV_ZERO;
  
  normal[0] /= length;
  normal[1] /= length;
  normal[2] /= length;

  length = ref_math_dot(normal,normal);
  RAS( (ABS(length-1.0) < 1.0e-13), "vector length not unity"); 

  return REF_SUCCESS;
}
