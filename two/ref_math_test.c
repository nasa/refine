#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_math.h"

#include "ref_test.h"

int main( void )
{

  { /* same 1 */
    REF_DBL vect1[3]={1.0,0.0,0.0};
    REF_DBL vect2[3]={1.0,0.0,0.0};
    TWDS(1.0,ref_math_dot(vect1,vect2),-1.0,"dot");
  }

  { /* orth 1 */
    REF_DBL vect1[3]={1.0,0.0,0.0};
    REF_DBL vect2[3]={0.0,1.0,0.0};
    TWDS(0.0,ref_math_dot(vect1,vect2),-1.0,"dot");
  }

  { /* two */
    REF_DBL vect1[3]={1.0,0.0,0.0};
    REF_DBL vect2[3]={2.0,1.0,0.0};
    TWDS(2.0,ref_math_dot(vect1,vect2),-1.0,"dot");
  }

  { /* normal zero */
    REF_DBL vect[3]={0.0,0.0,0.0};
    TFS(ref_math_normalize( vect ), "expect fail");
    TWDS(0.0,vect[0],-1.0,"same");
    TWDS(0.0,vect[1],-1.0,"same");
    TWDS(0.0,vect[2],-1.0,"same");
  }

  { /* normal one */
    REF_DBL vect[3]={1.0,0.0,0.0};
    TSS(ref_math_normalize( vect ), "expect success");
    TWDS(1.0,vect[0],-1.0,"same");
    TWDS(0.0,vect[1],-1.0,"same");
    TWDS(0.0,vect[2],-1.0,"same");
  }

  { /* normal two */
    REF_DBL vect[3]={2.0,0.0,0.0};
    TSS(ref_math_normalize( vect ), "expect success");
    TWDS(1.0,vect[0],-1.0,"same");
    TWDS(0.0,vect[1],-1.0,"same");
    TWDS(0.0,vect[2],-1.0,"same");
  }

  return 0;
}
