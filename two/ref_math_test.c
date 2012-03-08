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
    RWDS(1.0,ref_math_dot(vect1,vect2),-1.0,"dot");
  }

  { /* orth 1 */
    REF_DBL vect1[3]={1.0,0.0,0.0};
    REF_DBL vect2[3]={0.0,1.0,0.0};
    RWDS(0.0,ref_math_dot(vect1,vect2),-1.0,"dot");
  }

  { /* two */
    REF_DBL vect1[3]={1.0,0.0,0.0};
    REF_DBL vect2[3]={2.0,1.0,0.0};
    RWDS(2.0,ref_math_dot(vect1,vect2),-1.0,"dot");
  }

  { /* normal zero */
    REF_DBL vect[3]={0.0,0.0,0.0};
    REIS(REF_DIV_ZERO,ref_math_normalize( vect ), "expect fail");
    RWDS(0.0,vect[0],-1.0,"same");
    RWDS(0.0,vect[1],-1.0,"same");
    RWDS(0.0,vect[2],-1.0,"same");
  }

  { /* normal one */
    REF_DBL vect[3]={1.0,0.0,0.0};
    RSS(ref_math_normalize( vect ), "expect success");
    RWDS(1.0,vect[0],-1.0,"same");
    RWDS(0.0,vect[1],-1.0,"same");
    RWDS(0.0,vect[2],-1.0,"same");
  }

  { /* normal two */
    REF_DBL vect[3]={2.0,0.0,0.0};
    RSS(ref_math_normalize( vect ), "expect success");
    RWDS(1.0,vect[0],-1.0,"same");
    RWDS(0.0,vect[1],-1.0,"same");
    RWDS(0.0,vect[2],-1.0,"same");
  }

  { /* acos */
    RWDS(0.0,acos( 1.0),-1.0,"acos");
    RWDS(1.04719755119660,acos( 0.5),-1.0,"acos");
    RWDS(1.57079632679490,acos( 0.0),-1.0,"acos");
    RWDS(2.09439510239320,acos(-0.5),-1.0,"acos");
    RWDS(3.14159265358979,acos(-1.0),-1.0,"acos");
  }

  return 0;
}
