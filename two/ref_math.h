
#ifndef REF_MATH_H
#define REF_MATH_H

#include <math.h>

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define ref_math_dot(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
#define ref_math_pi (3.14159265358979)
#define ref_math_in_degrees(radians) ((radians)*180.0/3.14159265358979)
#define ref_math_in_radians(degrees) ((degrees)/180.0*3.14159265358979)

/* tester : printf(" %e / %e = %e \n",n,d,n/d) */
#define ref_math_divisible(n,d) (ABS(1.0e20*d)>ABS(n))

#define ref_math_cross_product(v0,v1,product)	   \
(product)[0] = (v0)[1]*(v1)[2] - (v0)[2]*(v1)[1]; \
(product)[1] = (v0)[2]*(v1)[0] - (v0)[0]*(v1)[2]; \
(product)[2] = (v0)[0]*(v1)[1] - (v0)[1]*(v1)[0]; 

REF_STATUS ref_math_normalize( REF_DBL *normal );

END_C_DECLORATION

#endif /* REF_MATH_H */
