/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDMATH_H
#define GRIDMATH_H

#include "master_header.h"

BEGIN_C_DECLORATION

#define gridPI (3.14159265358979)
#define gridConvertRadianToDegree(radian) ((radian)*57.2957795130823)
#define gridConvertDegreeToRadian(degree) ((degree)*0.0174532925199433)

#define gridSubtractVector(v1,v2,result) \
result[0] = v1[0] - v2[0]; result[1] = v1[1] - v2[1]; result[2] = v1[2] - v2[2];

#define gridDotProduct(v1, v2) ( v1[0]*v2[0] +  v1[1]*v2[1] + v1[2]*v2[2] )

#define gridCrossProduct(edge1,edge2,norm) \
norm[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1]; \
norm[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2]; \
norm[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0]; 

#define gridVectorCopy(a,b) a[0]=b[0];a[1]=b[1];a[2]=b[2];

double gridVectorLength(double *norm);
void gridVectorNormalize(double *norm);

void gridTriDiag3x3( double *m, double *d, double *e, 
		     double *q0, double *q1, double *q2 );

#define gridSign(a,b) (b>=0?ABS(a):-ABS(a))
bool gridEigTriDiag3x3( double *d, double *e,
		        double *q0, double *q1, double *q2 );

void gridEigSort3x3( double *eigenValues, double *v0, double *v1, double *v2 );
void gridEigOrtho3x3( double *v0, double *v1, double *v2 );

END_C_DECLORATION

#endif /* GRIDMATH_H */
