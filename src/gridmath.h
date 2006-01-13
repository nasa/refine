/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDMATH_H
#define GRIDMATH_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

#define gridPI (3.14159265358979)
#define gridConvertRadianToDegree(radian) ((radian)*57.2957795130823)
#define gridConvertDegreeToRadian(degree) ((degree)*0.0174532925199433)

#define gridSubtractVector(v1,v2,result) \
(result)[0] = (v1)[0] - (v2)[0]; \
(result)[1] = (v1)[1] - (v2)[1]; \
(result)[2] = (v1)[2] - (v2)[2];

#define gridAverageVector(v1,v2,result) \
(result)[0] = 0.5*((v1)[0] + (v2)[0]); \
(result)[1] = 0.5*((v1)[1] + (v2)[1]); \
(result)[2] = 0.5*((v1)[2] + (v2)[2]);

#define gridAddToVector(v,add) \
(v)[0] += (add)[0]; \
(v)[1] += (add)[1]; \
(v)[2] += (add)[2];

#define gridDotProduct(v1, v2) \
( (v1)[0]*(v2)[0] +  (v1)[1]*(v2)[1] + (v1)[2]*(v2)[2] )

#define gridCrossProduct(edge1,edge2,norm) \
(norm)[0] = (edge1)[1]*(edge2)[2] - (edge1)[2]*(edge2)[1]; \
(norm)[1] = (edge1)[2]*(edge2)[0] - (edge1)[0]*(edge2)[2]; \
(norm)[2] = (edge1)[0]*(edge2)[1] - (edge1)[1]*(edge2)[0]; 

#define gridVectorCopy(a,b) (a)[0]=(b)[0];(a)[1]=(b)[1];(a)[2]=(b)[2];
#define gridVectorMirror(a,b) (a)[0]=-(b)[0];(a)[1]=-(b)[1];(a)[2]=-(b)[2];
#define gridVectorScale(v,s) (v)[0]*=(s); (v)[1]*=(s); (v)[2]*=(s);

double gridVectorLength(double *norm);
void gridVectorNormalize(double *norm);
void gridVectorOrthogonalize(double *norm, double *axle);

void gridProjectToTriangle(double *projected_target, 
			   double *xyz0, double *xyz1, double *xyz2 );

void gridRotateDirection(double *v0, double *v1, 
			 double *axle, double rotation, double *result);

void gridTriDiag3x3( double *m, double *d, double *e, 
		     double *q0, double *q1, double *q2 );

#define gridSign(a,b) (b>=0?ABS(a):-ABS(a))
GridBool gridEigTriDiag3x3( double *d, double *e,
		        double *q0, double *q1, double *q2 );

void gridEigSort3x3( double *eigenValues, double *v0, double *v1, double *v2 );
void gridEigOrtho3x3( double *v0, double *v1, double *v2 );

void gridLU3x3( double *a, double *lu );
void gridBackSolve3x3( double *lu, double *b );

#define gridMatrixDeterminate(m) ( \
(m)[0]*(m)[4]*(m)[8] + \
(m)[1]*(m)[5]*(m)[6] + \
(m)[2]*(m)[3]*(m)[7] - \
(m)[0]*(m)[5]*(m)[7] - \
(m)[1]*(m)[3]*(m)[8] - \
(m)[2]*(m)[4]*(m)[6] ) 

void gridTriangularBarycentricCoordinate3D( double *xyz0, double *xyz1,
					    double *xyz2, double *xyz, 
					    double *bary );
END_C_DECLORATION

#endif /* GRIDMATH_H */
