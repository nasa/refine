
#include "ruby.h"
#include "gridmath.h"

static VALUE grid_subtractVector( VALUE self, VALUE rb_v1, VALUE rb_v2)
{
  int i;
  double v1[3], v2[3], result[3];
  VALUE rb_result;
  for (i=0;i<3;i++){
    v1[i] = NUM2DBL(rb_ary_entry(rb_v1, i));
    v2[i] = NUM2DBL(rb_ary_entry(rb_v2, i));
  }
  gridSubtractVector(v1,v2,result);
  rb_result = rb_ary_new2(3);
  for (i=0;i<3;i++) rb_ary_store(rb_result, i, rb_float_new(result[i]));
  return rb_result;
}

static VALUE grid_dotProduct( VALUE self, VALUE rb_v1, VALUE rb_v2)
{
  int i;
  double v1[3], v2[3];
  for (i=0;i<3;i++){
    v1[i] = NUM2DBL(rb_ary_entry(rb_v1, i));
    v2[i] = NUM2DBL(rb_ary_entry(rb_v2, i));
  }
  return rb_float_new(gridDotProduct(v1,v2));
}

static VALUE grid_crossProduct( VALUE self, VALUE rb_v1, VALUE rb_v2)
{
  int i;
  double v1[3], v2[3], result[3];
  VALUE rb_result;
  for (i=0;i<3;i++){
    v1[i] = NUM2DBL(rb_ary_entry(rb_v1, i));
    v2[i] = NUM2DBL(rb_ary_entry(rb_v2, i));
  }
  gridCrossProduct(v1,v2,result);
  rb_result = rb_ary_new2(3);
  for (i=0;i<3;i++) rb_ary_store(rb_result, i, rb_float_new(result[i]));
  return rb_result;
}

static VALUE grid_vectorLength( VALUE self, VALUE rb_vect )
{
  int i;
  double vect[3];
  for (i=0;i<3;i++) vect[i] = NUM2DBL(rb_ary_entry(rb_vect, i));
  return rb_float_new(gridVectorLength(vect));
}

static VALUE grid_vectorNormalize( VALUE self, VALUE rb_vect )
{
  int i;
  double vect[3];
  for (i=0;i<3;i++) vect[i] = NUM2DBL(rb_ary_entry(rb_vect, i));
  gridVectorNormalize(vect);
  for (i=0;i<3;i++) rb_ary_store(rb_vect, i, rb_float_new(vect[i]));
  return rb_vect;
}

static VALUE grid_vectorOrthogonalize( VALUE self, VALUE rb_vect, VALUE rb_axle)
{
  int i;
  double vect[3], axle[3];
  for (i=0;i<3;i++) vect[i] = NUM2DBL(rb_ary_entry(rb_vect, i));
  for (i=0;i<3;i++) axle[i] = NUM2DBL(rb_ary_entry(rb_axle, i));
  gridVectorOrthogonalize(vect,axle);
  for (i=0;i<3;i++) rb_ary_store(rb_vect, i, rb_float_new(vect[i]));
  return rb_vect;
}

static VALUE grid_barycentricCoordinate( VALUE self, 
					VALUE rb_xyz0, VALUE rb_xyz1, 
					VALUE rb_xyz2, VALUE rb_xyz3, 
					VALUE rb_target )
{
  int i;
  double xyz0[3], xyz1[3], xyz2[3], xyz3[3];
  double target[3], bary[4];
  VALUE rb_bary;
  for (i=0;i<3;i++){
    xyz0[i] = NUM2DBL(rb_ary_entry(rb_xyz0, i));
    xyz1[i] = NUM2DBL(rb_ary_entry(rb_xyz1, i));
    xyz2[i] = NUM2DBL(rb_ary_entry(rb_xyz2, i));
    xyz3[i] = NUM2DBL(rb_ary_entry(rb_xyz3, i));
    target[i] = NUM2DBL(rb_ary_entry(rb_target, i));
  }
  gridBarycentricCoordinate( xyz0, xyz1,xyz2, xyz3, target, bary );
  rb_bary = rb_ary_new2(4);
  for(i=0;i<4;i++) rb_ary_store( rb_bary, i, rb_float_new(bary[i]) );
  return rb_bary;
}

static VALUE grid_barycentricCoordinateTri( VALUE self,
					    VALUE rb_xyz0,
					    VALUE rb_xyz1,
					    VALUE rb_xyz2,
					    VALUE rb_target )
{
  int i;
  double xyz0[3], xyz1[3], xyz2[3], target[3], bary[3];
  VALUE rb_bary;
  for (i=0;i<3;i++) xyz0[i] = NUM2DBL(rb_ary_entry(rb_xyz0, i));
  for (i=0;i<3;i++) xyz1[i] = NUM2DBL(rb_ary_entry(rb_xyz1, i));
  for (i=0;i<3;i++) xyz2[i] = NUM2DBL(rb_ary_entry(rb_xyz2, i));
  for (i=0;i<3;i++) target[i]  = NUM2DBL(rb_ary_entry(rb_target,  i));
  gridBarycentricCoordinateTri( xyz0, xyz1, xyz2, target, bary );
  rb_bary = rb_ary_new2(3);
  for (i=0;i<3;i++) rb_ary_store(rb_bary, i, rb_float_new(bary[i]));
  return rb_bary;
}

static VALUE grid_rotateDirection( VALUE self, VALUE rb_v0, VALUE rb_v1,
				VALUE rb_axle, VALUE rotation )
{
  int i;
  double v0[3], v1[3], axle[3], vect[3];
  VALUE rb_vect;
  for (i=0;i<3;i++) v0[i] = NUM2DBL(rb_ary_entry(rb_v0, i));
  for (i=0;i<3;i++) v1[i] = NUM2DBL(rb_ary_entry(rb_v1, i));
  for (i=0;i<3;i++) axle[i] = NUM2DBL(rb_ary_entry(rb_axle, i));
  gridRotateDirection(v0,v1,axle,NUM2DBL(rotation),vect);
  rb_vect = rb_ary_new2(3);
  for (i=0;i<3;i++) rb_ary_store(rb_vect, i, rb_float_new(vect[i]));
  return rb_vect;
}

static VALUE grid_triDiag( VALUE self, VALUE rb_m )
{
  int i;
  double m[6], d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_d;
  for (i=0;i<6;i++) m[i] = NUM2DBL(rb_ary_entry(rb_m,i));
  gridTriDiag3x3(m,d,e,q0,q1,q2);
  rb_d = rb_ary_new2(3);
  for(i=0;i<3;i++) rb_ary_store( rb_d, i, rb_float_new(d[i]) );
  return rb_d;
}

static VALUE grid_triOffDiag( VALUE self, VALUE rb_m )
{
  int i;
  double m[6], d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_e;
  for (i=0;i<6;i++) m[i] = NUM2DBL(rb_ary_entry(rb_m,i));
  gridTriDiag3x3(m,d,e,q0,q1,q2);
  rb_e = rb_ary_new2(3);
  for(i=0;i<3;i++) rb_ary_store( rb_e, i, rb_float_new(e[i]) );
  return rb_e;
}

static VALUE grid_triDiagTransform0( VALUE self, VALUE rb_m )
{
  int i;
  double m[6], d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_q;
  for (i=0;i<6;i++) m[i] = NUM2DBL(rb_ary_entry(rb_m,i));
  gridTriDiag3x3(m,d,e,q0,q1,q2);
  rb_q = rb_ary_new2(3);
  for(i=0;i<3;i++) rb_ary_store( rb_q, i, rb_float_new(q0[i]) );
  return rb_q;
}

static VALUE grid_triDiagTransform1( VALUE self, VALUE rb_m )
{
  int i;
  double m[6], d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_q;
  for (i=0;i<6;i++) m[i] = NUM2DBL(rb_ary_entry(rb_m,i));
  gridTriDiag3x3(m,d,e,q0,q1,q2);
  rb_q = rb_ary_new2(3);
  for(i=0;i<3;i++) rb_ary_store( rb_q, i, rb_float_new(q1[i]) );
  return rb_q;
}

static VALUE grid_triDiagTransform2( VALUE self, VALUE rb_m )
{
  int i;
  double m[6], d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_q;
  for (i=0;i<6;i++) m[i] = NUM2DBL(rb_ary_entry(rb_m,i));
  gridTriDiag3x3(m,d,e,q0,q1,q2);
  rb_q = rb_ary_new2(3);
  for(i=0;i<3;i++) rb_ary_store( rb_q, i, rb_float_new(q2[i]) );
  return rb_q;
}

static VALUE grid_eigTriDiag( VALUE self, VALUE rb_d, VALUE rb_e )
{
  int i;
  double d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_eig;
  for (i=0;i<3;i++) d[i] = NUM2DBL(rb_ary_entry(rb_d,i));
  for (i=0;i<3;i++) e[i] = NUM2DBL(rb_ary_entry(rb_e,i));
  q0[0] = 1.0; q1[0] = 0.0; q2[0] = 0.0;
  q0[1] = 0.0; q1[1] = 1.0; q2[1] = 0.0;
  q0[2] = 0.0; q1[2] = 0.0; q2[2] = 1.0;
  if ( gridEigTriDiag3x3( d, e, q0, q1, q2 ) ) {
    rb_eig = rb_ary_new2(3);
    for(i=0;i<3;i++) rb_ary_store( rb_eig, i, rb_float_new(d[i]) );
    return rb_eig;
  } else {
    return Qnil;
  }
}

static VALUE grid_vectTriDiag0( VALUE self, VALUE rb_d, VALUE rb_e, 
				VALUE rb_q0, VALUE rb_q1,  VALUE rb_q2 )
{
  int i;
  double d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_vect;
  for (i=0;i<3;i++) d[i] = NUM2DBL(rb_ary_entry(rb_d,i));
  for (i=0;i<3;i++) e[i] = NUM2DBL(rb_ary_entry(rb_e,i));
  for (i=0;i<3;i++) q0[i] = NUM2DBL(rb_ary_entry(rb_q0,i));
  for (i=0;i<3;i++) q1[i] = NUM2DBL(rb_ary_entry(rb_q1,i));
  for (i=0;i<3;i++) q2[i] = NUM2DBL(rb_ary_entry(rb_q2,i));
  if ( gridEigTriDiag3x3( d, e, q0, q1, q2 ) ) {
    rb_vect = rb_ary_new2(3);
    for(i=0;i<3;i++) rb_ary_store( rb_vect, i, rb_float_new(q0[i]) );
    return rb_vect;
  } else {
    return Qnil;
  }
}

static VALUE grid_vectTriDiag1( VALUE self, VALUE rb_d, VALUE rb_e, 
				VALUE rb_q0, VALUE rb_q1,  VALUE rb_q2 )
{
  int i;
  double d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_vect;
  for (i=0;i<3;i++) d[i] = NUM2DBL(rb_ary_entry(rb_d,i));
  for (i=0;i<3;i++) e[i] = NUM2DBL(rb_ary_entry(rb_e,i));
  for (i=0;i<3;i++) q0[i] = NUM2DBL(rb_ary_entry(rb_q0,i));
  for (i=0;i<3;i++) q1[i] = NUM2DBL(rb_ary_entry(rb_q1,i));
  for (i=0;i<3;i++) q2[i] = NUM2DBL(rb_ary_entry(rb_q2,i));
  if ( gridEigTriDiag3x3( d, e, q0, q1, q2 ) ) {
    rb_vect = rb_ary_new2(3);
    for(i=0;i<3;i++) rb_ary_store( rb_vect, i, rb_float_new(q1[i]) );
    return rb_vect;
  } else {
    return Qnil;
  }
}

static VALUE grid_vectTriDiag2( VALUE self, VALUE rb_d, VALUE rb_e, 
				VALUE rb_q0, VALUE rb_q1,  VALUE rb_q2 )
{
  int i;
  double d[3], e[3], q0[3], q1[3], q2[3];
  VALUE rb_vect;
  for (i=0;i<3;i++) d[i] = NUM2DBL(rb_ary_entry(rb_d,i));
  for (i=0;i<3;i++) e[i] = NUM2DBL(rb_ary_entry(rb_e,i));
  for (i=0;i<3;i++) q0[i] = NUM2DBL(rb_ary_entry(rb_q0,i));
  for (i=0;i<3;i++) q1[i] = NUM2DBL(rb_ary_entry(rb_q1,i));
  for (i=0;i<3;i++) q2[i] = NUM2DBL(rb_ary_entry(rb_q2,i));
  if ( gridEigTriDiag3x3( d, e, q0, q1, q2 ) ) {
    rb_vect = rb_ary_new2(3);
    for(i=0;i<3;i++) rb_ary_store( rb_vect, i, rb_float_new(q2[i]) );
    return rb_vect;
  } else {
    return Qnil;
  }
}

static VALUE grid_lu3x3( VALUE self, VALUE rb_a )
{
  int i;
  double a[9], lu[9];
  VALUE rb_lu;
  for (i=0;i<9;i++) a[i] = NUM2DBL(rb_ary_entry(rb_a,i));
  gridLU3x3( a, lu );
  rb_lu = rb_ary_new2(9);
  for(i=0;i<9;i++) rb_ary_store( rb_lu, i, rb_float_new(lu[i]) );
  return rb_lu;
}

static VALUE grid_backsolve3x3( VALUE self, VALUE rb_lu, VALUE rb_b )
{
  int i;
  double lu[9], b[3];
  VALUE rb_s;
  for (i=0;i<9;i++) lu[i] = NUM2DBL(rb_ary_entry(rb_lu,i));
  for (i=0;i<3;i++)  b[i] = NUM2DBL(rb_ary_entry(rb_b,i));
  gridBackSolve3x3( lu, b );
  rb_s = rb_ary_new2(3);
  for(i=0;i<3;i++) rb_ary_store( rb_s, i, rb_float_new(b[i]) );
  return rb_s;
}

static VALUE grid_matrixDeterminate( VALUE self, VALUE rb_m )
{
  int i;
  double m[9];
  for (i=0;i<9;i++) m[i] = NUM2DBL(rb_ary_entry(rb_m,i));
  return rb_float_new(gridMatrixDeterminate( m ));
}

static VALUE grid_gaussianElimination( VALUE self, 
				       VALUE rb_m, VALUE rb_n, VALUE rb_a )
{
  int i, m, n;
  double *a;
  m = NUM2INT(rb_m);
  n = NUM2INT(rb_n);
  a = (double *)malloc( m*n*sizeof(double) );
  for (i=0;i<m*n;i++) a[i] = NUM2DBL(rb_ary_entry(rb_a,i));
  if ( !gridGaussianElimination( m, n, a ) ) {
    free(a);
    return Qnil;
  }
  for (i=0;i<m*n;i++) rb_ary_store( rb_a, i, rb_float_new(a[i]) );
  free(a);
  return rb_a;
}

static VALUE grid_gaussianBacksolve( VALUE self, 
				       VALUE rb_m, VALUE rb_n, VALUE rb_a )
{
  int i, m, n;
  double *a;
  m = NUM2INT(rb_m);
  n = NUM2INT(rb_n);
  a = (double *)malloc( m*n*sizeof(double) );
  for (i=0;i<m*n;i++) a[i] = NUM2DBL(rb_ary_entry(rb_a,i));
  if ( !gridGaussianBacksolve( m, n, a ) ) {
    free(a);
    return Qnil;
  }
  for (i=0;i<m*n;i++) rb_ary_store( rb_a, i, rb_float_new(a[i]) );
  free(a);
  return rb_a;
}

VALUE cGridMath;

void Init_GridMath(  )
{
  cGridMath = rb_define_class( "GridMath", rb_cObject );
  rb_define_method( cGridMath, "subtractVector", grid_subtractVector, 2 );
  rb_define_method( cGridMath, "dotProduct", grid_dotProduct, 2 );
  rb_define_method( cGridMath, "crossProduct", grid_crossProduct, 2 );
  rb_define_method( cGridMath, "vectorLength", grid_vectorLength, 1 );
  rb_define_method( cGridMath, "vectorNormalize", grid_vectorNormalize, 1 );
  rb_define_method( cGridMath, "vectorOrthogonalize", 
                                grid_vectorOrthogonalize, 2 );
  rb_define_method( cGridMath, "barycentricCoordinate", 
		                grid_barycentricCoordinate, 5 );
  rb_define_method( cGridMath, "barycentricCoordinateTri",
		    grid_barycentricCoordinateTri, 4);
  rb_define_method( cGridMath, "rotateDirection", grid_rotateDirection, 4 );
  rb_define_method( cGridMath, "triDiag", grid_triDiag, 1 );
  rb_define_method( cGridMath, "triOffDiag", grid_triOffDiag, 1 );
  rb_define_method( cGridMath, "triDiagTransform0", grid_triDiagTransform0, 1 );
  rb_define_method( cGridMath, "triDiagTransform1", grid_triDiagTransform1, 1 );
  rb_define_method( cGridMath, "triDiagTransform2", grid_triDiagTransform2, 1 );
  rb_define_method( cGridMath, "eigTriDiag", grid_eigTriDiag, 2 );
  rb_define_method( cGridMath, "vectTriDiag0", grid_vectTriDiag0, 5 );
  rb_define_method( cGridMath, "vectTriDiag1", grid_vectTriDiag1, 5 );
  rb_define_method( cGridMath, "vectTriDiag2", grid_vectTriDiag2, 5 );
  rb_define_method( cGridMath, "lu3x3", grid_lu3x3, 1 );
  rb_define_method( cGridMath, "backsolve3x3", grid_backsolve3x3, 2 );
  rb_define_method( cGridMath, "matrixDeterminate", grid_matrixDeterminate, 1);
  rb_define_method( cGridMath, "gaussianElimination", grid_gaussianElimination, 3 );
  rb_define_method( cGridMath, "gaussianBacksolve", grid_gaussianBacksolve, 3 );
}
