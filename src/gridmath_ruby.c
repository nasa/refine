
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
  q0[0] = 0.0; q1[0] = 1.0; q2[0] = 0.0;
  q0[0] = 0.0; q1[0] = 0.0; q2[0] = 1.0;
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

VALUE cGridMath;

void Init_GridMath(  )
{
  cGridMath = rb_define_class( "GridMath", rb_cObject );
  rb_define_method( cGridMath, "subtractVector", grid_subtractVector, 2 );
  rb_define_method( cGridMath, "dotProduct", grid_dotProduct, 2 );
  rb_define_method( cGridMath, "crossProduct", grid_crossProduct, 2 );
  rb_define_method( cGridMath, "vectorLength", grid_vectorLength, 1 );
  rb_define_method( cGridMath, "vectorNormalize", grid_vectorNormalize, 1 );
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
}
