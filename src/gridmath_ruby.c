
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
  VALUE rb_result;
  for (i=0;i<3;i++) vect[i] = NUM2DBL(rb_ary_entry(rb_vect, i));
  gridVectorNormalize(vect);
  for (i=0;i<3;i++) rb_ary_store(rb_vect, i, rb_float_new(vect[i]));
  return rb_vect;
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
}
