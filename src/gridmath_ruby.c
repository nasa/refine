
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

VALUE cGridMath;

void Init_GridMath(  )
{
  cGridMath = rb_define_class( "GridMath", rb_cObject );
  rb_define_method( cGridMath, "subtractVector", grid_subtractVector, 2 );
}
