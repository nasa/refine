
#include "ruby.h"
#include "intersect.h"

static VALUE intersect_side( VALUE self, 
			     VALUE tri0, VALUE tri1, VALUE tri2, VALUE node )
{
  int i;
  double t0[3], t1[3], t2[3], n[3];
  for (i=0;i<3;i++){
    t0[i] = NUM2DBL(rb_ary_entry(tri0, i));
    t1[i] = NUM2DBL(rb_ary_entry(tri1, i));
    t2[i] = NUM2DBL(rb_ary_entry(tri2, i));
    n[i]  = NUM2DBL(rb_ary_entry(node, i));
  }
  return INT2FIX( intersectSide( t0, t1, t2, n ) );
}

VALUE cIntersect;

void Init_Intersect(  )
{
  cIntersect = rb_define_class( "Intersect", rb_cObject );
  rb_define_method( cIntersect, "side", intersect_side, 4 );
}
