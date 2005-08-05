
#include "ruby.h"
#include "ring.h"

#define GET_RING_FROM_SELF Ring *ring; Data_Get_Struct( self, Ring, ring );

static void ring_free( void *ring )
{
  ringFree( ring );
}

VALUE ring_new( VALUE class )
{
  Ring *ring;
  VALUE obj;
  ring = ringCreate(  );
  obj = Data_Wrap_Struct( class, 0, ring_free, ring );
  return obj;
}

VALUE ring_segments( VALUE self )
{
  GET_RING_FROM_SELF;
  return INT2NUM( ringSegments(ring) );
}

VALUE ring_addSegment( VALUE self,
		       VALUE n0, VALUE n1, 
		       VALUE rb_uv0, VALUE rb_uv1 )
{
  double uv0[2], uv1[2];
  GET_RING_FROM_SELF;

  uv0[0] = NUM2DBL(rb_ary_entry(rb_uv0,0));
  uv0[1] = NUM2DBL(rb_ary_entry(rb_uv0,1));

  uv1[0] = NUM2DBL(rb_ary_entry(rb_uv1,0));
  uv1[1] = NUM2DBL(rb_ary_entry(rb_uv1,1));

  return (ringAddSegment(ring,NUM2INT(n0),NUM2INT(n1),uv0,uv1)==ring?self:Qnil);
}

VALUE ring_segment( VALUE self, VALUE segment )
{
  int node0, node1;
  double uv0[2], uv1[2];
  VALUE rb;
  VALUE rb_uv0, rb_uv1;
  GET_RING_FROM_SELF;

  if ( ring != ringSegment(ring, NUM2INT(segment),
			   &node0, &node1, uv0, uv1 ) ) return Qnil;

  rb = rb_ary_new2(4);

  rb_ary_store( rb, 0, INT2NUM(node0) );
  rb_ary_store( rb, 1, INT2NUM(node1) );

  rb_uv0 = rb_ary_new2(2);
  rb_ary_store( rb_uv0, 0, rb_float_new(uv0[0]) );
  rb_ary_store( rb_uv0, 1, rb_float_new(uv0[1]) );
  rb_ary_store( rb, 2, rb_uv0 );

  rb_uv1 = rb_ary_new2(2);
  rb_ary_store( rb_uv1, 0, rb_float_new(uv1[0]) );
  rb_ary_store( rb_uv1, 1, rb_float_new(uv1[1]) );
  rb_ary_store( rb, 3, rb_uv1 );

  return rb;
}

VALUE ring_triangles( VALUE self )
{
  GET_RING_FROM_SELF;
  return INT2NUM( ringTriangles(ring) );
}

VALUE ring_addTriangle( VALUE self,
			VALUE n0, VALUE n1, VALUE n2, 
			VALUE rb_uv2 )
{
  double uv2[2];
  GET_RING_FROM_SELF;

  uv2[0] = NUM2DBL(rb_ary_entry(rb_uv2,0));
  uv2[1] = NUM2DBL(rb_ary_entry(rb_uv2,1));

  return (ringAddTriangle(ring,
			  NUM2INT(n0),NUM2INT(n1),NUM2INT(n2),
			  uv2)==ring?self:Qnil);
}

VALUE cRing;

void Init_Ring() 
{
  cRing = rb_define_class( "Ring", rb_cObject );
  rb_define_singleton_method( cRing, "new", ring_new, 0 );
  rb_define_method( cRing, "segments", ring_segments, 0 );
  rb_define_method( cRing, "addSegment", ring_addSegment, 4 );
  rb_define_method( cRing, "segment", ring_segment, 1 );
  rb_define_method( cRing, "triangles", ring_triangles, 0 );
  rb_define_method( cRing, "addTriangle", ring_addTriangle, 4 );
}
