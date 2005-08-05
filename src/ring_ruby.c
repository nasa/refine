
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

VALUE ring_triangles( VALUE self )
{
  GET_RING_FROM_SELF;
  return INT2NUM( ringTriangles(ring) );
}

VALUE cRing;

void Init_Ring() 
{
  cRing = rb_define_class( "Ring", rb_cObject );
  rb_define_singleton_method( cRing, "new", ring_new, 0 );
  rb_define_method( cRing, "segments", ring_segments, 0 );
  rb_define_method( cRing, "triangles", ring_triangles, 0 );
}
