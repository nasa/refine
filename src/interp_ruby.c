
#include "ruby.h"
#include "interp.h"

#define GET_INTERP_FROM_SELF Interp *interp; Data_Get_Struct( self, Interp, interp );

static void interp_free( void *interp )
{
  interpFree( interp );
}

static VALUE interp_new( VALUE class, VALUE function_id )
{
  Interp *interp;
  VALUE obj;
  interp = interpCreate( NUM2INT(function_id) );
  obj = Data_Wrap_Struct( class, 0, interp_free, interp );
  return obj;
}

static VALUE interp_function( VALUE self, VALUE rb_xyz )
{
  int i;
  double xyz[3];
  double func;
  GET_INTERP_FROM_SELF;
  for (i=0;i<3;i++) xyz[i] = NUM2DBL(rb_ary_entry(rb_xyz,i));
  if ( interpFunction( interp, xyz, &func ) ) {
    return rb_float_new(func);
  }else{
    return Qnil;
  }
}

VALUE cInterp;

void Init_Interp() 
{
  cInterp = rb_define_class( "Interp", rb_cObject );
  rb_define_singleton_method( cInterp, "new", interp_new, 1 );
  rb_define_method( cInterp, "function", interp_function, 1 );
}
