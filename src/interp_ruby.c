
#include "ruby.h"
#include "interp.h"

#define GET_INTERP_FROM_SELF Interp *interp; Data_Get_Struct( self, Interp, interp );

static void interp_free( void *interp )
{
  interpFree( interp );
}

static VALUE interp_new( VALUE class, VALUE function_id, VALUE order )
{
  Interp *interp;
  VALUE obj;
  interp = interpCreate( NUM2INT(function_id), NUM2INT(order) );
  obj = Data_Wrap_Struct( class, 0, interp_free, interp );
  return obj;
}

VALUE interp_functionId( VALUE self )
{
  GET_INTERP_FROM_SELF;
  return INT2NUM( interpFunctionId(interp) );
}

VALUE interp_order( VALUE self )
{
  GET_INTERP_FROM_SELF;
  return INT2NUM( interpOrder(interp) );
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

static VALUE interp_error( VALUE self, VALUE rb_xyz0, VALUE rb_xyz1, VALUE rb_xyz2, VALUE rb_xyz3 )
{
  int i;
  double xyz0[3], xyz1[3], xyz2[3], xyz3[3];
  double error;
  GET_INTERP_FROM_SELF;
  for (i=0;i<3;i++) xyz0[i] = NUM2DBL(rb_ary_entry(rb_xyz0,i));
  for (i=0;i<3;i++) xyz1[i] = NUM2DBL(rb_ary_entry(rb_xyz1,i));
  for (i=0;i<3;i++) xyz2[i] = NUM2DBL(rb_ary_entry(rb_xyz2,i));
  for (i=0;i<3;i++) xyz3[i] = NUM2DBL(rb_ary_entry(rb_xyz3,i));
  if ( interpError( interp, xyz0, &error ) ) {
    return rb_float_new(error);
  }else{
    return Qnil;
  }
}

VALUE cInterp;

void Init_Interp() 
{
  cInterp = rb_define_class( "Interp", rb_cObject );
  rb_define_singleton_method( cInterp, "new", interp_new, 2 );
  rb_define_method( cInterp, "functionId", interp_functionId, 0 );
  rb_define_method( cInterp, "order", interp_order, 0 );
  rb_define_method( cInterp, "function", interp_function, 1 );
  rb_define_method( cInterp, "error", interp_error, 4 );
}
