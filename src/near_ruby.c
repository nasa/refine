
#include "ruby.h"
#include "near.h"

#define GET_NEAR_FROM_SELF Near *near; Data_Get_Struct( self, Near, near );

static void near_free( void *near )
{
  nearFree( near );
}

VALUE near_new( VALUE class, 
		VALUE index, VALUE x, VALUE y, VALUE z, VALUE radius )
{
  Near *near;
  VALUE obj;
  near = nearCreate( NUM2INT(index), 
		     NUM2DBL(x), NUM2DBL(y), NUM2DBL(z), NUM2DBL(radius) );
  obj = Data_Wrap_Struct( class, 0, near_free, near );
  return obj;
}

VALUE near_index( VALUE self )
{
  GET_NEAR_FROM_SELF;
  return INT2NUM(nearIndex(near));
}

VALUE near_clearance( VALUE rb_self, VALUE rb_other )
{
  Near *self, *other; 
  Data_Get_Struct( rb_self,  Near, self );
  Data_Get_Struct( rb_other, Near, other );
  return(rb_float_new(nearClearance(self, other)));
}

VALUE near_rightIndex( VALUE self )
{
  GET_NEAR_FROM_SELF;
  return INT2NUM(nearRightIndex(near));
}

VALUE near_leftIndex( VALUE self )
{
  GET_NEAR_FROM_SELF;
  return INT2NUM(nearLeftIndex(near));
}

VALUE near_insert( VALUE rb_self, VALUE rb_child )
{
  Near *self, *child; 
  Data_Get_Struct( rb_self,  Near, self );
  Data_Get_Struct( rb_child, Near, child );
  return (self==nearInsert(self,child)?rb_self:Qnil);
}

VALUE cNear;

void Init_Near() 
{
  cNear = rb_define_class( "Near", rb_cObject );
  rb_define_singleton_method( cNear, "new", near_new, 5 );
  rb_define_method( cNear, "index", near_index, 0 );
  rb_define_method( cNear, "clearance", near_clearance, 1 );
  rb_define_method( cNear, "rightIndex", near_rightIndex, 0 );
  rb_define_method( cNear, "leftIndex", near_leftIndex, 0 );
  rb_define_method( cNear, "insert", near_insert, 1 );
}
