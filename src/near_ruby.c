
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

VALUE near_leftIndex( VALUE self )
{
  GET_NEAR_FROM_SELF;
  return INT2NUM(nearLeftIndex(near));
}

VALUE near_rightIndex( VALUE self )
{
  GET_NEAR_FROM_SELF;
  return INT2NUM(nearRightIndex(near));
}

VALUE near_insert( VALUE rb_self, VALUE rb_child )
{
  Near *self, *child; 
  Data_Get_Struct( rb_self,  Near, self );
  Data_Get_Struct( rb_child, Near, child );
  return (self==nearInsert(self,child)?rb_self:Qnil);
}

VALUE near_leftRadius( VALUE self )
{
  GET_NEAR_FROM_SELF;
  return rb_float_new(nearLeftRadius(near));
}

VALUE near_rightRadius( VALUE self )
{
  GET_NEAR_FROM_SELF;
  return rb_float_new(nearRightRadius(near));
}

VALUE near_visualize( VALUE self )
{
  GET_NEAR_FROM_SELF;
  return (near==nearVisualize(near)?self:Qnil);
}

VALUE near_collisions( VALUE self, VALUE rb_target )
{
  Near *tree, *target;
  Data_Get_Struct( self,      Near, tree );
  Data_Get_Struct( rb_target, Near, target );

  return INT2NUM(nearCollisions(tree,target));
}

VALUE near_touched( VALUE self, VALUE rb_target )
{
  Near *tree, *target;
  int i, collisions, found, maxfound, *list;
  VALUE array;
  Data_Get_Struct( self,      Near, tree );
  Data_Get_Struct( rb_target, Near, target );

  collisions = nearCollisions(tree,target);
  list = malloc(MAX(1,collisions)*sizeof(int));

  maxfound = collisions;
  found = 0;
  if (tree != nearTouched(tree,target,&found,maxfound,list)) {
    free(list);
    return Qnil;
  }

  array = rb_ary_new2(collisions);
  for(i=0;i<collisions;i++) rb_ary_store(array, i, INT2NUM(list[i]));

  free(list);
  return array;
}

VALUE near_nearestIndex( VALUE rb_root, VALUE rb_key )
{
  Near *root, *key;
  Data_Get_Struct( rb_root, Near, root );
  Data_Get_Struct( rb_key,  Near, key );

  return INT2NUM(nearNearestIndex(root,key));
}

VALUE near_nearestIndexAndDistance( VALUE rb_root, VALUE rb_key )
{
  Near *root, *key;
  int index;
  double distance;
  VALUE answer;
  Data_Get_Struct( rb_root, Near, root );
  Data_Get_Struct( rb_key,  Near, key );
  
  if (NULL == nearNearestIndexAndDistance(root,key,&index,&distance)) 
    return Qnil;

  answer = rb_ary_new2(2);
  rb_ary_store(answer, 0, INT2NUM(index));
  rb_ary_store(answer, 1, rb_float_new(distance));

  return answer;
}

VALUE cNear;

void Init_Near() 
{
  cNear = rb_define_class( "Near", rb_cObject );
  rb_define_singleton_method( cNear, "new", near_new, 5 );
  rb_define_method( cNear, "index", near_index, 0 );
  rb_define_method( cNear, "clearance", near_clearance, 1 );
  rb_define_method( cNear, "leftIndex", near_leftIndex, 0 );
  rb_define_method( cNear, "rightIndex", near_rightIndex, 0 );
  rb_define_method( cNear, "insert", near_insert, 1 );
  rb_define_method( cNear, "leftRadius", near_leftRadius, 0 );
  rb_define_method( cNear, "rightRadius", near_rightRadius, 0 );
  rb_define_method( cNear, "visualize", near_visualize, 0 );
  rb_define_method( cNear, "collisions", near_collisions, 1 );
  rb_define_method( cNear, "touched", near_touched, 1 );
  rb_define_method( cNear, "nearestIndex", near_nearestIndex, 1 );
  rb_define_method( cNear, "nearestIndexAndDistance", 
		    near_nearestIndexAndDistance, 1 );
}
