
#include "ruby.h"
#include "octree.h"

#define GET_OCTREE_FROM_SELF Octree *octree; Data_Get_Struct( self, Octree, octree );

static void octree_free( void *octree )
{
  octreeFree( octree );
}

VALUE octree_new( VALUE class, VALUE boundingBox )
{
  Octree *octree;
  VALUE obj;
  octree = octreeCreate( NUM2DBL(rb_ary_entry(boundingBox,0)),
			 NUM2DBL(rb_ary_entry(boundingBox,1)),
			 NUM2DBL(rb_ary_entry(boundingBox,2)),
			 NUM2DBL(rb_ary_entry(boundingBox,3)),
			 NUM2DBL(rb_ary_entry(boundingBox,4)),
			 NUM2DBL(rb_ary_entry(boundingBox,5)) );

  obj = Data_Wrap_Struct( class, 0, octree_free, octree );
  return obj;
}

VALUE octree_boundingBox( VALUE self )
{
  VALUE rb_boundingBox;
  int i;
  double boundingBox[6];
  GET_OCTREE_FROM_SELF;
  octreeBoundingBox( octree, boundingBox );
  rb_boundingBox = rb_ary_new2(6);
  for(i=0;i<6;i++)
    rb_ary_store(rb_boundingBox, i, rb_float_new(boundingBox[i]) );
  return rb_boundingBox;
}

VALUE octree_nOctant( VALUE self )
{
  GET_OCTREE_FROM_SELF;
  return INT2NUM( octreeNOctant(octree) );
}

VALUE octree_insert( VALUE self, VALUE rb_location, VALUE rb_data )
{
  int i;
  double location[3], data[9];
  GET_OCTREE_FROM_SELF;
  for(i=0;i<3;i++) location[i] = NUM2DBL(rb_ary_entry(rb_location,i));
  for(i=0;i<9;i++) data[i] = NUM2DBL(rb_ary_entry(rb_data,i));
  return (octree == octreeInsert(octree,location,data)?self:Qnil);
}

VALUE octree_query( VALUE self, VALUE rb_location )
{
  int i;
  double location[3], data[9];
  VALUE rb_data;
  GET_OCTREE_FROM_SELF;
  for(i=0;i<3;i++) location[i] = NUM2DBL(rb_ary_entry(rb_location,i));
  if (octree != octreeQuery(octree,location,data) ) return Qnil;
  rb_data = rb_ary_new2(9);
  for(i=0;i<9;i++) rb_ary_store(rb_data,i,rb_float_new(data[i]));
  return rb_data;
}

VALUE cOctree;

void Init_Octree() 
{
  cOctree = rb_define_class( "Octree", rb_cObject );
  rb_define_singleton_method( cOctree, "new", octree_new, 1 );
  rb_define_method( cOctree, "boundingBox", octree_boundingBox, 0 );
  rb_define_method( cOctree, "nOctant", octree_nOctant, 0 );
  rb_define_method( cOctree, "insert", octree_insert, 2 );
  rb_define_method( cOctree, "query", octree_query, 1 );
}
