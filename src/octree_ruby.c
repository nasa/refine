
#include "ruby.h"
#include "octree.h"

#define GET_OCTREE_FROM_SELF Octree *octree; Data_Get_Struct( self, Octree, octree );

static void octree_free( void *octree )
{
  octreeFree( octree );
}

VALUE octree_new( VALUE class, 
		  VALUE xMin, VALUE xMax,
		  VALUE yMin, VALUE yMax,
		  VALUE zMin, VALUE zMax )
{
  Octree *octree;
  VALUE obj;
  octree = octreeCreate( NUM2DBL(xMin), NUM2DBL(xMax), 
			 NUM2DBL(yMin), NUM2DBL(yMax), 
			 NUM2DBL(zMin), NUM2DBL(zMax) );
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


VALUE cOctree;

void Init_Octree() 
{
  cOctree = rb_define_class( "Octree", rb_cObject );
  rb_define_singleton_method( cOctree, "new", octree_new, 6 );
  rb_define_method( cOctree, "boundingBox", octree_boundingBox, 0 );
  rb_define_method( cOctree, "nOctant", octree_nOctant, 0 );
}
