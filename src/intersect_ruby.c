
#include "ruby.h"
#include "intersect.h"

static VALUE intersect_above( VALUE self,
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
  return (intersectAbove( t0, t1, t2, n )?Qtrue:Qfalse);
}

static VALUE intersect_triangleNode( VALUE self, 
				     VALUE tri0, VALUE tri1, VALUE tri2, 
				     VALUE node )
{
  int i;
  double t0[3], t1[3], t2[3], n[3];
  for (i=0;i<3;i++){
    t0[i] = NUM2DBL(rb_ary_entry(tri0, i));
    t1[i] = NUM2DBL(rb_ary_entry(tri1, i));
    t2[i] = NUM2DBL(rb_ary_entry(tri2, i));
    n[i]  = NUM2DBL(rb_ary_entry(node, i));
  }
  return (intersectTriangleNode(t0,t1,t2,n)?Qtrue:Qfalse);
}

static VALUE intersect_triangleSegment( VALUE self, 
				       VALUE tri0, VALUE tri1, VALUE tri2, 
				       VALUE node0, VALUE node1 )
{
  int i;
  double t0[3], t1[3], t2[3], n0[3], n1[3];
  double ratio;
  for (i=0;i<3;i++){
    t0[i] = NUM2DBL(rb_ary_entry(tri0, i));
    t1[i] = NUM2DBL(rb_ary_entry(tri1, i));
    t2[i] = NUM2DBL(rb_ary_entry(tri2, i));
    n0[i] = NUM2DBL(rb_ary_entry(node0, i));
    n1[i] = NUM2DBL(rb_ary_entry(node1, i));
  }
  return (intersectTriangleSegment(t0,t1,t2,n0,n1,&ratio)?Qtrue:Qfalse);
}

static VALUE intersect_insideTet( VALUE self, 
				  VALUE vert0, VALUE vert1, 
				  VALUE vert2, VALUE vert3, 
				  VALUE node )
{
  int i;
  double v0[3], v1[3], v2[3], v3[3], n[3];
  for (i=0;i<3;i++){
    v0[i] = NUM2DBL(rb_ary_entry(vert0, i));
    v1[i] = NUM2DBL(rb_ary_entry(vert1, i));
    v2[i] = NUM2DBL(rb_ary_entry(vert2, i));
    v3[i] = NUM2DBL(rb_ary_entry(vert3, i));
    n[i]  = NUM2DBL(rb_ary_entry(node, i));
  }
  return (intersectInsideTet(v0,v1,v2,v3,n)?Qtrue:Qfalse);
}

static VALUE intersect_tetSegment( VALUE self, 
				   VALUE vert0, VALUE vert1, 
				   VALUE vert2, VALUE vert3, 
				   VALUE node0, VALUE node1 )
{
  int i;
  double v0[3], v1[3], v2[3], v3[3], n0[3], n1[3];
  for (i=0;i<3;i++){
    v0[i] = NUM2DBL(rb_ary_entry(vert0, i));
    v1[i] = NUM2DBL(rb_ary_entry(vert1, i));
    v2[i] = NUM2DBL(rb_ary_entry(vert2, i));
    v3[i] = NUM2DBL(rb_ary_entry(vert3, i));
    n0[i] = NUM2DBL(rb_ary_entry(node0, i));
    n1[i] = NUM2DBL(rb_ary_entry(node1, i));
  }
  return (intersectTetSegment(v0,v1,v2,v3,n0,n1)?Qtrue:Qfalse);
}

static VALUE intersect_tetTet( VALUE self, 
			       VALUE vert0, VALUE vert1, 
			       VALUE vert2, VALUE vert3, 
			       VALUE node0, VALUE node1,
			       VALUE node2, VALUE node3 )
{
  int i;
  double v0[3], v1[3], v2[3], v3[3];
  double n0[3], n1[3], n2[3], n3[3];
  for (i=0;i<3;i++){
    v0[i] = NUM2DBL(rb_ary_entry(vert0, i));
    v1[i] = NUM2DBL(rb_ary_entry(vert1, i));
    v2[i] = NUM2DBL(rb_ary_entry(vert2, i));
    v3[i] = NUM2DBL(rb_ary_entry(vert3, i));
    n0[i] = NUM2DBL(rb_ary_entry(node0, i));
    n1[i] = NUM2DBL(rb_ary_entry(node1, i));
    n2[i] = NUM2DBL(rb_ary_entry(node2, i));
    n3[i] = NUM2DBL(rb_ary_entry(node3, i));
  }
  return (intersectTetTet(v0,v1,v2,v3,n0,n1,n2,n3)?Qtrue:Qfalse);
}

VALUE cIntersect;

void Init_Intersect(  )
{
  cIntersect = rb_define_class( "Intersect", rb_cObject );
  rb_define_method( cIntersect, "above", intersect_above, 4 );
  rb_define_method( cIntersect, "triangleNode", intersect_triangleNode, 4 );
  rb_define_method( cIntersect, "triangleSegment", intersect_triangleSegment, 5 );
  rb_define_method( cIntersect, "insideTet", intersect_insideTet, 5 );
  rb_define_method( cIntersect, "tetSegment", intersect_tetSegment, 6 );
  rb_define_method( cIntersect, "tetTet", intersect_tetTet, 8 );
}
