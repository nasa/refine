
#include "ruby.h"
#include "gridshape.h"
 
#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_shapeJacobian1( VALUE self, 
			   VALUE rb_n0, VALUE rb_n1, VALUE rb_n2, VALUE rb_n3,
			   VALUE rb_w )
{
  int i;
  double n0[3], n1[3], n2[3], n3[3], w[3], j[9];
  VALUE rb_j;
  GET_GRID_FROM_SELF;
  for (i=0;i<3;i++) n0[i] = NUM2DBL(rb_ary_entry(rb_n0,i));
  for (i=0;i<3;i++) n1[i] = NUM2DBL(rb_ary_entry(rb_n1,i));
  for (i=0;i<3;i++) n2[i] = NUM2DBL(rb_ary_entry(rb_n2,i));
  for (i=0;i<3;i++) n3[i] = NUM2DBL(rb_ary_entry(rb_n3,i));
  for (i=0;i<3;i++)  w[i] = NUM2DBL(rb_ary_entry(rb_w,i));
  if (grid!=gridShapeJacobian1(grid,n0,n1,n2,n3,w,j)) return Qnil;
  rb_j = rb_ary_new2(9);
  for(i=0;i<9;i++) rb_ary_store( rb_j, i, rb_float_new(j[i]) );
  return rb_j;
}

VALUE grid_shapeJacobian2( VALUE self, 
			   VALUE rb_n0, VALUE rb_n1, VALUE rb_n2, VALUE rb_n3,
			   VALUE rb_e01, VALUE rb_e02, VALUE rb_e03,
			   VALUE rb_e12, VALUE rb_e13, VALUE rb_e23,
			   VALUE rb_w )
{
  int i;
  double n0[3], n1[3], n2[3], n3[3];
  double e01[3], e02[3], e03[3];
  double e12[3], e13[3], e23[3];
  double w[3], j[9];
  VALUE rb_j;
  GET_GRID_FROM_SELF;
  for (i=0;i<3;i++) n0[i] = NUM2DBL(rb_ary_entry(rb_n0,i));
  for (i=0;i<3;i++) n1[i] = NUM2DBL(rb_ary_entry(rb_n1,i));
  for (i=0;i<3;i++) n2[i] = NUM2DBL(rb_ary_entry(rb_n2,i));
  for (i=0;i<3;i++) n3[i] = NUM2DBL(rb_ary_entry(rb_n3,i));

  for (i=0;i<3;i++) e01[i] = NUM2DBL(rb_ary_entry(rb_e01,i));
  for (i=0;i<3;i++) e02[i] = NUM2DBL(rb_ary_entry(rb_e02,i));
  for (i=0;i<3;i++) e03[i] = NUM2DBL(rb_ary_entry(rb_e03,i));

  for (i=0;i<3;i++) e12[i] = NUM2DBL(rb_ary_entry(rb_e12,i));
  for (i=0;i<3;i++) e13[i] = NUM2DBL(rb_ary_entry(rb_e13,i));
  for (i=0;i<3;i++) e23[i] = NUM2DBL(rb_ary_entry(rb_e23,i));

  for (i=0;i<3;i++)  w[i] = NUM2DBL(rb_ary_entry(rb_w,i));
  if (grid!=gridShapeJacobian2(grid,n0,n1,n2,n3, 
			       e01, e02, e03, 
			       e12, e13, e23, 
			       w,j)) return Qnil;
  rb_j = rb_ary_new2(9);
  for(i=0;i<9;i++) rb_ary_store( rb_j, i, rb_float_new(j[i]) );
  return rb_j;
}

VALUE cGridShape;
 
void Init_GridShape()
{
   cGridShape = rb_define_module( "GridShape" );
   rb_define_method( cGridShape, "shapeJacobian1", grid_shapeJacobian1, 5 );
   rb_define_method( cGridShape, "shapeJacobian2", grid_shapeJacobian1, 11 );
 }

