
#include "ruby.h"
#include "gridcad.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_projectNodeToEdge( VALUE self, VALUE node, VALUE edgeId )
{
  GET_GRID_FROM_SELF;
  return (gridProjectNodeToEdge( grid, NUM2INT(node), NUM2INT(edgeId) )==grid?self:Qnil);
}

VALUE grid_projectNodeToFace( VALUE self, VALUE node, VALUE faceId )
{
  GET_GRID_FROM_SELF;
  return (gridProjectNodeToFace( grid, NUM2INT(node), NUM2INT(faceId) )==grid?self:Qnil);
}

VALUE grid_safeProjectNode( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridSafeProjectNode( grid, NUM2INT(node) )==grid?self:Qnil);
}

VALUE grid_smooth( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridSmooth( grid )==grid?self:Qnil);
}

VALUE grid_smoothNode( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridSmoothNode( grid, NUM2INT(node) )==grid?self:Qnil);
}

VALUE grid_optimizeUV( VALUE self, VALUE node, VALUE rb_dudv )
{
  double dudv[2];
  GET_GRID_FROM_SELF;
  dudv[0] = NUM2DBL(rb_ary_entry(rb_dudv,0));
  dudv[1] = NUM2DBL(rb_ary_entry(rb_dudv,1));
  return (gridOptimizeUV( grid, NUM2INT(node), dudv )==grid?self:Qnil);
}

VALUE cGridCAD;

void Init_GridCAD() 
{
  cGridCAD = rb_define_module( "GridCAD" );
  rb_define_method( cGridCAD, "projectNodeToEdge", grid_projectNodeToEdge, 2 );
  rb_define_method( cGridCAD, "projectNodeToFace", grid_projectNodeToFace, 2 );
  rb_define_method( cGridCAD, "safeProjectNode", grid_safeProjectNode, 1 );
  rb_define_method( cGridCAD, "smooth", grid_smooth, 0 );
  rb_define_method( cGridCAD, "smoothNode", grid_smoothNode, 1 );
  rb_define_method( cGridCAD, "optimizeUV", grid_optimizeUV, 2 );
}
