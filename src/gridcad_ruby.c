
#include "ruby.h"
#include "gridcad.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_projectNodeToEdge( VALUE self, VALUE node, VALUE edgeId )
{
  GET_GRID_FROM_SELF;
  return (gridProjectNodeToEdge( grid, NUM2INT(node), NUM2INT(edgeId) )==grid?self:Qnil);
}

VALUE grid_safeProject( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridSafeProject( grid, NUM2INT(node) )==grid?self:Qnil);
}

VALUE cGridCAD;

void Init_GridCAD() 
{
  cGridCAD = rb_define_module( "GridCAD" );
  rb_define_method( cGridCAD, "projectNodeToEdge", grid_projectNodeToEdge, 2 );
  rb_define_method( cGridCAD, "safeProject", grid_safeProject, 1 );
}
