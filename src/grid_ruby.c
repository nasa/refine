
#include "ruby.h"
#include "grid.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

static void grid_free( void *grid )
{
  gridFree( grid );
}

VALUE grid_init( VALUE self ) // not needed but for example
{
  return self;
}

VALUE grid_new( VALUE class, VALUE nnode, VALUE ncell, VALUE nlist )
{
  VALUE *argv;
  Grid *grid;
  VALUE obj;
  grid = gridCreate( NUM2INT(nnode), NUM2INT(ncell), NUM2INT(nlist) );
  obj = Data_Wrap_Struct( class, 0, grid_free, grid );
  rb_obj_call_init( obj, 0, argv ); // not needed but for example
  return obj;
}

VALUE grid_nnode( VALUE self )
{
  GET_GRID_FROM_SELF;
  return INT2NUM( gridNNode(grid) );
}

VALUE grid_ncell( VALUE self )
{
  GET_GRID_FROM_SELF;
  return INT2NUM( gridNCell(grid) );
}

VALUE grid_nodeDeg( VALUE self, VALUE nodeId )
{
  GET_GRID_FROM_SELF;
  return INT2NUM( gridNodeDeg(grid, NUM2INT(nodeId) ) );
}

VALUE grid_registerNodeCell( VALUE self, VALUE nodeId, VALUE cellId )
{
  GET_GRID_FROM_SELF;
  gridRegisterNodeCell(grid, NUM2INT(nodeId), NUM2INT(cellId) );
  return self;
}

VALUE grid_validNodeCell( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridValidNodeCell(grid)?Qtrue:Qfalse);
}

VALUE grid_firstNodeCell( VALUE self, VALUE nodeId )
{
  GET_GRID_FROM_SELF;
  gridFirstNodeCell(grid, NUM2INT(nodeId) );
  return Qnil;
}

VALUE cGrid;

void Init_Grid() 
{
  printf("in Init_Grid\n");
  cGrid = rb_define_class( "Grid", rb_cObject );
  rb_define_singleton_method( cGrid, "new", grid_new, 3 );
  rb_define_method( cGrid, "initialize", grid_init, 0 );
  rb_define_method( cGrid, "nnode", grid_nnode, 0 );
  rb_define_method( cGrid, "ncell", grid_ncell, 0 );
  rb_define_method( cGrid, "nodeDeg", grid_nodeDeg, 1 );
  rb_define_method( cGrid, "registerNodeCell", grid_registerNodeCell, 2 );
  rb_define_method( cGrid, "validNodeCell", grid_validNodeCell, 0 );
  rb_define_method( cGrid, "firstNodeCell", grid_firstNodeCell, 1 );
}
