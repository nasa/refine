
#include "ruby.h"
#include "grid.h"

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
  Grid *grid;
  Data_Get_Struct( self, Grid, grid );
  return INT2NUM( gridNNode(grid) );
}

VALUE grid_ncell( VALUE self )
{
  Grid *grid;
  Data_Get_Struct( self, Grid, grid );
  return INT2NUM( gridNCell(grid) );
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
}
