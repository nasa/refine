
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

VALUE grid_dump( VALUE self )
{
  GET_GRID_FROM_SELF;
  gridDump( grid );
  return self;
}

VALUE grid_nnode( VALUE self )
{
  GET_GRID_FROM_SELF;
  return INT2NUM( gridNNode(grid) );
}

VALUE grid_maxcell( VALUE self )
{
  GET_GRID_FROM_SELF;
  return INT2NUM( gridMaxCell(grid) );
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
  Grid *returnedGrid;
  GET_GRID_FROM_SELF;
  returnedGrid = 
    gridRegisterNodeCell( grid, NUM2INT(nodeId), NUM2INT(cellId) );
  return (returnedGrid==NULL?Qnil:self);
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

VALUE grid_currentNodeCell( VALUE self )
{
  GET_GRID_FROM_SELF;
  return INT2NUM( gridCurrentNodeCell(grid) );
}

VALUE grid_moreNodeCell( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridMoreNodeCell(grid)?Qtrue:Qfalse);
}

VALUE grid_nextNodeCell( VALUE self )
{
  GET_GRID_FROM_SELF;
  gridNextNodeCell(grid);
  return self;
}

VALUE grid_removeNodeCell( VALUE self, VALUE nodeId, VALUE cellId )
{
  Grid *returnedGrid;
  GET_GRID_FROM_SELF;
  returnedGrid = gridRemoveNodeCell( grid, NUM2INT(nodeId), NUM2INT(cellId) );
  return (returnedGrid==NULL?Qnil:self);
}

VALUE grid_cellExists( VALUE self, VALUE nodeId, VALUE cellId )
{
  GET_GRID_FROM_SELF;
  return
    ( gridCellExists( grid, NUM2INT(nodeId), NUM2INT(cellId) )?Qtrue:Qfalse );
}

VALUE grid_addCell( VALUE self, VALUE n0, VALUE n1, VALUE n2, VALUE n3 )
{
  Grid *returnedGrid;
  GET_GRID_FROM_SELF;
  returnedGrid = 
    gridAddCell( grid, NUM2INT(n0),  NUM2INT(n1),  NUM2INT(n2),  NUM2INT(n3) );
  return (returnedGrid==NULL?Qnil:self);
}

VALUE grid_pack( VALUE self )
{
  GET_GRID_FROM_SELF;
  gridPack( grid );
  return self;
}

VALUE grid_gem( VALUE self, VALUE n0, VALUE n1 )
{
  VALUE rb_gem;
  int i, ngem;
  GET_GRID_FROM_SELF;
  gridMakeGem( grid, NUM2INT(n0), NUM2INT(n1) );
  ngem = gridNGem(grid);
  rb_gem = rb_ary_new2(ngem);
  for ( i=0 ; i < ngem ; i++ ) 
    rb_ary_store( rb_gem, i, INT2NUM(gridGem(grid,i)) );
  return rb_gem;
}

VALUE grid_equator( VALUE self, VALUE n0, VALUE n1 )
{
#define MAXEQU 200
  int i, nequ;
  int equ[MAXEQU];
  VALUE rb_equ; // bug, what is return if this is cleaned up?
  GET_GRID_FROM_SELF;
  gridEquator( grid, NUM2INT(n0), NUM2INT(n1), MAXEQU, &nequ, equ );
  rb_equ = rb_ary_new();
  for ( i=0 ; i < nequ ; i++ ) rb_equ = rb_ary_push( rb_equ, INT2NUM(equ[i]) );
  return rb_equ;
}

VALUE grid_orient( VALUE self, VALUE c0, VALUE c1 , VALUE c2, VALUE c3,
		   VALUE n0, VALUE n1 )
{
  VALUE rb_result;
  int i, cell[4], result[4];
  GET_GRID_FROM_SELF;

  cell[0] = NUM2INT(c0);
  cell[1] = NUM2INT(c1);
  cell[2] = NUM2INT(c2);
  cell[3] = NUM2INT(c3);

  result[0] = NUM2INT(n0);
  result[1] = NUM2INT(n1);

  if ( NULL == gridOrient( grid, cell, result ) ) return Qnil;

  rb_result = rb_ary_new2(4);
  for ( i=0 ; i<4 ; i++ ) rb_ary_store( rb_result, i, INT2NUM(result[i]) );
  return rb_result;
}

VALUE cGrid;

void Init_Grid() 
{
  cGrid = rb_define_class( "Grid", rb_cObject );
  rb_define_singleton_method( cGrid, "new", grid_new, 3 );
  rb_define_method( cGrid, "initialize", grid_init, 0 );
  rb_define_method( cGrid, "dump", grid_dump, 0 );
  rb_define_method( cGrid, "nnode", grid_nnode, 0 );
  rb_define_method( cGrid, "maxcell", grid_maxcell, 0 );
  rb_define_method( cGrid, "ncell", grid_ncell, 0 );
  rb_define_method( cGrid, "nodeDeg", grid_nodeDeg, 1 );
  rb_define_method( cGrid, "pack", grid_pack, 0 );
  rb_define_method( cGrid, "registerNodeCell", grid_registerNodeCell, 2 );
  rb_define_method( cGrid, "validNodeCell", grid_validNodeCell, 0 );
  rb_define_method( cGrid, "firstNodeCell", grid_firstNodeCell, 1 );
  rb_define_method( cGrid, "currentNodeCell", grid_currentNodeCell, 0 );
  rb_define_method( cGrid, "moreNodeCell", grid_moreNodeCell, 0 );
  rb_define_method( cGrid, "nextNodeCell", grid_nextNodeCell, 0 );
  rb_define_method( cGrid, "removeNodeCell", grid_removeNodeCell, 2 );
  rb_define_method( cGrid, "cellExists", grid_cellExists, 2 );
  rb_define_method( cGrid, "addCell", grid_addCell, 4 );
  rb_define_method( cGrid, "gem", grid_gem, 2 );
  rb_define_method( cGrid, "equator", grid_equator, 2 );
  rb_define_method( cGrid, "orient", grid_orient, 6 );
}
