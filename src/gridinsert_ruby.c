
#include "ruby.h"
#include "gridinsert.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_adapt( VALUE self, VALUE minLength, VALUE maxLength)
{
  GET_GRID_FROM_SELF;
  return (gridAdapt( grid, 
		     NUM2DBL(minLength), NUM2DBL(maxLength) )==grid?self:Qnil);
}

VALUE grid_splitEdge( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return INT2NUM(gridSplitEdge( grid, NUM2INT(n0),  NUM2INT(n1) ));
}

VALUE grid_collapseEdge( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return (gridCollapseEdge( grid, NUM2INT(n0),  NUM2INT(n1) )==grid?self:Qnil);
}
VALUE cGridInsert;

void Init_GridInsert() 
{
  cGridInsert = rb_define_module( "GridInsert" );
  rb_define_method( cGridInsert, "adapt", grid_adapt, 2 );
  rb_define_method( cGridInsert, "splitEdge", grid_splitEdge, 2 );
  rb_define_method( cGridInsert, "collapseEdge", grid_collapseEdge, 2 );
}
