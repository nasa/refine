
#include "ruby.h"
#include "gridinsert.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_splitEdge( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return (gridSplitEdge( grid, NUM2INT(n0),  NUM2INT(n1) )==grid?self:Qnil);
}

VALUE cGridInsert;

void Init_GridInsert() 
{
  cGridInsert = rb_define_module( "GridInsert" );
  rb_define_method( cGridInsert, "splitEdge", grid_splitEdge, 2 );
}
