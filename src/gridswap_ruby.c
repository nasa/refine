
#include "ruby.h"
#include "gridswap.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_swapEdge( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return (gridSwapEdge( grid, NUM2INT(n0),  NUM2INT(n1) )==grid?self:Qnil);
}

VALUE grid_swap( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridSwap( grid )==grid?self:Qnil);
}

VALUE cGridSwap;

void Init_GridSwap() 
{
  cGridSwap = rb_define_module( "GridSwap" );
  rb_define_method( cGridSwap, "swapEdge", grid_swapEdge, 2 );
  rb_define_method( cGridSwap, "swap", grid_swap, 0 );
}
