
#include "ruby.h"
#include "gridswap.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_swapFace( VALUE self, VALUE n0, VALUE n1, VALUE n2)
{
  GET_GRID_FROM_SELF;
  return (gridSwapFace( grid, NULL, NUM2INT(n0),  
			NUM2INT(n1), NUM2INT(n2) )==grid?self:Qnil);
}

VALUE grid_swapEdge( VALUE self, VALUE n0, VALUE n1 )
{
  GET_GRID_FROM_SELF;
  return (gridSwapEdge( grid, NULL, NUM2INT(n0), NUM2INT(n1) )==grid?self:Qnil);
}

VALUE grid_swap( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridSwap( grid )==grid?self:Qnil);
}

VALUE grid_removeTwoFaceCell( VALUE self, VALUE cell )
{
  GET_GRID_FROM_SELF;
  return (gridRemoveTwoFaceCell( grid, NULL, NUM2INT(cell) )==grid?self:Qnil);
}

VALUE cGridSwap;

void Init_GridSwap() 
{
  cGridSwap = rb_define_module( "GridSwap" );
  rb_define_method( cGridSwap, "swapFace", grid_swapFace, 3 );
  rb_define_method( cGridSwap, "swapEdge", grid_swapEdge, 2 );
  rb_define_method( cGridSwap, "swap", grid_swap, 0 );
  rb_define_method( cGridSwap, "removeTwoFaceCell", grid_removeTwoFaceCell, 1 );
}
