
#include "ruby.h"
#include "gridmpi.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_identityGlobal( VALUE self )
{
  GET_GRID_FROM_SELF;
  return( grid == gridIdentityGlobal( grid )?self:Qnil);
}

VALUE grid_setAllLocal( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridSetAllLocal(grid)==grid?self:Qnil);
}

VALUE cGridMPI;

void Init_GridMPI() 
{
  cGridMPI = rb_define_module( "GridMPI" );
  rb_define_method( cGridMPI, "setAllLocal", grid_setAllLocal, 0 );
  rb_define_method( cGridMPI, "identityGlobal", grid_identityGlobal, 0 );
}
