
#include "ruby.h"
#include "gridmpi.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_identityGlobal( VALUE self, VALUE offset )
{
  GET_GRID_FROM_SELF;
  return( grid == gridIdentityGlobal( grid, NUM2INT(offset) )?self:Qnil);
}

VALUE grid_setAllLocal( VALUE self )
{
  GET_GRID_FROM_SELF;
  return (gridSetAllLocal(grid)==grid?self:Qnil);
}

VALUE grid_setGhost( VALUE self, VALUE node )
{
  GET_GRID_FROM_SELF;
  return (gridSetGhost(grid, NUM2INT(node))==grid?self:Qnil);
}

VALUE grid_parallelEdgeSplit( VALUE self, VALUE rb_queue, 
			      VALUE node0, VALUE node1 )
{
  Queue *queue;
  GET_GRID_FROM_SELF;
  Data_Get_Struct( rb_queue, Queue, queue );
  return INT2NUM( gridParallelEdgeSplit( grid, queue, 
					 NUM2INT(node0), NUM2INT(node1) ) );
}

VALUE cGridMPI;

void Init_GridMPI() 
{
  cGridMPI = rb_define_module( "GridMPI" );
  rb_define_method( cGridMPI, "setAllLocal", grid_setAllLocal, 0 );
  rb_define_method( cGridMPI, "identityGlobal", grid_identityGlobal, 1 );
  rb_define_method( cGridMPI, "setGhost", grid_setGhost, 1 );
  rb_define_method( cGridMPI, "parallelEdgeSplit", grid_parallelEdgeSplit, 3 );
}
