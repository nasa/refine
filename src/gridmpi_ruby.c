
#include "ruby.h"
#include "gridmpi.h"

#define GET_GRID_FROM_SELF Grid *grid; Data_Get_Struct( self, Grid, grid );

VALUE grid_identityNodeGlobal( VALUE self, VALUE offset )
{
  GET_GRID_FROM_SELF;
  return( grid == gridIdentityNodeGlobal( grid, NUM2INT(offset) )?self:Qnil);
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

VALUE grid_parallelEdgeCollapse( VALUE self, VALUE rb_queue, 
				 VALUE node0, VALUE node1 )
{
  Queue *queue;
  GET_GRID_FROM_SELF;
  Data_Get_Struct( rb_queue, Queue, queue );
  return ( grid==gridParallelEdgeCollapse( grid, queue, 
					   NUM2INT(node0),
					   NUM2INT(node1) )?self:Qnil );
}

VALUE grid_parallelEdgeSwap( VALUE self, VALUE rb_queue, 
			     VALUE node0, VALUE node1 )
{
  Queue *queue;
  GET_GRID_FROM_SELF;
  Data_Get_Struct( rb_queue, Queue, queue );
  return ( grid==gridParallelEdgeSwap( grid, queue, 
				       NUM2INT(node0), 
				       NUM2INT(node1) )?self:Qnil );
}

VALUE grid_applyQueue( VALUE self, VALUE rb_queue )
{
  Queue *queue;
  GET_GRID_FROM_SELF;
  Data_Get_Struct( rb_queue, Queue, queue );
  return ( grid==gridApplyQueue( grid, queue )?self:Qnil );
}


VALUE grid_ghostDataCountByPartition( VALUE self, VALUE total_number_of_partitions )
{
  VALUE array;
  int i, *partition_nodes;
  GET_GRID_FROM_SELF;

  partition_nodes = malloc(NUM2INT(total_number_of_partitions) * sizeof(int));

  if (grid != gridGhostDataCountByPartition(grid,
					    NUM2INT(total_number_of_partitions),
					    partition_nodes)) {
    free(partition_nodes);
    return Qnil;
  }

  array = rb_ary_new2(NUM2INT(total_number_of_partitions));
  for (i=0;i<NUM2INT(total_number_of_partitions);i++) 
    rb_ary_store(array,i,INT2NUM(partition_nodes[i]));
  return array;
}

VALUE cGridMPI;

void Init_GridMPI() 
{
  cGridMPI = rb_define_module( "GridMPI" );
  rb_define_method( cGridMPI, "setAllLocal", grid_setAllLocal, 0 );
  rb_define_method( cGridMPI, "identityNodeGlobal", grid_identityNodeGlobal, 1);
  rb_define_method( cGridMPI, "setGhost", grid_setGhost, 1 );
  rb_define_method( cGridMPI, "parallelEdgeSplit", grid_parallelEdgeSplit, 3 );
  rb_define_method( cGridMPI, "parallelEdgeCollapse", grid_parallelEdgeCollapse, 3 );
  rb_define_method( cGridMPI, "parallelEdgeSwap", grid_parallelEdgeSwap, 3 );
  rb_define_method( cGridMPI, "applyQueue", grid_applyQueue, 1 );

  rb_define_method( cGridMPI, "ghostDataCountByPartition", grid_ghostDataCountByPartition, 1 );
}
