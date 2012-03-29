
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_fixture.h"
#include "ref_mpi.h"
#include "ref_part.h"

REF_STATUS ref_fixture_tet_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[4] = {0,1,2,3};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS( ref_node_initialize_n_global( ref_node, 4 ), "init glob" );

  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");

  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pyr_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[5] = {0,1,2,3,4};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,4,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS( ref_node_initialize_n_global( ref_node, 5 ), "init glob" );

  RSS(ref_cell_add(ref_grid_pyr(ref_grid),nodes,&cell),"add pyr");

  nodes[0] = 0; nodes[1] = 3; nodes[2] = 4; nodes[3] = 1; nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add qua");

  nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 1; nodes[1] = 4; nodes[2] = 2; nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 4; nodes[1] = 3; nodes[2] = 2; nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 0; nodes[1] = 2; nodes[2] = 3; nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  return REF_SUCCESS;
}

#define add_that_node( node, x, y, z )					\
  RSS(ref_node_add(ref_node,global[(node)],&(local[(node)])),"add node"); \
  ref_node_xyz(ref_node,0,local[(node)]) = (x);				\
  ref_node_xyz(ref_node,1,local[(node)]) = (y);				\
  ref_node_xyz(ref_node,2,local[(node)]) = (z);				\
  ref_node_part(ref_node,local[(node)]) =				\
    ref_part_implicit( nnodesg, ref_mpi_n,				\
		       ref_node_global(ref_node,local[(node)]) );
   

REF_STATUS ref_fixture_pri_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 6;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  
  global[0] = 0; global[1] = 1; global[2] = 2;
  global[3] = 3; global[4] = 4; global[5] = 5;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[4] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[5] ) )
    {
      add_that_node(0,0.0,0.0,0.0);
      add_that_node(1,1.0,0.0,0.0);
      add_that_node(2,0.0,1.0,0.0);
      add_that_node(3,0.0,0.0,1.0);
      add_that_node(4,1.0,0.0,1.0);
      add_that_node(5,0.0,1.0,1.0);

      RSS(ref_cell_add(ref_grid_pri(ref_grid),local,&cell),"add prism");
    }

  RSS( ref_node_initialize_n_global( ref_node, nnodesg ), "init glob" );

  global[0] = 0; global[1] = 3; global[2] = 4; global[3] = 1; global[4] = 10;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      RSS( ref_node_local(ref_node,global[3], &(local[3])),"loc");
      local[4]=global[4];
      RSS(ref_cell_add(ref_grid_qua(ref_grid),local,&cell),"add quad");
    }

  global[0] = 3; global[1] = 5; global[2] = 4; global[3] = 100;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      local[3]=global[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid),local,&cell),"add tri");
    }

  global[0] = 0; global[1] = 1; global[2] = 2; global[3] = 101;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      local[3]=global[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid),local,&cell),"add tri");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pri_tet_cap_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 6;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  
  global[0] = 0; global[1] = 1; global[2] = 2;
  global[3] = 3; global[4] = 4; global[5] = 5;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[4] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[5] ) )
    {
      add_that_node(0,0.0,0.0,0.0);
      add_that_node(1,1.0,0.0,0.0);
      add_that_node(2,0.0,1.0,0.0);
      add_that_node(3,0.0,0.0,1.0);
      add_that_node(4,1.0,0.0,1.0);
      add_that_node(5,0.0,1.0,1.0);

      RSS(ref_cell_add(ref_grid_pri(ref_grid),local,&cell),"add prism");
    }

  global[0] = 3; global[1] = 4; global[2] = 5;
  global[3] = 6;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) )
    {
      add_that_node(0,0.0,0.0,1.0);
      add_that_node(1,1.0,0.0,1.0);
      add_that_node(2,0.0,1.0,1.0);
      add_that_node(3,0.3,0.3,1.0);

      RSS(ref_cell_add(ref_grid_tet(ref_grid),local,&cell),"add prism");
    }

  RSS( ref_node_initialize_n_global( ref_node, nnodesg ), "init glob" );

  global[0] = 0; global[1] = 3; global[2] = 4; global[3] = 1; global[4] = 10;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      RSS( ref_node_local(ref_node,global[3], &(local[3])),"loc");
      local[4]=global[4];
      RSS(ref_cell_add(ref_grid_qua(ref_grid),local,&cell),"add quad");
    }

  global[0] = 3; global[1] = 5; global[2] = 4; global[3] = 100;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      local[3]=global[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid),local,&cell),"add tri");
    }

  global[0] = 0; global[1] = 1; global[2] = 2; global[3] = 101;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      local[3]=global[3];
      RSS(ref_cell_add(ref_grid_tri(ref_grid),local,&cell),"add tri");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pri_stack_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnodesg = 12;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  global[0] = 0; global[1] = 1; global[2] = 2;
  global[3] = 3; global[4] = 4; global[5] = 5;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[4] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[5] ) )
    {
      add_that_node(0,0.0,0.0,0.0);
      add_that_node(1,1.0,0.0,0.0);
      add_that_node(2,0.0,1.0,0.0);
      add_that_node(3,0.0,0.0,1.0);
      add_that_node(4,1.0,0.0,1.0);
      add_that_node(5,0.0,1.0,1.0);

      RSS(ref_cell_add(ref_grid_pri(ref_grid),local,&cell),"add prism");
    }

  global[0] = 3; global[1] = 4; global[2] = 5;
  global[3] = 6; global[4] = 7; global[5] = 8;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[4] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[5] ) )
    {
      add_that_node(0,0.0,0.0,1.0);
      add_that_node(1,1.0,0.0,1.0);
      add_that_node(2,0.0,1.0,1.0);
      add_that_node(3,0.0,0.0,2.0);
      add_that_node(4,1.0,0.0,2.0);
      add_that_node(5,0.0,1.0,2.0);

      RSS(ref_cell_add(ref_grid_pri(ref_grid),local,&cell),"add prism");
    }

  global[0] = 6; global[1] = 7; global[2] = 8;
  global[3] = 9; global[4] =10; global[5] =11;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[4] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[5] ) )
    {
      add_that_node(0,0.0,0.0,2.0);
      add_that_node(1,1.0,0.0,2.0);
      add_that_node(2,0.0,1.0,2.0);
      add_that_node(3,0.0,0.0,3.0);
      add_that_node(4,1.0,0.0,3.0);
      add_that_node(5,0.0,1.0,3.0);

      RSS(ref_cell_add(ref_grid_pri(ref_grid),local,&cell),"add prism");
    }

  RSS( ref_node_initialize_n_global(ref_node,nnodesg), "glob" );

  global[0] = 1; global[1] = 0; global[2] = 3; global[3] = 4; global[4] = 20;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      RSS( ref_node_local(ref_node,global[3], &(local[3])),"loc");
      local[4]=global[4];
      RSS(ref_cell_add(ref_grid_qua(ref_grid),local,&cell),"add quad");
    }

  global[0] = 4; global[1] = 3; global[2] = 6; global[3] = 7; global[4] = 20;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      RSS( ref_node_local(ref_node,global[3], &(local[3])),"loc");
      local[4]=global[4];
      RSS(ref_cell_add(ref_grid_qua(ref_grid),local,&cell),"add quad");
    }

  global[0] = 7; global[1] = 6; global[2] = 9; global[3] =10; global[4] = 20;
  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) )
    {
      RSS( ref_node_local(ref_node,global[0], &(local[0])),"loc");
      RSS( ref_node_local(ref_node,global[1], &(local[1])),"loc");
      RSS( ref_node_local(ref_node,global[2], &(local[2])),"loc");
      RSS( ref_node_local(ref_node,global[3], &(local[3])),"loc");
      local[4]=global[4];
      RSS(ref_cell_add(ref_grid_qua(ref_grid),local,&cell),"add quad");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_hex_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[8] = {0,1,2,3,4,5,6,7};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,4,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,5,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,6,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,7,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS( ref_node_initialize_n_global( ref_node, 8 ), "init glob" );

  RSS(ref_cell_add(ref_grid_hex(ref_grid),nodes,&cell),"add prism");

  nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3; nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 4; nodes[1] = 7; nodes[2] = 6; nodes[3] = 5; nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 1; nodes[1] = 5; nodes[2] = 6; nodes[3] = 2; nodes[4] = 20;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 0; nodes[1] = 3; nodes[2] = 7; nodes[3] = 4; nodes[4] = 20;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 0; nodes[1] = 4; nodes[2] = 5; nodes[3] = 1; nodes[4] = 30;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 2; nodes[1] = 6; nodes[2] = 7; nodes[3] = 3; nodes[4] = 30;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_brick_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], cell;


  REF_INT l=5,m=3,n=7;
  REF_INT i, j, k;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 0.1;

  REF_DBL dx, dy;

  REF_DBL dz0 = 0.1;
  REF_DBL r = 1.2;

  dx = (x1-x0)/((REF_DBL)(l-1));
  dy = (y1-y0)/((REF_DBL)(m-1));

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

#define ijk2node(i,j,k) ((i) + (j)*l + (k)*l*m)

  for ( k = 0 ; k < n ; k++ )
    for ( j = 0 ; j < m ; j++ )
      for ( i = 0 ; i < l ; i++ )
	{
	  global = ijk2node(i,j,k);
	  RSS( ref_node_add( ref_node, global, &node ), "node");
	  ref_node_xyz(ref_node, 0, node ) = x0 + dx*(REF_DBL)i;
	  ref_node_xyz(ref_node, 1, node ) = y0 + dy*(REF_DBL)j;
	  ref_node_xyz(ref_node, 2, node ) = dz0*(1.0-pow(r,k))/(1.0-r);
	}

#define ijk2hex(i,j,k,hex)			  \
  (hex)[0] = ijk2node((i)-1,(j)-1,(k)-1);	  \
  (hex)[1] = ijk2node((i)  ,(j)-1,(k)-1);	  \
  (hex)[2] = ijk2node((i)  ,(j)  ,(k)-1);	  \
  (hex)[3] = ijk2node((i)-1,(j)  ,(k)-1);	  \
  (hex)[4] = ijk2node((i)-1,(j)-1,(k)  );	  \
  (hex)[5] = ijk2node((i)  ,(j)-1,(k)  );	  \
  (hex)[6] = ijk2node((i)  ,(j)  ,(k)  );	  \
  (hex)[7] = ijk2node((i)-1,(j)  ,(k)  );	  \
  

  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      for ( i = 1 ; i < l ; i++ )
	{
	  ijk2hex(i,j,k,hex);
	  RSS( ref_cell_add(ref_grid_hex(ref_grid),hex, &cell),"hex");
	}

  return REF_SUCCESS;
}
