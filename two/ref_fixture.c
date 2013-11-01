
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_fixture.h"
#include "ref_mpi.h"
#include "ref_part.h"
#include "ref_math.h"

#define add_that_node( node, x, y, z )					\
  RSS(ref_node_add(ref_node,global[(node)],&(local[(node)])),"add node"); \
  ref_node_xyz(ref_node,0,local[(node)]) = (x);				\
  ref_node_xyz(ref_node,1,local[(node)]) = (y);				\
  ref_node_xyz(ref_node,2,local[(node)]) = (z);				\
  ref_node_part(ref_node,local[(node)]) =				\
    ref_part_implicit( nnodesg, ref_mpi_n,				\
		       ref_node_global(ref_node,local[(node)]) );

REF_STATUS ref_fixture_tet_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global[REF_CELL_MAX_SIZE_PER];
  REF_INT local[REF_CELL_MAX_SIZE_PER];
  REF_INT nnodesg = 4;
  REF_INT cell;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);


  global[0]=0;global[1]=1;global[2]=2;global[3]=3;

  if ( ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[0] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[1] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[2] ) ||
       ref_mpi_id == ref_part_implicit( nnodesg, ref_mpi_n, global[3] ) )
    {
      add_that_node(0,0.0,0.0,0.0);
      add_that_node(1,1.0,0.0,0.0);
      add_that_node(2,0.0,1.0,0.0);
      add_that_node(3,0.0,0.0,1.0);

      RSS(ref_cell_add(ref_grid_tet(ref_grid),local,&cell),"add tet");
    }

  RSS( ref_node_initialize_n_global( ref_node, nnodesg ), "init glob" );

  global[0]=0;global[1]=1;global[2]=2;global[3]=10;
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

  ref_grid_twod(ref_grid) = REF_TRUE;
 
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
      add_that_node(1,0.0,0.0,1.0);
      add_that_node(2,1.0,0.0,0.0);
      add_that_node(3,0.0,1.0,0.0);
      add_that_node(4,0.0,1.0,1.0);
      add_that_node(5,1.0,1.0,0.0);

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
      add_that_node(3,0.3,0.3,2.0);

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

REF_STATUS ref_fixture_hex_brick_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], cell;
  REF_INT quad[5];


  REF_INT l=5,m=3,n=9;
  REF_INT i, j, k;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 0.1;

  REF_DBL dx, dy;

  REF_DBL dz0 = 0.0005;
  REF_DBL r = 1.5;

  dx = (x1-x0)/((REF_DBL)(l-1));
  dy = (y1-y0)/((REF_DBL)(m-1));

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

#define ijk2node(i,j,k,l,m,n) ((i) + (j)*(l) + (k)*(l)*(m))

  for ( k = 0 ; k < n ; k++ )
    for ( j = 0 ; j < m ; j++ )
      for ( i = 0 ; i < l ; i++ )
	{
	  global = ijk2node(i,j,k,l,m,n);
	  RSS( ref_node_add( ref_node, global, &node ), "node");
	  ref_node_xyz(ref_node, 0, node ) = x0 + dx*(REF_DBL)i;
	  ref_node_xyz(ref_node, 1, node ) = y0 + dy*(REF_DBL)j;
	  ref_node_xyz(ref_node, 2, node ) = dz0*(1.0-pow(r,k))/(1.0-r);
	}

#define ijk2hex(i,j,k,l,m,n,hex)			\
  (hex)[0] = ijk2node((i)-1,(j)-1,(k)-1,(l),(m),(n));	\
  (hex)[1] = ijk2node((i)  ,(j)-1,(k)-1,(l),(m),(n));	  \
  (hex)[2] = ijk2node((i)  ,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[3] = ijk2node((i)-1,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[4] = ijk2node((i)-1,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[5] = ijk2node((i)  ,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[6] = ijk2node((i)  ,(j)  ,(k)  ,(l),(m),(n));	  \
  (hex)[7] = ijk2node((i)-1,(j)  ,(k)  ,(l),(m),(n));	  \
  

  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      for ( i = 1 ; i < l ; i++ )
	{
	  ijk2hex(i,j,k,l,m,n,hex);
	  RSS( ref_cell_add(ref_grid_hex(ref_grid),hex, &cell),"hex");
	}

  quad[4]=1;
  i = 1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[3];
	quad[2]=hex[7];
	quad[3]=hex[4];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=2;
  i = l-1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[2];
	quad[1]=hex[1];
	quad[2]=hex[5];
	quad[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=3;
  j=1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[1];
	quad[1]=hex[0];
	quad[2]=hex[4];
	quad[3]=hex[5];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=4;
  j=m-1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[3];
	quad[1]=hex[2];
	quad[2]=hex[6];
	quad[3]=hex[7];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=5;
  k=1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[1];
	quad[2]=hex[2];
	quad[3]=hex[3];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=6;
  k=n-1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[5];
	quad[1]=hex[4];
	quad[2]=hex[7];
	quad[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pri_brick_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], pri[6], cell;
  REF_INT quad[5], tri[4];


  REF_INT l=21,m=3,n=31;
  REF_INT i, j, k;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 0.1;

  REF_DBL dx, dy;

  REF_DBL dz0 = 0.0001;
  REF_DBL r = 1.2;

  dx = (x1-x0)/((REF_DBL)(l-1));
  dy = (y1-y0)/((REF_DBL)(m-1));

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

#define ijk2node(i,j,k,l,m,n) ((i) + (j)*(l) + (k)*(l)*(m))

  for ( k = 0 ; k < n ; k++ )
    for ( j = 0 ; j < m ; j++ )
      for ( i = 0 ; i < l ; i++ )
	{
	  global = ijk2node(i,j,k,l,m,n);
	  RSS( ref_node_add( ref_node, global, &node ), "node");
	  ref_node_xyz(ref_node, 0, node ) = x0 + dx*(REF_DBL)i;
	  ref_node_xyz(ref_node, 1, node ) = y0 + dy*(REF_DBL)j;
	  ref_node_xyz(ref_node, 2, node ) = dz0*(1.0-pow(r,k))/(1.0-r);
	}

#define ijk2hex(i,j,k,l,m,n,hex)			\
  (hex)[0] = ijk2node((i)-1,(j)-1,(k)-1,(l),(m),(n));	\
  (hex)[1] = ijk2node((i)  ,(j)-1,(k)-1,(l),(m),(n));	  \
  (hex)[2] = ijk2node((i)  ,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[3] = ijk2node((i)-1,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[4] = ijk2node((i)-1,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[5] = ijk2node((i)  ,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[6] = ijk2node((i)  ,(j)  ,(k)  ,(l),(m),(n));	  \
  (hex)[7] = ijk2node((i)-1,(j)  ,(k)  ,(l),(m),(n));	  \
  

  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      for ( i = 1 ; i < l ; i++ )
	{
	  ijk2hex(i,j,k,l,m,n,hex);
	  pri[0] = hex[0];
	  pri[1] = hex[1];
	  pri[2] = hex[2];
	  pri[3] = hex[4];
	  pri[4] = hex[5];
	  pri[5] = hex[6];
	  RSS( ref_cell_add(ref_grid_pri(ref_grid),pri, &cell),"pri");
	  pri[0] = hex[0];
	  pri[1] = hex[2];
	  pri[2] = hex[3];
	  pri[3] = hex[4];
	  pri[4] = hex[6];
	  pri[5] = hex[7];
	  RSS( ref_cell_add(ref_grid_pri(ref_grid),pri, &cell),"pri");
	}

  quad[4]=1;
  i = 1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[3];
	quad[2]=hex[7];
	quad[3]=hex[4];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=2;
  i = l-1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[2];
	quad[1]=hex[1];
	quad[2]=hex[5];
	quad[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=3;
  j=1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[1];
	quad[1]=hex[0];
	quad[2]=hex[4];
	quad[3]=hex[5];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=4;
  j=m-1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[3];
	quad[1]=hex[2];
	quad[2]=hex[6];
	quad[3]=hex[7];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  tri[3]=5;
  k=1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	tri[0]=hex[0];
	tri[1]=hex[1];
	tri[2]=hex[2];
	RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");
	tri[0]=hex[0];
	tri[1]=hex[2];
	tri[2]=hex[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");
      }

  tri[3]=6;
  k=n-1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	tri[0]=hex[4];
	tri[1]=hex[6];
	tri[2]=hex[5];
	RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");
	tri[0]=hex[4];
	tri[1]=hex[7];
	tri[2]=hex[6];
	RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");
      }

  return REF_SUCCESS;
}

/*
  REF_DBL t0 = ref_math_pi;
  REF_DBL t1 = 0.0;
  REF_DBL r0 = 0.10;
  REF_DBL r1 = 0.15;
  REF_DBL y0 = 0.0;
  REF_DBL y1 = 0.5;

  dt = (t1-t0)/((REF_DBL)(l-1));
  dr = (r1-r0)/((REF_DBL)(m-1));
  dy = (y1-y0)/((REF_DBL)(n-1));

	  t = t0 + dt*(REF_DBL)i;
	  y = y0 + dy*(REF_DBL)j;
	  r = r0 + dr*(REF_DBL)k;
	  ref_node_xyz(ref_node, 0, node ) = r*cos(t);
	  ref_node_xyz(ref_node, 1, node ) = y;
	  ref_node_xyz(ref_node, 2, node ) = r*sin(t);
 */

REF_STATUS ref_fixture_tet_brick_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], tet[4], cell;
  REF_INT quad[5], tri[4];


  REF_INT l=4,m=4,n=4;
  REF_INT i, j, k;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 1.0;

  REF_DBL z0 = 0.0;
  REF_DBL z1 = 1.0;

  REF_DBL dx, dy, dz;

  dx = (x1-x0)/((REF_DBL)(l-1));
  dy = (y1-y0)/((REF_DBL)(m-1));
  dz = (z1-z0)/((REF_DBL)(n-1));

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

#define ijk2node(i,j,k,l,m,n) ((i) + (j)*(l) + (k)*(l)*(m))

  for ( k = 0 ; k < n ; k++ )
    for ( j = 0 ; j < m ; j++ )
      for ( i = 0 ; i < l ; i++ )
	{
	  global = ijk2node(i,j,k,l,m,n);
	  RSS( ref_node_add( ref_node, global, &node ), "node");
	  ref_node_xyz(ref_node, 0, node ) = x0 + dx*(REF_DBL)i;
	  ref_node_xyz(ref_node, 1, node ) = y0 + dy*(REF_DBL)j;
	  ref_node_xyz(ref_node, 2, node ) = z0 + dz*(REF_DBL)k;
	}

#define ijk2hex(i,j,k,l,m,n,hex)			\
  (hex)[0] = ijk2node((i)-1,(j)-1,(k)-1,(l),(m),(n));	\
  (hex)[1] = ijk2node((i)  ,(j)-1,(k)-1,(l),(m),(n));	  \
  (hex)[2] = ijk2node((i)  ,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[3] = ijk2node((i)-1,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[4] = ijk2node((i)-1,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[5] = ijk2node((i)  ,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[6] = ijk2node((i)  ,(j)  ,(k)  ,(l),(m),(n));	  \
  (hex)[7] = ijk2node((i)-1,(j)  ,(k)  ,(l),(m),(n));	  \
  
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      for ( i = 1 ; i < l ; i++ )
	{
	  ijk2hex(i,j,k,l,m,n,hex);
	  tet[0] = hex[0];
	  tet[1] = hex[1];
	  tet[2] = hex[3];
	  tet[3] = hex[4];
	  RSS( ref_cell_add(ref_grid_tet(ref_grid),tet, &cell),"tet");
	  tet[0] = hex[1];
	  tet[1] = hex[3];
	  tet[2] = hex[4];
	  tet[3] = hex[5];
	  RSS( ref_cell_add(ref_grid_tet(ref_grid),tet, &cell),"tet");
	  tet[0] = hex[3];
	  tet[1] = hex[4];
	  tet[2] = hex[5];
	  tet[3] = hex[7];
	  RSS( ref_cell_add(ref_grid_tet(ref_grid),tet, &cell),"tet");
	  tet[0] = hex[1];
	  tet[1] = hex[2];
	  tet[2] = hex[3];
	  tet[3] = hex[5];
	  RSS( ref_cell_add(ref_grid_tet(ref_grid),tet, &cell),"tet");
	  tet[0] = hex[2];
	  tet[1] = hex[3];
	  tet[2] = hex[5];
	  tet[3] = hex[7];
	  RSS( ref_cell_add(ref_grid_tet(ref_grid),tet, &cell),"tet");
	  tet[0] = hex[2];
	  tet[1] = hex[7];
	  tet[2] = hex[5];
	  tet[3] = hex[6];
	  RSS( ref_cell_add(ref_grid_tet(ref_grid),tet, &cell),"tet");
	}

  quad[4]=1;tri[3]=quad[4];
  i = 1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[3];
	quad[2]=hex[7];
	quad[3]=hex[4];
	tri[0] = quad[0];
	tri[1] = quad[1];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
	tri[0] = quad[1];
	tri[1] = quad[2];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
      }

  quad[4]=2;tri[3]=quad[4];
  i = l-1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[2];
	quad[1]=hex[1];
	quad[2]=hex[5];
	quad[3]=hex[6];
	tri[0] = quad[0];
	tri[1] = quad[1];
	tri[2] = quad[2];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
	tri[0] = quad[0];
	tri[1] = quad[2];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
       }

  quad[4]=3;tri[3]=quad[4];
  j=1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[1];
	quad[1]=hex[0];
	quad[2]=hex[4];
	quad[3]=hex[5];
	tri[0] = quad[0];
	tri[1] = quad[1];
	tri[2] = quad[2];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
	tri[0] = quad[0];
	tri[1] = quad[2];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
      }

  quad[4]=4;tri[3]=quad[4];
  j=m-1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[3];
	quad[1]=hex[2];
	quad[2]=hex[6];
	quad[3]=hex[7];
	tri[0] = quad[0];
	tri[1] = quad[1];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
	tri[0] = quad[1];
	tri[1] = quad[2];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
      }

  quad[4]=5;tri[3]=quad[4];
  k=1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[1];
	quad[2]=hex[2];
	quad[3]=hex[3];
	tri[0] = quad[0];
	tri[1] = quad[1];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
	tri[0] = quad[1];
	tri[1] = quad[2];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
      }

  quad[4]=6;tri[3]=quad[4];
  k=n-1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[5];
	quad[1]=hex[4];
	quad[2]=hex[7];
	quad[3]=hex[6];
	tri[0] = quad[0];
	tri[1] = quad[1];
	tri[2] = quad[2];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
	tri[0] = quad[0];
	tri[1] = quad[2];
	tri[2] = quad[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
      }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_twod_brick_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], pri[6], cell;
  REF_INT qua[5], tri[4];


  REF_INT l=4,m=2,n=4;
  REF_INT i, j, k;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 1.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 1.0;

  REF_DBL z0 = 0.0;
  REF_DBL z1 = 1.0;

  REF_DBL dx, dy, dz;

  dx = (x1-x0)/((REF_DBL)(l-1));
  dy = (y1-y0)/((REF_DBL)(m-1));
  dz = (z1-z0)/((REF_DBL)(n-1));

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  ref_grid_twod(ref_grid) = REF_TRUE;

  /*
                               inode7-----11-----inode6
                                 /.                /|
                                / .               / |
                               /  .              /  |
                              /   .             /   |
                             9    .        /  10    6
                            /     7           /     |
                           /      .          /      |
                          /       .         /       |
                       inode4-8----------inode5     |
                         |      inode3.....|...5..inode2
                         |       .         |       /
                         |      .          |      /
                         |     .           |     /
                         2    1    /       4    3
 z                       |   .             |   /
 ^   y                   |  .              |  /
 |  /                    | .               | /
 | /                     |.                |/
 |/                    inode0------0-----inode1
 +----> x
*/

#define ijk2node(i,j,k,l,m,n) ((i) + (j)*(l) + (k)*(l)*(m))

  for ( k = 0 ; k < n ; k++ )
    for ( j = 0 ; j < m ; j++ )
      for ( i = 0 ; i < l ; i++ )
	{
	  global = ijk2node(i,j,k,l,m,n);
	  RSS( ref_node_add( ref_node, global, &node ), "node");
	  ref_node_xyz(ref_node, 0, node ) = x0 + dx*(REF_DBL)i;
	  ref_node_xyz(ref_node, 1, node ) = y0 + dy*(REF_DBL)j;
	  ref_node_xyz(ref_node, 2, node ) = z0 + dz*(REF_DBL)k;
	}

#define ijk2hex(i,j,k,l,m,n,hex)			\
  (hex)[0] = ijk2node((i)-1,(j)-1,(k)-1,(l),(m),(n));	\
  (hex)[1] = ijk2node((i)  ,(j)-1,(k)-1,(l),(m),(n));	  \
  (hex)[2] = ijk2node((i)  ,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[3] = ijk2node((i)-1,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[4] = ijk2node((i)-1,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[5] = ijk2node((i)  ,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[6] = ijk2node((i)  ,(j)  ,(k)  ,(l),(m),(n));	  \
  (hex)[7] = ijk2node((i)-1,(j)  ,(k)  ,(l),(m),(n));	  \
  
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      for ( i = 1 ; i < l ; i++ )
	{
	  ijk2hex(i,j,k,l,m,n,hex);
	  pri[0] = hex[0];
	  pri[1] = hex[4];
	  pri[2] = hex[5];
	  pri[3] = hex[3];
	  pri[4] = hex[7];
	  pri[5] = hex[6];
	  RSS( ref_cell_add(ref_grid_pri(ref_grid),pri, &cell),"pri");
	  pri[0] = hex[0];
	  pri[1] = hex[5];
	  pri[2] = hex[1];
	  pri[3] = hex[3];
	  pri[4] = hex[6];
	  pri[5] = hex[2];
	  RSS( ref_cell_add(ref_grid_pri(ref_grid),pri, &cell),"pri");
	}

  qua[4]=3;
  i = 1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	qua[0]=hex[0];
	qua[1]=hex[3];
	qua[2]=hex[7];
	qua[3]=hex[4];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),qua, &cell),"qua");
      }

  qua[4]=4;
  i = l-1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	qua[0]=hex[2];
	qua[1]=hex[1];
	qua[2]=hex[5];
	qua[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),qua, &cell),"qua");
       }

  qua[4]=1;tri[3]=qua[4];
  j=1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	qua[0]=hex[1];
	qua[1]=hex[0];
	qua[2]=hex[4];
	qua[3]=hex[5];
	tri[0] = qua[0];
	tri[1] = qua[1];
	tri[2] = qua[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
	tri[0] = qua[1];
	tri[1] = qua[2];
	tri[2] = qua[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
      }

  qua[4]=2;tri[3]=qua[4];
  j=m-1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	qua[0]=hex[3];
	qua[1]=hex[2];
	qua[2]=hex[6];
	qua[3]=hex[7];
	tri[0] = qua[0];
	tri[1] = qua[2];
	tri[2] = qua[3];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
	tri[0] = qua[0];
	tri[1] = qua[1];
	tri[2] = qua[2];
	RSS( ref_cell_add(ref_grid_tri(ref_grid),tri, &cell),"qua");
      }

  qua[4]=5;
  k=1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	qua[0]=hex[0];
	qua[1]=hex[1];
	qua[2]=hex[2];
	qua[3]=hex[3];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),qua, &cell),"qua");
      }

  qua[4]=6;
  k=n-1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	qua[0]=hex[5];
	qua[1]=hex[4];
	qua[2]=hex[7];
	qua[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),qua, &cell),"qua");
      }

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_clock( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT node, nodes[8], cell;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  /*
#!/usr/bin/env ruby

$:.push('../util')

require 'pslg'

nodes = Array.new

r = 1.0
(0..11).each do |t|
  th=2.0*Math::PI*(1.0/12.0)*t
  x = -r*Math.sin(th)
  y =  r*Math.cos(th)
  nodes << Node.new(x,y)
end

msh=Msh.new
11.times do |s|
  msh << Segment.new(3,nodes[s],nodes[s+1])
end
msh << Segment.new(3,nodes[11],nodes[0])

msh.export('clock_g.msh')
   */

  RSS(ref_node_add(ref_node,0,&node),"node");
  ref_node_xyz(ref_node,0,node) = -0.000000000000000e+00;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 1.000000000000000e+00;
  RSS(ref_node_add(ref_node,1,&node),"node");
  ref_node_xyz(ref_node,0,node) = -5.000000000000000e-01;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 8.660254037844390e-01;
  RSS(ref_node_add(ref_node,2,&node),"node");
  ref_node_xyz(ref_node,0,node) = -8.660254037844390e-01;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 5.000000000000000e-01;
  RSS(ref_node_add(ref_node,3,&node),"node");
  ref_node_xyz(ref_node,0,node) = -1.000000000000000e+00;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 6.123233995736770e-17;
  RSS(ref_node_add(ref_node,4,&node),"node");
  ref_node_xyz(ref_node,0,node) = -8.660254037844390e-01;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -5.000000000000000e-01;
  RSS(ref_node_add(ref_node,5,&node),"node");
  ref_node_xyz(ref_node,0,node) = -5.000000000000000e-01;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -8.660254037844380e-01;
  RSS(ref_node_add(ref_node,6,&node),"node");
  ref_node_xyz(ref_node,0,node) = -1.224646799147350e-16;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -1.000000000000000e+00;
  RSS(ref_node_add(ref_node,7,&node),"node");
  ref_node_xyz(ref_node,0,node) = 5.000000000000000e-01;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -8.660254037844390e-01;
  RSS(ref_node_add(ref_node,8,&node),"node");
  ref_node_xyz(ref_node,0,node) = 8.660254037844380e-01;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -5.000000000000000e-01;
  RSS(ref_node_add(ref_node,9,&node),"node");
  ref_node_xyz(ref_node,0,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -1.836970198721030e-16;
  RSS(ref_node_add(ref_node,10,&node),"node");
  ref_node_xyz(ref_node,0,node) = 8.660254037844390e-01;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 4.999999999999990e-01;
  RSS(ref_node_add(ref_node,11,&node),"node");
  ref_node_xyz(ref_node,0,node) = 5.000000000000000e-01;
  ref_node_xyz(ref_node,1,node) = 0.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 8.660254037844380e-01;
  RSS(ref_node_add(ref_node,12,&node),"node");
  ref_node_xyz(ref_node,0,node) = -0.000000000000000e+00;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 1.000000000000000e+00;
  RSS(ref_node_add(ref_node,13,&node),"node");
  ref_node_xyz(ref_node,0,node) = -5.000000000000000e-01;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 8.660254037844390e-01;
  RSS(ref_node_add(ref_node,14,&node),"node");
  ref_node_xyz(ref_node,0,node) = -8.660254037844390e-01;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 5.000000000000000e-01;
  RSS(ref_node_add(ref_node,15,&node),"node");
  ref_node_xyz(ref_node,0,node) = -1.000000000000000e+00;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 6.123233995736770e-17;
  RSS(ref_node_add(ref_node,16,&node),"node");
  ref_node_xyz(ref_node,0,node) = -8.660254037844390e-01;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -5.000000000000000e-01;
  RSS(ref_node_add(ref_node,17,&node),"node");
  ref_node_xyz(ref_node,0,node) = -5.000000000000000e-01;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -8.660254037844380e-01;
  RSS(ref_node_add(ref_node,18,&node),"node");
  ref_node_xyz(ref_node,0,node) = -1.224646799147350e-16;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -1.000000000000000e+00;
  RSS(ref_node_add(ref_node,19,&node),"node");
  ref_node_xyz(ref_node,0,node) = 5.000000000000000e-01;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -8.660254037844390e-01;
  RSS(ref_node_add(ref_node,20,&node),"node");
  ref_node_xyz(ref_node,0,node) = 8.660254037844380e-01;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -5.000000000000000e-01;
  RSS(ref_node_add(ref_node,21,&node),"node");
  ref_node_xyz(ref_node,0,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = -1.836970198721030e-16;
  RSS(ref_node_add(ref_node,22,&node),"node");
  ref_node_xyz(ref_node,0,node) = 8.660254037844390e-01;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 4.999999999999990e-01;
  RSS(ref_node_add(ref_node,23,&node),"node");
  ref_node_xyz(ref_node,0,node) = 5.000000000000000e-01;
  ref_node_xyz(ref_node,1,node) = 1.000000000000000e+00;
  ref_node_xyz(ref_node,2,node) = 8.660254037844380e-01;
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 13;
  nodes[3] = 12;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 1;
  nodes[1] = 2;
  nodes[2] = 14;
  nodes[3] = 13;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 2;
  nodes[1] = 3;
  nodes[2] = 15;
  nodes[3] = 14;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 3;
  nodes[1] = 4;
  nodes[2] = 16;
  nodes[3] = 15;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 4;
  nodes[1] = 5;
  nodes[2] = 17;
  nodes[3] = 16;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 5;
  nodes[1] = 6;
  nodes[2] = 18;
  nodes[3] = 17;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 6;
  nodes[1] = 7;
  nodes[2] = 19;
  nodes[3] = 18;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 7;
  nodes[1] = 8;
  nodes[2] = 20;
  nodes[3] = 19;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 8;
  nodes[1] = 9;
  nodes[2] = 21;
  nodes[3] = 20;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 9;
  nodes[1] = 10;
  nodes[2] = 22;
  nodes[3] = 21;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 10;
  nodes[1] = 11;
  nodes[2] = 23;
  nodes[3] = 22;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");
  nodes[0] = 11;
  nodes[1] = 0;
  nodes[2] = 12;
  nodes[3] = 23;
  nodes[4] = 3;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"qua");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_boom2d_grid( REF_GRID *ref_grid_ptr, 
				    REF_DBL theta_deg,
				    REF_DBL beta_deg,
				    REF_INT nx, REF_INT nz )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], cell;
  REF_INT quad[5];

  REF_INT l=8*nx+1,m=2,n=10*nz+1;
  REF_INT i, j, k;

  REF_DBL x0 = -4.0;
  REF_DBL x1 =  4.0;

  REF_DBL y0 = 0.0;
  REF_DBL y1 = 1.0;

  REF_DBL z0 = 0.0;
  REF_DBL z1 =20.0;

  REF_DBL dx, dy, dz;

  REF_DBL tan_beta;

  tan_beta = tan(ref_math_in_radians(theta_deg));
  tan_beta = tan(ref_math_in_radians(beta_deg));
  
  dx = (x1-x0)/((REF_DBL)(l-1));
  dy = (y1-y0)/((REF_DBL)(m-1));
  dz = (z1-z0)/((REF_DBL)(n-1));

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

#define ijk2node(i,j,k,l,m,n) ((i) + (j)*(l) + (k)*(l)*(m))

  for ( k = 0 ; k < n ; k++ )
    for ( j = 0 ; j < m ; j++ )
      for ( i = 0 ; i < l ; i++ )
	{
	  global = ijk2node(i,j,k,l,m,n);
	  RSS( ref_node_add( ref_node, global, &node ), "node");
	  ref_node_xyz(ref_node, 0, node ) = x0 + dx*(REF_DBL)i;
	  ref_node_xyz(ref_node, 1, node ) = y0 + dy*(REF_DBL)j;
	  ref_node_xyz(ref_node, 2, node ) = z0 + dz*(REF_DBL)k;
	  /* ramp
	  if ( ref_node_xyz(ref_node, 0, node ) > 0.0 )
	    ref_node_xyz(ref_node, 2, node )
	      += tan_theta * ref_node_xyz(ref_node, 0, node );
 */
	  /* cosine */
	  if ( ref_node_xyz(ref_node, 0, node ) > -2.0 &&
	       ref_node_xyz(ref_node, 0, node ) < 2.0 )
	    {
	      ref_node_xyz(ref_node, 2, node )
		-= 0.005 * (1+cos(ref_math_pi*ref_node_xyz(ref_node, 0, node )/2.0));
	    }
	  

	  /* shear */
	  if ( ABS(beta_deg) > 1.0e-8 )
	    ref_node_xyz(ref_node, 0, node ) 
	      += ref_node_xyz(ref_node, 2, node ) / tan_beta;
	}

#define ijk2hex(i,j,k,l,m,n,hex)			\
  (hex)[0] = ijk2node((i)-1,(j)-1,(k)-1,(l),(m),(n));	\
  (hex)[1] = ijk2node((i)  ,(j)-1,(k)-1,(l),(m),(n));	  \
  (hex)[2] = ijk2node((i)  ,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[3] = ijk2node((i)-1,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[4] = ijk2node((i)-1,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[5] = ijk2node((i)  ,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[6] = ijk2node((i)  ,(j)  ,(k)  ,(l),(m),(n));	  \
  (hex)[7] = ijk2node((i)-1,(j)  ,(k)  ,(l),(m),(n));	  \
  

  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      for ( i = 1 ; i < l ; i++ )
	{
	  ijk2hex(i,j,k,l,m,n,hex);
	  RSS( ref_cell_add(ref_grid_hex(ref_grid),hex, &cell),"hex");
	}

  quad[4]=1;
  i = 1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[3];
	quad[2]=hex[7];
	quad[3]=hex[4];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=2;
  i = l-1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[2];
	quad[1]=hex[1];
	quad[2]=hex[5];
	quad[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=3;
  j=1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[1];
	quad[1]=hex[0];
	quad[2]=hex[4];
	quad[3]=hex[5];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=4;
  j=m-1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[3];
	quad[1]=hex[2];
	quad[2]=hex[6];
	quad[3]=hex[7];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=5;
  k=1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l/2 ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[1];
	quad[2]=hex[2];
	quad[3]=hex[3];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=6;
  k=1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = l/2 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[1];
	quad[2]=hex[2];
	quad[3]=hex[3];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=7;
  k=n-1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[5];
	quad[1]=hex[4];
	quad[2]=hex[7];
	quad[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  return REF_SUCCESS;
}


REF_STATUS ref_fixture_boom3d_grid( REF_GRID *ref_grid_ptr,
				    REF_INT nx,
				    REF_INT nt,
				    REF_INT nr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT global, node, hex[8], cell;
  REF_INT quad[5];

  REF_INT l=nx,m=nt,n=nr;
  REF_INT i, j, k;

  REF_DBL x0 = 0.0;
  REF_DBL x1 = 4.0;
  REF_DBL t0 = 0.0;
  REF_DBL t1 = -0.5*ref_math_pi;
  REF_DBL r0 = 1.0;
  REF_DBL r1 = 4.0;

  REF_DBL mach, mu, tan_mu;

  REF_DBL x, t, r;
  REF_DBL dx, dt, dr;

  mach = 1.6;
  mu = asin(1.0/mach);
  tan_mu = tan(mu);
  printf("mach %f mu %f tan %f\n",mach,mu,tan_mu);

  dx = (x1-x0)/((REF_DBL)(l-1));
  dt = (t1-t0)/((REF_DBL)(m-1));
  dr = (r1-r0)/((REF_DBL)(n-1));

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

#define ijk2node(i,j,k,l,m,n) ((i) + (j)*(l) + (k)*(l)*(m))

  for ( k = 0 ; k < n ; k++ )
    for ( j = 0 ; j < m ; j++ )
      for ( i = 0 ; i < l ; i++ )
	{
	  global = ijk2node(i,j,k,l,m,n);
	  RSS( ref_node_add( ref_node, global, &node ), "node");
	  x = x0 + dx*(REF_DBL)i;
	  t = t0 + dt*(REF_DBL)j;
	  r = r0 + dr*(REF_DBL)k;

	  /* cosine */
	  if ( x > 1.0 &&
	       x < 3.0 &&
	       k == 0 )
	    {
	      r -= 0.005 * (1+cos(2.0*ref_math_pi*((x-2.0)/2.0)));
	    }

	  ref_node_xyz(ref_node, 0, node ) = x;
	  ref_node_xyz(ref_node, 1, node ) = r*cos(t);
	  ref_node_xyz(ref_node, 2, node ) = r*sin(t);

	  /* shear */
	  ref_node_xyz(ref_node, 0, node ) += r / tan_mu;
	}

#define ijk2hex(i,j,k,l,m,n,hex)			\
  (hex)[0] = ijk2node((i)-1,(j)-1,(k)-1,(l),(m),(n));	\
  (hex)[1] = ijk2node((i)  ,(j)-1,(k)-1,(l),(m),(n));	  \
  (hex)[2] = ijk2node((i)  ,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[3] = ijk2node((i)-1,(j)  ,(k)-1,(l),(m),(n));	  \
  (hex)[4] = ijk2node((i)-1,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[5] = ijk2node((i)  ,(j)-1,(k)  ,(l),(m),(n));	  \
  (hex)[6] = ijk2node((i)  ,(j)  ,(k)  ,(l),(m),(n));	  \
  (hex)[7] = ijk2node((i)-1,(j)  ,(k)  ,(l),(m),(n));	  \
  

  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      for ( i = 1 ; i < l ; i++ )
	{
	  ijk2hex(i,j,k,l,m,n,hex);
	  RSS( ref_cell_add(ref_grid_hex(ref_grid),hex, &cell),"hex");
	}

  quad[4]=1;
  i = 1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[3];
	quad[2]=hex[7];
	quad[3]=hex[4];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=2;
  i = l-1;
  for ( k = 1 ; k < n ; k++ )
    for ( j = 1 ; j < m ; j++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[2];
	quad[1]=hex[1];
	quad[2]=hex[5];
	quad[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=3;
  j=1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[1];
	quad[1]=hex[0];
	quad[2]=hex[4];
	quad[3]=hex[5];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=4;
  j=m-1;
  for ( k = 1 ; k < n ; k++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[3];
	quad[1]=hex[2];
	quad[2]=hex[6];
	quad[3]=hex[7];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=5;
  k=1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[0];
	quad[1]=hex[1];
	quad[2]=hex[2];
	quad[3]=hex[3];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  quad[4]=6;
  k=n-1;
  for ( j = 1 ; j < m ; j++ )
    for ( i = 1 ; i < l ; i++ )
      {
	ijk2hex(i,j,k,l,m,n,hex);
	quad[0]=hex[5];
	quad[1]=hex[4];
	quad[2]=hex[7];
	quad[3]=hex[6];
	RSS( ref_cell_add(ref_grid_qua(ref_grid),quad, &cell),"qua");
      }

  return REF_SUCCESS;
}
