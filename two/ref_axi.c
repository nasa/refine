
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_axi.h"

#include "ref_math.h"

#include "ref_node.h"
#include "ref_cell.h"

#include "ref_sort.h"

REF_STATUS ref_axi_wedge( REF_GRID ref_grid )
{
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_DBL pole_tol;
  REF_INT *o2n, node;
  REF_INT nhalf;

  REF_INT cell, nodes[8], new_nodes[8], pyr_nodes[8];
  REF_INT nunique, unique[8];
  REF_INT new_cell;
  
  REF_DBL radius, wedge_angle;

  ref_node = ref_grid_node(ref_grid);

  o2n = (REF_INT *)malloc( ref_node_n(ref_node) * sizeof(REF_INT));
  RNS(o2n,"malloc o2n NULL");

  pole_tol = 1.0e-6;
  wedge_angle = ref_math_in_radians(5.0);

  nhalf = ref_node_n(ref_node)/2;
  each_ref_node_valid_node( ref_node, node )
    {
      if ( ABS(ref_node_xyz(ref_node,2,node)) < pole_tol &&
	   ABS(ref_node_xyz(ref_node,1,node)) > 0.5 )
	{
	  o2n[node] = node-nhalf;
	  RSS( ref_node_remove( ref_node, node ), "remove" );
	}
      else
	{
	  o2n[node] = node;
	}
      if (  ABS(ref_node_xyz(ref_node,1,node)) > 0.5 )
	{
	  radius = ref_node_xyz(ref_node,2,node);
	  ref_node_xyz(ref_node,1,node) = radius * sin( wedge_angle );
	  ref_node_xyz(ref_node,2,node) = radius * cos( wedge_angle );
	}
    }

  ref_cell = ref_grid_tri(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell,cell,nodes)
    {
      for (node=0;node<3;node++)
	new_nodes[node] = o2n[nodes[node]];
      new_nodes[3] = nodes[3];
      RSS( ref_cell_replace( ref_cell, cell, new_nodes ), "renum" );
    }

  ref_cell = ref_grid_qua(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell,cell,nodes)
    {
      for (node=0;node<4;node++)
	new_nodes[node] = o2n[nodes[node]];
      new_nodes[4] = nodes[4];
      RSS( ref_sort_unique_int( 4, new_nodes, 
				&nunique, unique), "uniq" );
      if ( 4 > nunique ) RSS( ref_cell_remove( ref_cell, cell ), "rm qua" );
      if ( 3 == nunique )
	{
	  if ( new_nodes[2] == new_nodes[3] )
	    {
	      new_nodes[3]=new_nodes[4];
	    }
	  if ( new_nodes[1] == new_nodes[2] )
	    {
	      new_nodes[2]=new_nodes[3];
	      new_nodes[3]=new_nodes[4];
	    }
	  if ( new_nodes[0] == new_nodes[1] )
	    {
	      new_nodes[1]=new_nodes[2];
	      new_nodes[2]=new_nodes[3];
	      new_nodes[3]=new_nodes[4];
	    }
	  if ( new_nodes[3] == new_nodes[0] )
	    {
	      new_nodes[0]=new_nodes[1];
	      new_nodes[1]=new_nodes[2];
	      new_nodes[2]=new_nodes[3];
	      new_nodes[3]=new_nodes[4];
	    }
	  RSS( ref_cell_add( ref_grid_tri(ref_grid), 
			     new_nodes, &new_cell ), "new cell" );
	}
    }

  ref_cell = ref_grid_pri(ref_grid);

  each_ref_cell_valid_cell_with_nodes(ref_cell,cell,nodes)
    {
      for (node=0;node<6;node++)
	new_nodes[node] = o2n[nodes[node]];
      RSS( ref_sort_unique_int( 6, new_nodes, 
				&nunique, unique), "uniq" );
      if ( 6 > nunique ) RSS( ref_cell_remove( ref_cell, cell ), "rm qua" );
      if ( 5 == nunique )
	{
	  if ( new_nodes[0] == new_nodes[3] )
	    {
	      pyr_nodes[0] = new_nodes[1];
	      pyr_nodes[1] = new_nodes[2];
	      pyr_nodes[2] = new_nodes[0];
	      pyr_nodes[3] = new_nodes[4];
	      pyr_nodes[4] = new_nodes[5];
	    }
	  if ( new_nodes[1] == new_nodes[4] )
	    {
	      pyr_nodes[0] = new_nodes[2];
	      pyr_nodes[1] = new_nodes[0];
	      pyr_nodes[2] = new_nodes[1];
	      pyr_nodes[3] = new_nodes[5];
	      pyr_nodes[4] = new_nodes[3];
	    }
	  if ( new_nodes[2] == new_nodes[5] )
	    {
	      pyr_nodes[0] = new_nodes[0];
	      pyr_nodes[1] = new_nodes[1];
	      pyr_nodes[2] = new_nodes[2];
	      pyr_nodes[3] = new_nodes[3];
	      pyr_nodes[4] = new_nodes[4];
	    }
	  RSS( ref_cell_add( ref_grid_pyr(ref_grid), 
			     pyr_nodes, &new_cell ), "new cell" );
	}
      if ( 4 == nunique )
	{
	  if ( new_nodes[1] != new_nodes[4] ) new_nodes[3] = new_nodes[4];
	  if ( new_nodes[2] != new_nodes[5] ) new_nodes[3] = new_nodes[5];
	  RSS( ref_cell_add( ref_grid_tet(ref_grid), 
			     new_nodes, &new_cell ), "new cell" );
	}
    }

  free(o2n);

  return REF_SUCCESS;
}

