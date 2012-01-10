
#include <stdlib.h>
#include <stdio.h>

#include "ref_swap.h"

REF_STATUS ref_swap_remove_two_face_cell( REF_GRID ref_grid, REF_INT cell )
{
  REF_CELL ref_cell;
  REF_INT cell_face;
  REF_INT node;
  REF_INT cell_nodes[4];
  REF_INT face_nodes[4];
  REF_INT found;
  REF_INT face0, face1;
  REF_INT cell_face0, cell_face1;
  REF_INT faceid0, faceid1;
  REF_INT temp, face;

  ref_cell = ref_grid_tet(ref_grid);
 
  face0 = REF_EMPTY;
  face1 = REF_EMPTY;
  faceid0 = REF_EMPTY;
  faceid1 = REF_EMPTY;
  cell_face0 = REF_EMPTY;
  cell_face1 = REF_EMPTY;
  for (cell_face=0;cell_face<4;cell_face++)
    {
      for(node=0;node<4;node++)
	cell_nodes[node]=ref_cell_f2n(ref_cell,node,cell,cell_face);
      if ( REF_SUCCESS == ref_cell_with( ref_grid_tri( ref_grid ), 
					 cell_nodes, &found ) )
	{
	  RSS( ref_cell_nodes( ref_grid_tri( ref_grid ), 
			       found, face_nodes), "tri");
	  if ( REF_EMPTY == face0 )
	    {
	      face0 = found;
	      faceid0 = face_nodes[3];
	      cell_face0 = cell_face;
	    }
	  else
	    {
	      if ( REF_EMPTY != face1 )
		{
		  RSS(REF_INVALID,"three or more faces detected");
		}
	      face1 = found;
	      faceid1 = face_nodes[3];
	      cell_face1 = cell_face;
	    }
	}
    }
  
  if ( REF_EMPTY == face0 ) return REF_INVALID;
  if ( REF_EMPTY == face1 ) return REF_INVALID;
  if ( faceid0 != faceid1 ) return REF_INVALID;

  RSS( ref_cell_remove( ref_grid_tri( ref_grid ), face0 ), "remove tri0" );
  RSS( ref_cell_remove( ref_grid_tri( ref_grid ), face1 ), "remove tri1" );

  for (cell_face=0;cell_face<4;cell_face++)
    {
      for(node=0;node<4;node++)
	cell_nodes[node]=ref_cell_f2n(ref_cell,node,cell,cell_face);
      if ( cell_face != cell_face0 && cell_face != cell_face1 )
	{
	  cell_nodes[3] = faceid0;
	  temp = cell_nodes[0]; 
	  cell_nodes[0] = cell_nodes[1];
	  cell_nodes[1] = temp;
	  RSS( ref_cell_add( ref_grid_tri( ref_grid ), cell_nodes, &face ), 
	       "add tri" );
	}
    }

  RSS( ref_cell_remove( ref_cell, cell ), "remove tet" );

  return REF_SUCCESS;
}
