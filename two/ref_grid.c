
#include <stdlib.h>
#include <stdio.h>

#include "ref_grid.h"
#include "ref_adj.h"

#include "ref_malloc.h"
#include "ref_matrix.h"

REF_STATUS ref_grid_create( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;

  ref_malloc( *ref_grid_ptr, 1, REF_GRID_STRUCT );

  ref_grid = *ref_grid_ptr;

  RSS( ref_node_create( &ref_grid_node(ref_grid) ), "node create" );

  RSS( ref_cell_create( &ref_grid_tet(ref_grid), 4, REF_FALSE ), "tet create" );
  RSS( ref_cell_create( &ref_grid_pyr(ref_grid), 5, REF_FALSE ), "pyr create" );
  RSS( ref_cell_create( &ref_grid_pri(ref_grid), 6, REF_FALSE ), "pri create" );
  RSS( ref_cell_create( &ref_grid_hex(ref_grid), 8, REF_FALSE ), "hex create" );

  ref_grid_cell(ref_grid,4) = NULL;

  RSS( ref_cell_create( &ref_grid_edg(ref_grid), 2, REF_TRUE ), "edg create" );
  RSS( ref_cell_create( &ref_grid_tri(ref_grid), 3, REF_TRUE ), "tri create" );
  RSS( ref_cell_create( &ref_grid_qua(ref_grid), 4, REF_TRUE ), "qua create" );

  RSS( ref_geom_create( &ref_grid_geom(ref_grid) ), "geom create" );

  ref_grid_twod(ref_grid) = REF_FALSE;

  return REF_SUCCESS;
}

REF_STATUS ref_grid_deep_copy( REF_GRID *ref_grid_ptr, REF_GRID original )
{
  REF_GRID ref_grid;

  ref_malloc( *ref_grid_ptr, 1, REF_GRID_STRUCT );

  ref_grid = *ref_grid_ptr;

  RSS( ref_node_deep_copy( &ref_grid_node(ref_grid),
			   ref_grid_node(original) ), "node deep copy" );

  RSS( ref_cell_deep_copy( &ref_grid_tet(ref_grid),
			   ref_grid_tet(original) ), "tet deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_pyr(ref_grid),
			   ref_grid_pyr(original) ), "pyr deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_pri(ref_grid),
			   ref_grid_pri(original) ), "pri deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_hex(ref_grid),
			   ref_grid_hex(original) ), "hex deep copy" );

  ref_grid_cell(ref_grid,4) = NULL;

  RSS( ref_cell_deep_copy( &ref_grid_edg(ref_grid),
			   ref_grid_edg(original) ), "edg deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_tri(ref_grid),
			   ref_grid_tri(original) ), "tri deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_qua(ref_grid),
			   ref_grid_qua(original) ), "qua deep copy" );

  RSS( ref_geom_deep_copy( &ref_grid_geom(ref_grid),
			   ref_grid_geom(original) ), "geom deep copy" );

  ref_grid_twod(ref_grid) = ref_grid_twod(original);

  return REF_SUCCESS;
}

REF_STATUS ref_grid_empty_cell_clone( REF_GRID *ref_grid_ptr, REF_GRID parent )
{
  REF_GRID ref_grid;

  ref_malloc( *ref_grid_ptr, 1, REF_GRID_STRUCT );

  ref_grid = *ref_grid_ptr;

  ref_grid_node(ref_grid) = ref_grid_node(parent);

  RSS( ref_cell_create( &ref_grid_tet(ref_grid), 4, REF_FALSE ), "tet create" );
  RSS( ref_cell_create( &ref_grid_pyr(ref_grid), 5, REF_FALSE ), "pyr create" );
  RSS( ref_cell_create( &ref_grid_pri(ref_grid), 6, REF_FALSE ), "pri create" );
  RSS( ref_cell_create( &ref_grid_hex(ref_grid), 8, REF_FALSE ), "hex create" );

  ref_grid_cell(ref_grid,4) = NULL;

  RSS( ref_cell_create( &ref_grid_edg(ref_grid), 2, REF_TRUE ), "edg create" );
  RSS( ref_cell_create( &ref_grid_tri(ref_grid), 3, REF_TRUE ), "tri create" );
  RSS( ref_cell_create( &ref_grid_qua(ref_grid), 4, REF_TRUE ), "qua create" );

  RSS( ref_geom_create( &ref_grid_geom(ref_grid) ), "geom create" );

  return REF_SUCCESS;
}

REF_STATUS ref_grid_free( REF_GRID ref_grid )
{
  if ( NULL == (void *)ref_grid ) return REF_NULL;

  RSS( ref_geom_free( ref_grid_geom(ref_grid) ), "geom free");
  
  RSS( ref_cell_free( ref_grid_qua(ref_grid) ), "qua free");
  RSS( ref_cell_free( ref_grid_tri(ref_grid) ), "tri free");
  RSS( ref_cell_free( ref_grid_edg(ref_grid) ), "edg free");

  RSS( ref_cell_free( ref_grid_hex(ref_grid) ), "hex free");
  RSS( ref_cell_free( ref_grid_pri(ref_grid) ), "pri free");
  RSS( ref_cell_free( ref_grid_pyr(ref_grid) ), "pyr free");
  RSS( ref_cell_free( ref_grid_tet(ref_grid) ), "tet free");

  RSS( ref_node_free( ref_grid_node(ref_grid) ), "node free");

  ref_free( ref_grid );
  return REF_SUCCESS;
}

REF_STATUS ref_grid_free_cell_clone( REF_GRID ref_grid )
{
  if ( NULL == (void *)ref_grid ) return REF_NULL;

  RSS( ref_geom_free( ref_grid_geom(ref_grid) ), "geom free");
  
  RSS( ref_cell_free( ref_grid_qua(ref_grid) ), "qua free");
  RSS( ref_cell_free( ref_grid_tri(ref_grid) ), "tri free");
  RSS( ref_cell_free( ref_grid_edg(ref_grid) ), "edg free");

  RSS( ref_cell_free( ref_grid_hex(ref_grid) ), "hex free");
  RSS( ref_cell_free( ref_grid_pri(ref_grid) ), "pri free");
  RSS( ref_cell_free( ref_grid_pyr(ref_grid) ), "pyr free");
  RSS( ref_cell_free( ref_grid_tet(ref_grid) ), "tet free");

  ref_free( ref_grid );
  return REF_SUCCESS;
}

REF_STATUS ref_grid_inspect( REF_GRID ref_grid )
{
  printf("ref_grid = %p\n",(void *)ref_grid);
  printf(" %d node\n",ref_node_n(ref_grid_node(ref_grid)));
  printf(" %d tet\n",ref_cell_n(ref_grid_tet(ref_grid)));
  printf(" %d pyr\n",ref_cell_n(ref_grid_pyr(ref_grid)));
  printf(" %d pri\n",ref_cell_n(ref_grid_pri(ref_grid)));
  printf(" %d hex\n",ref_cell_n(ref_grid_hex(ref_grid)));
  printf(" %d edg\n",ref_cell_n(ref_grid_edg(ref_grid)));
  printf(" %d tri\n",ref_cell_n(ref_grid_tri(ref_grid)));
  printf(" %d qua\n",ref_cell_n(ref_grid_qua(ref_grid)));
  printf(" %d geom\n",ref_cell_n(ref_grid_geom(ref_grid)));

  return REF_SUCCESS;
}

REF_STATUS ref_grid_cell_with( REF_GRID ref_grid, REF_INT node_per,
			       REF_CELL *ref_cell )
{

  *ref_cell = NULL;

  switch ( node_per )
    {
    case 4:
      *ref_cell = ref_grid_tet(ref_grid);
      break;
    case 5:
      *ref_cell = ref_grid_pyr(ref_grid);
      break;
    case 6:
      *ref_cell = ref_grid_pri(ref_grid);
      break;
    case 8:
      *ref_cell = ref_grid_hex(ref_grid);
      break;
    default:
      printf("node_per %d\n",node_per);
      RSS(REF_FAILURE, "unexpected node_per");
      break;    
    }

  return REF_SUCCESS;
}

REF_STATUS ref_grid_face_with( REF_GRID ref_grid, REF_INT node_per,
			       REF_CELL *ref_cell )
{

  *ref_cell = NULL;

  switch ( node_per )
    {
    case 3:
      *ref_cell = ref_grid_tri(ref_grid);
      break;
    case 4:
      *ref_cell = ref_grid_qua(ref_grid);
      break;
    default:
      printf("node_per %d\n",node_per);
      RSS(REF_FAILURE, "unexpected node_per");
      break;    
    }

  return REF_SUCCESS;
}

REF_STATUS ref_grid_cell_has_face( REF_GRID ref_grid, 
				   REF_INT *face_nodes,
				   REF_BOOL *has_face )
{
  REF_CELL ref_cell;
  REF_INT group;
  REF_INT cell0, cell1;

  *has_face = REF_FALSE;
  
  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      RSS( ref_cell_with_face( ref_cell, face_nodes, &cell0, &cell1),
	   "face search failed" );
      *has_face = ( REF_EMPTY != cell0 );
      if ( *has_face ) return REF_SUCCESS;
    }

 return REF_SUCCESS;
}

REF_STATUS ref_grid_boundary_nodes( REF_GRID ref_grid, 
				    REF_INT boundary_tag, 
				    REF_INT *nnode, REF_INT *nface, 
				    REF_INT **g2l, REF_INT **l2g )
{
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  ref_node = ref_grid_node(ref_grid);

  ref_malloc_init( *g2l, ref_node_max(ref_node), REF_INT, REF_EMPTY );

  (*nnode) = 0;
  (*nface) = 0;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( boundary_tag == nodes[3] )
      {
	(*nface)++;
	for ( node = 0; node < 3; node++ )
	  if ( REF_EMPTY == (*g2l)[nodes[node]] )
	    { (*g2l)[nodes[node]] = (*nnode); (*nnode)++; }
      }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( boundary_tag == nodes[4] )
      {
	(*nface)++;
	for ( node = 0; node < 4; node++ )
	  if ( REF_EMPTY == (*g2l)[nodes[node]] )
	    { (*g2l)[nodes[node]] = (*nnode); (*nnode)++; }
      }

  ref_malloc( *l2g, *nnode, REF_INT );

  (*nnode) = 0;
  for ( node = 0 ; node < ref_node_max(ref_node) ; node++ )
    if ( REF_EMPTY != (*g2l)[node] ) 
      {
	(*g2l)[node] = (*nnode);
	(*l2g)[(*nnode)] = node;
	(*nnode)++;
      }

  return REF_SUCCESS;

}

REF_STATUS ref_grid_replace_node( REF_GRID ref_grid, 
				  REF_INT old_node, REF_INT new_node )
{
  REF_CELL ref_cell;
  REF_INT group;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      RSS( ref_cell_replace_node( ref_cell, old_node, new_node ),
	   "cell node replace" );
    }

  ref_cell = ref_grid_tri(ref_grid);
  RSS( ref_cell_replace_node( ref_cell, old_node, new_node ),
       "cell node replace" );

  ref_cell = ref_grid_qua(ref_grid);
  RSS( ref_cell_replace_node( ref_cell, old_node, new_node ),
       "cell node replace" );

  return REF_SUCCESS;
}

REF_STATUS ref_update_guess( REF_CELL ref_cell, REF_INT node0, REF_INT node1,
			     REF_INT *guess)
{
  REF_INT ncell, max_cell = 2, cell_list[2];

  RSS( ref_cell_list_with2( ref_cell, node0, node1,
			    max_cell, &ncell, cell_list ), "next" );
  if ( ncell == 1 )
    THROW("bary update hit boundary");

  if ( *guess == cell_list[0] )
    {
      *guess = cell_list[1];
      return REF_SUCCESS;
    }
  if ( *guess == cell_list[1] )
    {
      *guess = cell_list[0];
      return REF_SUCCESS;
    }

  return REF_NOT_FOUND;
}

REF_STATUS ref_grid_enclosing_tri( REF_GRID ref_grid, REF_DBL *xyz,
                                  REF_INT *tri, REF_DBL *bary )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT guess;
  REF_DBL inside = -1.0e-12;
  REF_INT step, limit;
  
  guess = *tri;
  if ( ! ref_cell_valid(ref_cell,guess) )
    {
      each_ref_cell_valid_cell( ref_cell, guess)
	{
	  break;
	}
    }
  RAS( ref_cell_valid(ref_cell,guess), "unable to find start");

  limit = 1000; /* was 10e6^(1/3), required 108 for twod testcase  */

  for ( step=0; step < limit; step++)
    {
      RSS( ref_cell_nodes( ref_cell, guess, nodes), "cell" );
      RSS( ref_node_bary3( ref_node, nodes, xyz, bary ), "bary");

      if ( step > 990 )
	{
	  printf("step %d, tri %d, bary %f %f %f\n",
		 step,guess,bary[0],bary[1],bary[2]);
	}
      
      if ( bary[0] >= inside &&
	   bary[1] >= inside &&
	   bary[2] >= inside )
	{
	  *tri = guess;
	  return REF_SUCCESS;
	}
      
      if ( bary[0] < bary[1] && bary[0] < bary[2] )
	{
	  RSS( ref_update_guess( ref_cell, nodes[1], nodes[2], &guess ),
	       "update 1 2");
	  continue;
	}

      if ( bary[1] < bary[0] && bary[1] < bary[2] )
	{
	  RSS( ref_update_guess( ref_cell, nodes[2], nodes[0], &guess ),
	       "update 2 0");
	  continue;
	}

      if ( bary[2] < bary[0] && bary[2] < bary[1] )
	{
	  RSS( ref_update_guess( ref_cell, nodes[0], nodes[1], &guess ),
	       "update 0 1");
	  continue;
	}

      THROW("unable to find the next step");
    }

  THROW("max steps exceeded");

  return REF_SUCCESS;
}
