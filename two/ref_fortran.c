
#include <stdlib.h>
#include <stdio.h>

#include "ref_fortran.h"
#include "ref_dict.h"
#include "ref_grid.h"
#include "ref_grid_export.h"

static REF_GRID ref_grid = NULL;
static REF_DICT ref_bc_flag = NULL;

REF_STATUS ref_init_node__(REF_INT *nnodes, REF_INT *nnodesg,
			   REF_INT *l2g, REF_INT *part, REF_INT *partition, 
			   REF_DBL *x, REF_DBL *y, REF_DBL *z )
{
  return ref_init_node_( nnodes, nnodesg, l2g, part, partition, x, y, z );
}
REF_STATUS ref_init_node_(REF_INT *nnodes, REF_INT *nnodesg,
			  REF_INT *l2g, REF_INT *part, REF_INT *partition, 
			  REF_DBL *x, REF_DBL *y, REF_DBL *z )
{
  REF_NODE ref_node;
  int node, pos;
  RSS( ref_grid_create( &ref_grid ), "create grid" );
  RSS( ref_dict_create( &ref_bc_flag ), "create bc_flag dict" );
  ref_node = ref_grid_node(ref_grid);

  ref_node_partition(ref_node) = *partition;
  ref_node_n_global(ref_node) = *nnodesg;

  for (node=0;node<(*nnodes);node++)
    {
      RSS( ref_node_add( ref_node, l2g[node]-1, &pos ), 
	   "add node");
      ref_node_xyz(ref_node,0,pos) = x[node];
      ref_node_xyz(ref_node,1,pos) = y[node];
      ref_node_xyz(ref_node,2,pos) = z[node];
      ref_node_part(ref_node,pos) = part[node]-1;
    }


  return REF_SUCCESS;
}

REF_STATUS ref_import_cell__(REF_INT *node_per_cell, REF_INT *ncell,
			     REF_INT *c2n)
{
  return ref_import_cell_(node_per_cell, ncell, c2n);
}
REF_STATUS ref_import_cell_(REF_INT *node_per_cell, REF_INT *ncell,
			     REF_INT *c2n)
{
  REF_INT *nodes;
  REF_CELL ref_cell;
  REF_INT cell, node, new_cell;
  switch ( *node_per_cell )
    {
    case 4:
      ref_cell = ref_grid_tet(ref_grid);
      break;
    case 5:
      ref_cell = ref_grid_pyr(ref_grid);
      break;
    case 6:
      ref_cell = ref_grid_pri(ref_grid);
      break;
    case 8:
      ref_cell = ref_grid_hex(ref_grid);
      break;
    default:
      RSS(REF_IMPLEMENT, "unexpected node_per_cell");
      break;    
    }
  nodes = (REF_INT *)malloc( (*node_per_cell) * sizeof(REF_INT));
  RNS( nodes,"malloc nodes NULL");
  for ( cell = 0 ; cell < (*ncell) ; cell++ ) 
    {
      for ( node = 0 ; node < (*node_per_cell) ; node++ )
	nodes[node] = c2n[node+(*node_per_cell)*cell] - 1;
      RSS( ref_cell_add( ref_cell, nodes, &new_cell ), "add cell");
    }
  free(nodes);
  return REF_SUCCESS;
}

REF_STATUS ref_import_boundary_flag__(REF_INT *boundary_index, 
				      REF_INT *boundary_flag)
{
  return ref_import_boundary_flag_(boundary_index, 
				   boundary_flag);
}
REF_STATUS ref_import_boundary_flag_(REF_INT *boundary_index, 
				     REF_INT *boundary_flag)
{
  return ref_dict_store( ref_bc_flag, *boundary_index, *boundary_flag );
}

REF_STATUS ref_viz__( void )
{
  return ref_viz_( );
}
REF_STATUS ref_viz_( void )
{
  char filename[1024];
  sprintf(filename,"ref_viz%04d.vtk",
	  ref_node_partition(ref_grid_node(ref_grid)));
  RSS( ref_grid_export_vtk( ref_grid, filename ), "export vtk");
  sprintf(filename,"ref_viz%04d.tec",
	  ref_node_partition(ref_grid_node(ref_grid)));
  RSS( ref_grid_export_tec( ref_grid, filename ), "export tec");
  return REF_SUCCESS;
}

REF_STATUS ref_free__( void )
{
  return ref_free_( );
}
REF_STATUS ref_free_( void )
{
  RSS( ref_dict_free( ref_bc_flag ), "free bc_flag");
  ref_bc_flag = NULL;
  RSS( ref_grid_free( ref_grid ), "free grid");
  ref_grid = NULL;
  return REF_SUCCESS;
}

