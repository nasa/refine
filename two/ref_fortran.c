
#include <stdlib.h>
#include <stdio.h>

#include "ref_fortran.h"
#include "ref_metric.h"
#include "ref_grid.h"
#include "ref_grid_export.h"

static REF_GRID ref_grid = NULL;

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

REF_STATUS ref_import_boundary__(REF_INT *node_per_face, REF_INT *nface,
			     REF_INT *f2n, REF_INT *boundary_index )
{
  return ref_import_boundary_(node_per_face, nface, f2n, boundary_index);
}
REF_STATUS ref_import_boundary_(REF_INT *node_per_face, REF_INT *nface,
			     REF_INT *f2n, REF_INT *boundary_index )
{
  REF_INT *nodes;
  REF_CELL ref_cell;
  REF_INT face, node, new_face;
  switch ( *node_per_face )
    {
    case 3:
      ref_cell = ref_grid_tri(ref_grid);
      break;
    case 4:
      ref_cell = ref_grid_qua(ref_grid);
      break;
    default:
      RSS(REF_IMPLEMENT, "unexpected node_per_face");
      break;    
    }
  nodes = (REF_INT *)malloc( ((*node_per_face)+1) * sizeof(REF_INT));
  RNS( nodes,"malloc nodes NULL");
  for ( face = 0 ; face < (*nface) ; face++ ) 
    {
      for ( node = 0 ; node < (*node_per_face) ; node++ )
	nodes[node] = f2n[node+(*node_per_face)*face] - 1;
      nodes[(*node_per_face)] = (*boundary_index);
      RSS( ref_cell_add( ref_cell, nodes, &new_face ), "add face");
    }
  free(nodes);
  return REF_SUCCESS;
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

REF_STATUS ref_import_metric__(REF_INT *nnodes, REF_DBL *m )
{
  return ref_import_metric_( nnodes, m );
}
REF_STATUS ref_import_metric_(REF_INT *nnodes, REF_DBL *m )
{
  int node;
  REF_METRIC ref_metric;

  ref_metric = ref_grid_metric(ref_grid);

  for (node = 0; node < (*nnodes); node++)
    RSS(ref_metric_set(ref_metric,node,&(m[6*node])),"set metric");

  return REF_SUCCESS;
}

REF_STATUS ref_free__( void )
{
  return ref_free_( );
}
REF_STATUS ref_free_( void )
{
  RSS( ref_grid_free( ref_grid ), "free grid");
  ref_grid = NULL;
  return REF_SUCCESS;
}

