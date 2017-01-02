
#include <stdlib.h>
#include <stdio.h>

#include "ref_clump.h"

#include "ref_dict.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_malloc.h"
#include "ref_adapt.h"
#include "ref_matrix.h"

static REF_STATUS ref_clump_zone_around( FILE *f,
					 REF_CELL ref_cell,
					 REF_DICT ref_dict,
					 char *zonetype,
					 REF_DICT node_dict,
					 REF_NODE ref_node,
					 REF_INT node )
{
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL xyz_comp[3], xyz_phys[3];
  REF_INT local;
  REF_DBL jacob[9];

  if ( ref_dict_n(ref_dict) <= 0 ) return REF_SUCCESS;

  RSS( ref_matrix_jacob_m( ref_node_metric_ptr(ref_node,node),
                           jacob ), "jac");

  fprintf(f,
	  "zone t=%s, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  zonetype, ref_dict_n(node_dict), ref_dict_n(ref_dict),
	  "point", zonetype );

  for ( item = 0; item < ref_dict_n(node_dict); item++ )
    {
      local = ref_dict_key(node_dict,item);
      xyz_phys[0] = ref_node_xyz(ref_node,0,local)
	- ref_node_xyz(ref_node,0,node);
      xyz_phys[1] = ref_node_xyz(ref_node,1,local)
	- ref_node_xyz(ref_node,1,node);
      xyz_phys[2] = ref_node_xyz(ref_node,2,local)
	- ref_node_xyz(ref_node,2,node);
      RSS( ref_matrix_vect_mult( jacob, xyz_phys, xyz_comp ), "ax");
      fprintf(f, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
	      ref_node_xyz(ref_node,0,local),
	      ref_node_xyz(ref_node,1,local),
	      ref_node_xyz(ref_node,2,local),	      
	      xyz_comp[0], xyz_comp[1], xyz_comp[2]);
    }
  
  for ( item = 0; item < ref_dict_n(ref_dict); item++ )
    {
      cell = ref_dict_key(ref_dict,item);
      RSS( ref_cell_nodes(ref_cell,cell,nodes), "n");
      each_ref_cell_cell_node(ref_cell,cell_node)
	{
	  RSS( ref_dict_location( node_dict,
				  nodes[cell_node], &local), "ret");
	  fprintf(f," %d",local + 1);
	}
      fprintf(f,"\n");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_clump_around( REF_GRID ref_grid, REF_INT node,
			     char *filename )
{
  REF_DICT node_dict, tri_dict, tet_dict;
  REF_DICT ref_dict;
  REF_CELL ref_cell;
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  char *zonetype;

  FILE *f;

  RSS(ref_dict_create(&node_dict),"create nodes");
  RSS(ref_dict_create(&tri_dict),"create tris");
  RSS(ref_dict_create(&tet_dict),"create tets");

  ref_cell = ref_grid_tri(ref_grid);
  ref_dict = tri_dict;
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes(ref_cell,cell,nodes), "n");
      RSS( ref_dict_store( ref_dict, cell, 0 ), "store");
      each_ref_cell_cell_node(ref_cell,cell_node)
	RSS( ref_dict_store( node_dict, nodes[cell_node], 0 ), "store");
    }
  ref_cell = ref_grid_tet(ref_grid);
  ref_dict = tet_dict;
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes(ref_cell,cell,nodes), "n");
      RSS( ref_dict_store( ref_dict, cell, 0 ), "store");
      each_ref_cell_cell_node(ref_cell,cell_node)
	RSS( ref_dict_store( node_dict, nodes[cell_node], 0 ), "store");
    }


  f = fopen(filename,"w");
  if (NULL == (void *)f)
    printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  fprintf(f, "title=\"tecplot refine clump file\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\" \"xm\" \"ym\" \"zm\"\n");

  ref_cell = ref_grid_tri(ref_grid);
  ref_dict = tri_dict;
  zonetype="fetriangle";
  RSS( ref_clump_zone_around( f, ref_cell, ref_dict, zonetype,
			      node_dict,
			      ref_grid_node(ref_grid), node ), "zone" );
  
  ref_cell = ref_grid_tet(ref_grid);
  ref_dict = tet_dict;
  zonetype = "fetetrahedron";
       RSS( ref_clump_zone_around( f, ref_cell, ref_dict, zonetype,
			      node_dict,
			      ref_grid_node(ref_grid), node ), "zone" );

  fclose(f);

  RSS(ref_dict_free(tet_dict),"free tet");
  RSS(ref_dict_free(tri_dict),"free tris");
  RSS(ref_dict_free(node_dict),"free nodes");

  return REF_SUCCESS;
}

REF_STATUS ref_clump_tri_around( REF_GRID ref_grid, REF_INT node,
				 char *filename )
{
  REF_DICT node_dict, tri_dict;
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT local;
  REF_DBL xyz_comp[3], xyz_phys[3];
  REF_DBL jacob[9];

  FILE *f;

  RSS(ref_dict_create(&node_dict),"create nodes");
  RSS(ref_dict_create(&tri_dict),"create tris");

  each_ref_cell_having_node( ref_cell, node, item, cell )
  {
    RSS( ref_cell_nodes(ref_cell,cell,nodes), "n");
    RSS( ref_dict_store( tri_dict, cell, 0 ), "store");
    each_ref_cell_cell_node(ref_cell,cell_node)
    RSS( ref_dict_store( node_dict, nodes[cell_node], 0 ), "store");
  }

  f = fopen(filename,"w");
  if (NULL == (void *)f)
    printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  fprintf(f, "title=\"tecplot refine clump file\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  fprintf(f,
          "zone t=clump, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
          ref_dict_n(node_dict), ref_dict_n(tri_dict),
          "point", "fequadrilateral" );

  RSS( ref_matrix_jacob_m( ref_node_metric_ptr(ref_grid_node(ref_grid),node),
                           jacob ), "jac");
  for ( item = 0; item < ref_dict_n(node_dict); item++ )
    {
      local = ref_dict_key(node_dict,item);
      xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid),0,local)
                    - ref_node_xyz(ref_grid_node(ref_grid),0,node);
      xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid),1,local)
                    - ref_node_xyz(ref_grid_node(ref_grid),1,node);
      xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid),2,local)
                    - ref_node_xyz(ref_grid_node(ref_grid),2,node);
      RSS( ref_matrix_vect_mult( jacob, xyz_phys, xyz_comp ), "ax");
      fprintf(f, " %.16e %.16e %.16e\n", xyz_comp[0], xyz_comp[1], xyz_comp[2]);
    }

  for ( item = 0; item < ref_dict_n(tri_dict); item++ )
    {
      cell = ref_dict_key(tri_dict,item);
      RSS( ref_cell_nodes(ref_cell,cell,nodes), "n");
      each_ref_cell_cell_node(ref_cell,cell_node)
      {
        RSS( ref_dict_location( node_dict,
                                nodes[cell_node], &local), "ret");
        fprintf(f," %d",local + 1);
      }
      RSS( ref_dict_location( node_dict,
                              nodes[0], &local), "ret");
      fprintf(f," %d",local + 1);
      fprintf(f,"\n");
    }

  fclose(f);

  RSS(ref_dict_free(tri_dict),"free tris");
  RSS(ref_dict_free(node_dict),"free nodes");

  return REF_SUCCESS;
}

REF_STATUS ref_clump_stuck_edges( REF_GRID ref_grid, REF_DBL ratio_tol )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT ntarget;
  REF_INT node, node0, node1;
  REF_INT edge;
  REF_DBL edge_ratio;

  char filename[1024];

  RSS( ref_edge_create( &ref_edge, ref_grid ), "orig edges" );

  ref_malloc_init( ratio, ref_node_max(ref_node),
                   REF_DBL, 2.0*ref_adapt_collapse_ratio );

  for (edge = 0; edge<ref_edge_n(ref_edge); edge++)
    {
      node0 = ref_edge_e2n( ref_edge, 0, edge );
      node1 = ref_edge_e2n( ref_edge, 1, edge );
      RSS( ref_node_ratio( ref_node, node0, node1,
                           &edge_ratio ), "ratio");
      ratio[node0] = MIN( ratio[node0], edge_ratio );
      ratio[node1] = MIN( ratio[node1], edge_ratio );
    }

  ntarget = 0;
  for ( node = 0; node < ref_node_max(ref_node); node++ )
    if ( ratio[node] < ratio_tol*ref_adapt_collapse_ratio )
      {
        sprintf(filename,"clump%d.t",ntarget);
        RSS(ref_clump_around(ref_grid, node, filename ), "dump");
        ntarget++;
      }

  ref_free( ratio )
  RSS( ref_edge_free( ref_edge ), "free edges" );

  return REF_SUCCESS;
}

REF_STATUS ref_clump_stuck_edges_twod( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT ntarget;
  REF_INT node, node0, node1;
  REF_INT edge;
  REF_DBL edge_ratio;
  REF_BOOL active;

  char filename[1024];

  RSS( ref_edge_create( &ref_edge, ref_grid ), "orig edges" );

  ref_malloc_init( ratio, ref_node_max(ref_node),
                   REF_DBL, 2.0*ref_adapt_collapse_ratio );

  for (edge = 0; edge<ref_edge_n(ref_edge); edge++)
    {
      node0 = ref_edge_e2n( ref_edge, 0, edge );
      node1 = ref_edge_e2n( ref_edge, 1, edge );
      RSS( ref_node_edge_twod( ref_node, node0, node1,
                               &active ), "act" );
      if ( !active )
        continue;

      RSS( ref_node_ratio( ref_node, node0, node1,
                           &edge_ratio ), "ratio");
      ratio[node0] = MIN( ratio[node0], edge_ratio );
      ratio[node1] = MIN( ratio[node1], edge_ratio );
    }

  ntarget = 0;
  for ( node = 0; node < ref_node_max(ref_node); node++ )
    if ( ratio[node] < ref_adapt_collapse_ratio )
      {
        sprintf(filename,"clump%d.t",ntarget);
        RSS(ref_clump_tri_around(ref_grid, node, filename ), "dump");
        ntarget++;
      }

  return REF_SUCCESS;
}
