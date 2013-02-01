
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_project.h"

#include "ref_node.h"
#include "ref_cell.h"

#include "ref_math.h"

REF_STATUS ref_project_edge( REF_GRID ref_grid,
			     REF_INT node0, REF_INT node1,
			     REF_INT new_node )
{
  REF_CELL ref_cell;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT ncell;
  REF_INT cell_list[2];
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL norm[3], area;
  REF_DBL norm0[3], area0;
  REF_DBL norm1[3], area1;
  REF_DBL s[3], r0[3], r1[3], temp[3];
  REF_DBL len;
  REF_INT i;
  REF_INT faceid;
  REF_DBL volume, min_volume;

  ref_cell = ref_grid_tri(ref_grid);

  /* exclude interior edges */
  RSS( ref_cell_list_with( ref_cell, node0, new_node,
			   2, &ncell, cell_list), "more than two" );
  if ( 2 != ncell ) return REF_SUCCESS;

  faceid = ref_cell_c2n(ref_cell,ref_cell_node_per(ref_cell),
			cell_list[0]);

  if ( faceid == ref_cell_c2n(ref_cell,ref_cell_node_per(ref_cell),
			      cell_list[1]) )
    { /* inbetween faces */
    }
  else
    {
      /* replace area with angle weight */
      norm0[0] = 0.0; norm0[1] = 0.0; norm0[2] = 0.0; area0 = 0.0;
      each_ref_cell_having_node( ref_cell, node0, item, cell )
	{
	  RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
	  /* only use same faceid */
	  if ( faceid != nodes[ref_cell_node_per(ref_cell)] ) continue;
	  RSS( ref_node_tri_normal( ref_node, nodes, norm ), "norm" );
	  RSS( ref_node_tri_area( ref_node, nodes, &area ), "area" );
	  norm0[0] += area*norm[0];
	  norm0[1] += area*norm[1];
	  norm0[2] += area*norm[2];
	  area0 += area;
	}
      if ( !ref_math_divisible(norm0[0],area0) ||
	   !ref_math_divisible(norm0[1],area0) ||
	   !ref_math_divisible(norm0[2],area0) ) THROW("divide area0");
      norm0[0] /= area0;
      norm0[1] /= area0;
      norm0[2] /= area0;
      RSS( ref_math_normalize( norm0 ), "n0" );
      
      norm1[0] = 0.0; norm1[1] = 0.0; norm1[2] = 0.0; area1 = 0.0;
      each_ref_cell_having_node( ref_cell, node1, item, cell )
	{
	  RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
	  /* only use same faceid */
	  if ( faceid != nodes[ref_cell_node_per(ref_cell)] ) continue;
	  RSS( ref_node_tri_normal( ref_node, nodes, norm ), "norm" );
	  RSS( ref_node_tri_area( ref_node, nodes, &area ), "area" );
	  norm1[0] += area*norm[0];
	  norm1[1] += area*norm[1];
	  norm1[2] += area*norm[2];
	  area1 += area;
	}
      if ( !ref_math_divisible(norm1[0],area1) ||
	   !ref_math_divisible(norm1[1],area1) ||
	   !ref_math_divisible(norm1[2],area1) ) THROW("divide area1");
      norm1[0] /= area1;
      norm1[1] /= area1;
      norm1[2] /= area1;
      RSS( ref_math_normalize( norm1 ), "n1" );

      s[0] = ref_node_xyz(ref_node,0,node1) - ref_node_xyz(ref_node,0,node0);
      s[1] = ref_node_xyz(ref_node,1,node1) - ref_node_xyz(ref_node,1,node0);
      s[2] = ref_node_xyz(ref_node,2,node1) - ref_node_xyz(ref_node,2,node0);

      len = sqrt( ref_math_dot(s,s) );
      
      ref_math_cross_product(s,norm0,temp);
      ref_math_cross_product(norm0,temp,r0);
      RSS( ref_math_normalize( r0 ), "r0" );
      r0[0] *= len;
      r0[1] *= len;
      r0[2] *= len;

      ref_math_cross_product(s,norm1,temp);
      ref_math_cross_product(norm1,temp,r1);
      RSS( ref_math_normalize( r1 ), "r1" );
      r1[0] *= len;
      r1[1] *= len;
      r1[2] *= len;

      for (i=0;i<3;i++)
	ref_node_xyz(ref_node,i,new_node) = 
	  0.500 * ref_node_xyz(ref_node,i,node0) +
	  0.125 * r0[i] - 0.125 * r1[i] +
	  0.500 * ref_node_xyz(ref_node,i,node1);
    }

  ref_cell = ref_grid_tet(ref_grid);
  min_volume = 1.0;
  each_ref_cell_having_node( ref_cell, new_node, item, cell)
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
      RSS( ref_node_tet_vol( ref_node, nodes, &volume ), "vol" );
      min_volume = MIN( min_volume, volume);
    }
  if ( min_volume <= 0.0 )
    {
      printf("vol %e\n",min_volume);
    }

  return REF_SUCCESS;
}

