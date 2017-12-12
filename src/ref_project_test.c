
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_matrix.h"

#include "ref_cell.h"
#include "ref_sort.h"


#include "ref_project.h"

#include "ref_import.h"
#include "ref_export.h"

#include "ref_dict.h"
#include "ref_mpi.h"
#include "ref_edge.h"

#include "ref_math.h"

#include "ref_subdiv.h"
#include "ref_metric.h"

#include "ref_malloc.h"

int main( int argc, char *argv[] )
{
  REF_MPI ref_mpi;
  RSS( ref_mpi_create( &ref_mpi ), "create" );
 
  if (argc==2) 
    {
      REF_GRID ref_grid;
      REF_SUBDIV ref_subdiv;
      REF_EDGE ref_edge;
      REF_DBL *metric;
      REF_INT edge;

      printf("import from %s\n",argv[1]);
      RSS( ref_import_by_extension( &ref_grid, ref_mpi,
				    argv[1] ), "examine header" );
      printf(" complete.\n");

      printf("imply metric\n");
      ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
      RSS( ref_metric_imply_from( metric, ref_grid ), "imply" );
      RSS( ref_metric_to_node( metric, ref_grid_node(ref_grid) ), "to n" );
      ref_free( metric );
      printf(" complete.\n");

      printf("export\n");
      RSS(ref_export_tec_surf( ref_grid, "ref_project_test_orig.tec" ),"e");
      RSS(ref_export_html( ref_grid, "ref_project_test_orig.html" ),"e");
      printf(" complete.\n");

      printf("subdivide\n");
      RSS( ref_subdiv_create( &ref_subdiv, ref_grid ), "subdiv c");
      RSS( ref_subdiv_mark_all( ref_subdiv ), "mark all");
      RSS( ref_subdiv_split( ref_subdiv ), "split");
      printf(" complete.\n");

      printf("project.\n");
      ref_edge = ref_subdiv_edge( ref_subdiv );
      for (edge = 0 ; 
	   edge < ref_edge_n( ref_edge ) ; 
	   edge++ )
	{
	  RSS( ref_project_edge( ref_grid,
				 ref_edge_e2n( ref_edge, 0, edge ),
				 ref_edge_e2n( ref_edge, 1, edge ),
				 ref_subdiv_node( ref_subdiv, edge ) ), "proj");
	}
      printf(" complete.\n");

      printf("export\n");
      RSS(ref_export_tec_surf( ref_grid, "ref_project_test_embed.tec"),"e");
      RSS(ref_export_html( ref_grid, "ref_project_test_embed.html"),"e");
      printf(" complete.\n");

      RSS( ref_subdiv_free( ref_subdiv), "free" );
      RSS( ref_grid_free( ref_grid ), "free" );
      printf("done.\n");

      RSS( ref_mpi_free( ref_mpi ), "free" );
      return 0;
    }

  {
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node, tri[4], cell;
    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");
    ref_node = ref_grid_node(ref_grid);
    RSS( ref_node_add( ref_node, 0, &node ), "node");
    ref_node_xyz(ref_node, 0, node ) = 0.0;
    ref_node_xyz(ref_node, 1, node ) = 0.0;
    ref_node_xyz(ref_node, 2, node ) = 0.0;
    RSS( ref_node_add( ref_node, 1, &node ), "node");
    ref_node_xyz(ref_node, 0, node ) = 0.5;
    ref_node_xyz(ref_node, 1, node ) = 0.0;
    ref_node_xyz(ref_node, 2, node ) = 0.0;
    RSS( ref_node_add( ref_node, 2, &node ), "node");
    ref_node_xyz(ref_node, 0, node ) = 1.0;
    ref_node_xyz(ref_node, 1, node ) = 0.0;
    ref_node_xyz(ref_node, 2, node ) = 0.0;
    RSS( ref_node_add( ref_node, 3, &node ), "node");
    ref_node_xyz(ref_node, 0, node ) = 2.0;
    ref_node_xyz(ref_node, 1, node ) = 1.0;
    ref_node_xyz(ref_node, 2, node ) = 0.0;
    RSS( ref_node_add( ref_node, 4, &node ), "node");
    ref_node_xyz(ref_node, 0, node ) = 1.0;
    ref_node_xyz(ref_node, 1, node ) = 0.0;
    ref_node_xyz(ref_node, 2, node ) =-1.0;
    RSS( ref_node_add( ref_node, 5, &node ), "node");
    ref_node_xyz(ref_node, 0, node ) = 1.0;
    ref_node_xyz(ref_node, 1, node ) = 0.0;
    ref_node_xyz(ref_node, 2, node ) = 1.0;

    tri[3]=10;
    tri[0]=0;tri[1]=1;tri[2]=4;
    RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");
    tri[0]=1;tri[1]=2;tri[2]=4;
    RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");
    tri[0]=2;tri[1]=3;tri[2]=4;
    RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");

    tri[3]=20;
    tri[0]=3;tri[1]=2;tri[2]=5;
    RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");
    tri[0]=2;tri[1]=1;tri[2]=5;
    RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");
    tri[0]=1;tri[1]=0;tri[2]=5;
    RSS( ref_cell_add(ref_grid_tri(ref_grid), tri, &cell),"tri");

    RSS( ref_project_edge( ref_grid, 0, 2, 1 ), "proj");

    RSS( ref_grid_free( ref_grid ), "free" );
  }

  RSS( ref_mpi_free( ref_mpi ), "free" );
  return 0;
}
