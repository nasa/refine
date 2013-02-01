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

#include "ref_import.c"
#include "ref_export.c"

#include "ref_dict.h"
#include "ref_mpi.h"
#include "ref_edge.h"

#include "ref_math.h"

#include "ref_subdiv.h"
#include "ref_metric.h"

int main( int argc, char *argv[] )
{

  if (argc==2) 
    {
      REF_GRID ref_grid;
      REF_SUBDIV ref_subdiv;
      REF_EDGE ref_edge;
      REF_DBL *metric;
      REF_INT edge;

      printf("import from %s\n",argv[1]);
      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "examine header" );
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

      return 0;
    }

  return 0;
}
