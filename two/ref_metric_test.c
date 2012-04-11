#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_metric.h"

#include "ref_grid.h"
#include "ref_node.h"
#include "ref_cell.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_sort.h"
#include "ref_matrix.h"

#include "ref_fixture.h"
#include "ref_malloc.h"

#include "ref_mpi.h"

#include "ref_part.h"
#include "ref_import.h"
#include "ref_migrate.h"
#include "ref_dict.h"

int main(  int argc, char *argv[] )
{

  if ( argc > 2 )
    {
      REF_GRID ref_grid;
      REF_DBL *metric_file;
      REF_DBL *metric_imply;
      REF_DBL *metric;
      REF_INT node, im;

      RSS( ref_mpi_start( argc, argv ), "start" );

      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "examine header" );
      ref_malloc( metric_file, 
		  6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
      ref_malloc( metric_imply, 
		  6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
      ref_malloc( metric, 
		  6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );

      RSS( ref_part_metric( ref_grid_node(ref_grid), argv[2] ), "get metric");

      each_ref_node_valid_node( ref_grid_node(ref_grid), node )
	for(im=0;im<6;im++)
	  metric_file[im+6*node] = 
	    ref_node_metric(ref_grid_node(ref_grid),im,node);

      RSS( ref_metric_imply_from( metric_imply, ref_grid ), "imply" );

      RSS( ref_metric_smr( metric_imply, metric_file, metric, ref_grid ), "smr" );

      ref_free( metric );
      ref_free( metric_imply );
      ref_free( metric_file );
      RSS( ref_grid_free( ref_grid ), "free");

      RSS( ref_mpi_stop( ), "stop" );
    }

  {  /* imply metric right tet */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS( ref_fixture_tet_grid( &ref_grid ), "tet" );

    ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );

    RSS( ref_metric_imply_from( metric, ref_grid ), "imply" );

    each_ref_node_valid_node( ref_grid_node(ref_grid), node )
      {
	RWDS( 1.0, metric[0+6*node], tol, "m[0]");
	RWDS( 0.5, metric[1+6*node], tol, "m[1]");
	RWDS( 0.5, metric[2+6*node], tol, "m[2]");
	RWDS( 1.0, metric[3+6*node], tol, "m[3]");
	RWDS( 0.5, metric[4+6*node], tol, "m[4]");
	RWDS( 1.0, metric[5+6*node], tol, "m[5]");
      }

    ref_free( metric );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  {  /* imply metric right prism */
    REF_DBL tol = 0.00001;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS( ref_fixture_pri_grid( &ref_grid ), "tet" );

    ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );

    RSS( ref_metric_imply_from( metric, ref_grid ), "imply" );

    each_ref_node_valid_node( ref_grid_node(ref_grid), node )
      {
	RWDS( 0.97569,  metric[0+6*node], tol, "m[0]");
	RWDS( 0.864269, metric[1+6*node], tol, "m[1]");
	RWDS(-0.43566,  metric[2+6*node], tol, "m[2]");
	RWDS( 1.714524, metric[3+6*node], tol, "m[3]");
	RWDS(-0.864269, metric[4+6*node], tol, "m[4]");
	RWDS( 0.97569,  metric[5+6*node], tol, "m[5]");
      }

    ref_free( metric );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  {  /* imply metric non tet prism */
    REF_DBL tol = 0.00001;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS( ref_fixture_pri_grid( &ref_grid ), "tet" );

    ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );

    RSS( ref_metric_imply_non_tet( metric, ref_grid ), "imply" );

    each_ref_node_valid_node( ref_grid_node(ref_grid), node )
      {
	RWDS( 0.97569,  metric[0+6*node], tol, "m[0]");
	RWDS( 0.864269, metric[1+6*node], tol, "m[1]");
	RWDS(-0.43566,  metric[2+6*node], tol, "m[2]");
	RWDS( 1.714524, metric[3+6*node], tol, "m[3]");
	RWDS(-0.864269, metric[4+6*node], tol, "m[4]");
	RWDS( 0.97569,  metric[5+6*node], tol, "m[5]");
      }

    ref_free( metric );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  {  /* imply metric tet prism */
    REF_DBL tol = 0.00001;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS( ref_fixture_pri_tet_cap_grid( &ref_grid ), "tet" );

    ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );

    RSS( ref_metric_imply_from( metric, ref_grid ), "imply" );

    node = 6;
    RWDS( 1.00, metric[0+6*node], tol, "m[0]");
    RWDS( 0.50, metric[1+6*node], tol, "m[1]");
    RWDS( 0.05, metric[2+6*node], tol, "m[2]");
    RWDS( 1.00, metric[3+6*node], tol, "m[3]");
    RWDS( 0.05, metric[4+6*node], tol, "m[4]");
    RWDS( 0.67, metric[5+6*node], tol, "m[5]");

    RSS( ref_metric_imply_non_tet( metric, ref_grid ), "imply" );

    node = 6;
    RWDS( 1.0,  metric[0+6*node], tol, "m[0]");
    RWDS( 0.5, metric[1+6*node], tol, "m[1]");
    RWDS( 0.05,  metric[2+6*node], tol, "m[2]");
    RWDS( 1.0, metric[3+6*node], tol, "m[3]");
    RWDS( 0.05, metric[4+6*node], tol, "m[4]");
    RWDS( 0.67,  metric[5+6*node], tol, "m[5]");

    ref_free( metric );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  return 0;
}
