#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_layer.h"
#include "ref_mpi.h"

#include "ref_fixture.h"
#include "ref_grid.h"
#include  "ref_cell.h"
#include  "ref_node.h"
#include   "ref_adj.h"
#include   "ref_sort.h"
#include   "ref_matrix.h"
#include   "ref_list.h"
#include   "ref_geom.h"
#include    "ref_math.h"

#include  "ref_export.h"
#include   "ref_edge.h"
#include   "ref_metric.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  {  /* make layers */
    REF_LAYER ref_layer;

    RSS(ref_layer_create(&ref_layer),"create");

    REIS( 0, ref_layer_n(ref_layer), "check total layers");

    RSS(ref_layer_free(ref_layer),"layer");
  }

  {  /* add layers to tet fixture */
    REF_LAYER ref_layer;
    REF_GRID ref_grid;
    REF_INT faceid;
    
    RSS(ref_fixture_tet_brick_grid( &ref_grid ), "tet brick");
    RSS(ref_metric_unit_node( ref_grid_node(ref_grid) ),"unit metric");
    RSS(ref_layer_create( &ref_layer ),"create");

    faceid = 6; RSS(ref_layer_attach( ref_layer, ref_grid, faceid ),"attach");
    faceid = 4; RSS(ref_layer_attach( ref_layer, ref_grid, faceid ),"attach");
    faceid = 3; RSS(ref_layer_attach( ref_layer, ref_grid, faceid ),"attach");
    faceid = 5; RSS(ref_layer_attach( ref_layer, ref_grid, faceid ),"attach");
    RSS(ref_layer_puff( ref_layer, ref_grid ),"puff");
    RSS(ref_layer_insert( ref_layer, ref_grid ),"insert");

    if ( argc > 1 )
      {
	RSS( ref_export_by_extension( ref_layer_grid(ref_layer),
				      argv[1] ), "tec" );
      }
    if ( argc > 2 )
      {
	RSS( ref_export_by_extension( ref_grid, argv[2] ), "tec" );
      }

    RSS(ref_layer_free(ref_layer),"layer");
    RSS(ref_grid_free(ref_grid),"grid");
  }

  RSS( ref_mpi_stop(  ), "stop" );
  
  return 0;
}
