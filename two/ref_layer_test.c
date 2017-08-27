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

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  {  /* make layers */
    REF_LAYER ref_layer;

    RSS(ref_layer_create(&ref_layer),"create");

    REIS( 0, ref_layer_n(ref_layer), "check total layers");

    RSS(ref_layer_free(ref_layer),"layer");
  }

  {  /* attach layers */
    REF_LAYER ref_layer;
    REF_GRID ref_grid;
    REF_DICT ref_dict;

    RSS(ref_fixture_tet_brick_grid( &ref_grid ), "tet brick");
    RSS(ref_dict_create( &ref_dict ), "dict");
    RSS(ref_layer_create( &ref_layer ),"create");
    RSS(ref_dict_store(ref_dict,6,REF_EMPTY),"mark top");

    RSS(ref_layer_attach( ref_layer, ref_grid, ref_dict ),"attach");
    
    RSS(ref_layer_free(ref_layer),"layer");
    RSS(ref_dict_free(ref_dict),"dict");
    RSS(ref_grid_free(ref_grid),"grid");
  }

  RSS( ref_mpi_stop(  ), "stop" );
  
  return 0;
}
