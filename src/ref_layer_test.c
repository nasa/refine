
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
    RSS(ref_layer_recon( ref_layer, ref_grid ),"insert");

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
