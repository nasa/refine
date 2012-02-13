#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_import.h"
#include "ref_export.h"

#include "ref_fixture.h"

#include "ref_adj.h"
#include "ref_sort.h"
#include "ref_dict.h"

#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_metric.h"

int main( int argc, char *argv[] )
{

  if (2 == argc) 
    {
      RSS( ref_import_examine_header( argv[1] ), "examine header" );
      return 0;
    }

  { /* export import .fgrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.fgrid";
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    TSS(ref_export_fgrid( export_grid, file ), "export" );
    TSS(ref_import_fgrid( &import_grid, file ), "import" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n(ref_grid_node(import_grid)), "node count" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n_global(ref_grid_node(import_grid)), "global node count" );
    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export import .ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.ugrid";
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    TSS(ref_export_ugrid( export_grid, file ), "export" );
    TSS(ref_import_ugrid( &import_grid, file ), "import" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n(ref_grid_node(import_grid)), "node count" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n_global(ref_grid_node(import_grid)), "global node count" );
    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export import .b8.ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.b8.ugrid";
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    TSS(ref_export_b8_ugrid( export_grid, file ), "export" );
    TSS(ref_import_b8_ugrid( &import_grid, file ), "import" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n(ref_grid_node(import_grid)), "node count" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n_global(ref_grid_node(import_grid)), "global node count" );
    TEIS( ref_cell_n(ref_grid_tri(export_grid)),
	  ref_cell_n(ref_grid_tri(import_grid)), "tri count" );
    TEIS( ref_cell_n(ref_grid_qua(export_grid)),
	  ref_cell_n(ref_grid_qua(import_grid)), "qua count" );
    TEIS( ref_cell_n(ref_grid_tet(export_grid)),
	  ref_cell_n(ref_grid_tet(import_grid)), "tet count" );
    TES( ref_node_xyz( ref_grid_node(export_grid),0,1),
	 ref_node_xyz( ref_grid_node(import_grid),0,1), "x 1" );
    TEIS( ref_cell_c2n(ref_grid_tet(export_grid),0,0),
	  ref_cell_c2n(ref_grid_tet(import_grid),0,0), "tet node0" );
    TEIS( ref_cell_c2n(ref_grid_tet(export_grid),1,0),
	  ref_cell_c2n(ref_grid_tet(import_grid),1,0), "tet node 1" );
    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export import .b8.ugrid tet by extension */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.b8.ugrid";
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    TSS(ref_export_b8_ugrid( export_grid, file ), "export" );
    TSS(ref_import_by_extension( &import_grid, file ), "import" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n(ref_grid_node(import_grid)), "node count" );
    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* import mapbc */
    REF_DICT ref_dict;
    char filename[] = "ref_import_test.mapbc";
    FILE *file;
    REF_INT key, value;

    file = fopen(filename,"w");
    RNS(file, "unable to open file" );

    fprintf(file, " 3\n");
    fprintf(file, " 1 4000\n");
    fprintf(file, " 2 5000 farfield_riem\n");
    fprintf(file, "   3 6662 \n");

    fclose(file);

    TSS(ref_import_mapbc( &ref_dict, filename ), "import" );

    key = 1;
    TSS(ref_dict_value(ref_dict,key,&value),"missing");
    TEIS(4000,value,"get value");

    key = 2;
    TSS(ref_dict_value(ref_dict,key,&value),"missing");
    TEIS(5000,value,"get value");

    key = 3;
    TSS(ref_dict_value(ref_dict,key,&value),"missing");
    TEIS(6662,value,"get value");

    TSS(ref_dict_free(ref_dict),"free");

    TEIS(0, remove( filename ), "test clean up");
  }

  return 0;
}
