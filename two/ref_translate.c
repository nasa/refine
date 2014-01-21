#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"
#include "ref_import.h"
#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_list.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_BOOL valid_inputs, shift;
  REF_INT node;
  REF_DBL dx, dy, dz;
  REF_NODE ref_node;

  valid_inputs = ( (3 == argc) || (7 == argc) );
  shift = (7 == argc);

  if ( !valid_inputs )
    {
      printf("usage: %s input_grid.extension output_grid.extension [--shift dx dy dz]\n",argv[0]);
      return 0;
    }

  RSS( ref_import_by_extension( &ref_grid, argv[1] ), "import" );
  if ( shift )
    {
      ref_node = ref_grid_node(ref_grid);
      dx = atof( argv[4] );
      dy = atof( argv[5] );
      dz = atof( argv[6] );
      each_ref_node_valid_node( ref_node, node )
	{
	  ref_node_xyz(ref_node,0,node) += dx;
	  ref_node_xyz(ref_node,1,node) += dy;
	  ref_node_xyz(ref_node,2,node) += dz;
	}
    }
  RSS( ref_export_by_extension( ref_grid, argv[2] ), "export" );

  RSS(ref_grid_free(ref_grid),"free");

  return 0;
}
