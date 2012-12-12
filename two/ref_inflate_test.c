#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_inflate.h"

#include "ref_math.h"

#include "ref_import.h"
#include "ref_export.h"

#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_sort.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_dict.h"
#include "ref_list.h"

#include "ref_edge.h"

int main( int argc, char *argv[] )
{

  if ( 2 == argc )
    {
      printf("usage: \n %s input.grid faceid nlayers total_thickness mach\n",
	     argv[0]);
    }

  if ( 6 == argc )
    {
      REF_GRID ref_grid;
      REF_INT faceid, nlayers;
      REF_DBL total_thickness, mach;
      REF_INT layer;
      REF_DBL thickness, xshift, mach_angle_rad;

      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "read grid" );

      faceid = atoi( argv[2] );
      nlayers = atoi( argv[3] );
      total_thickness = atof( argv[4] );
      mach = atof( argv[5] );

      thickness = total_thickness/(REF_DBL)nlayers;
      mach_angle_rad = asin(1/mach);
      xshift = thickness / tan(mach_angle_rad);

      printf("inflating face %d\n",faceid);
      printf("mach %f mach angle %f rad %f deg\n",
	     mach, mach_angle_rad,ref_math_in_degrees(mach_angle_rad));
      printf("%d layers of thickness %f xshift %f\n",
	     nlayers, thickness, xshift);
      printf("total thickness %f\n", total_thickness);

      for( layer=0;layer<nlayers;layer++)
	{
	  RSS( ref_inflate_normal( ref_grid, faceid, thickness, xshift ), 
	       "normals" );
	  RSS( ref_inflate_face( ref_grid, faceid, thickness, xshift ), 
	       "inflate" );
	  printf("layer %d of %d : %d nodes\n",
		 layer+1,nlayers,ref_node_n(ref_grid_node(ref_grid)));
	}

      RSS( ref_export_by_extension( ref_grid, "ref_inflate_test.tec" ), "tec" );
      RSS( ref_export_by_extension( ref_grid, "ref_inflate_test.b8.ugrid" ), "b8" );

      RSS(ref_grid_free(ref_grid),"free");
    }

  return 0;
}
