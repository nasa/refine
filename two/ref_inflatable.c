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
  REF_GRID ref_grid;
  REF_DICT faceids;
  REF_INT arg, faceid, nlayers;
  REF_DBL first_thickness, total_thickness, mach;
  REF_DBL rate, total;
  REF_INT layer;
  REF_DBL thickness, xshift, mach_angle_rad;
  REF_BOOL extrude_radially = REF_FALSE;

  if ( 7 > argc )
    {
      printf("usage: \n %s input.grid nlayers first_thickness total_thickness mach faceid [faceid...]\n",
	     argv[0]);
      printf("  when first_thickness <= 0, it is set to a uniform grid,\n");
      printf("    first_thickness = total_thickness/nlayers\n");
      printf("  when nlayers < 0, extrude radially\n");
      return 1;
    }


  RSS( ref_import_by_extension( &ref_grid, argv[1] ), "read grid" );

  nlayers = atoi( argv[2] );
  if ( nlayers < 0 )
    {
      nlayers = ABS(nlayers);
      extrude_radially = REF_TRUE;
    }
  first_thickness = atof( argv[3] );
  total_thickness = atof( argv[4] );
  mach = atof( argv[5] );

  RSS( ref_dict_create( &faceids ), "create" );
  for( arg=6;arg<argc;arg++)
    {
      faceid = atoi( argv[arg] );
      RSS( ref_dict_store( faceids, faceid, REF_EMPTY ), "store" );
    }

  if ( first_thickness <= 0.0 )
    {
      first_thickness = total_thickness/(REF_DBL)nlayers;
      rate = 1.0;
    }
  else
    {
      RSS( ref_inflate_rate( nlayers, first_thickness, total_thickness,
			     &rate ), "compute rate" );
    }

  mach_angle_rad = asin(1/mach);

  printf("inflating %d faces\n",ref_dict_n(faceids));
  printf("mach %f mach angle %f rad %f deg\n",
	 mach, mach_angle_rad,ref_math_in_degrees(mach_angle_rad));
  printf("first thickness %f\n", first_thickness);
  printf("total thickness %f\n", total_thickness);
  printf("rate %f\n", rate);

  total = 0.0;
  for( layer=0;layer<nlayers;layer++)
    {
      thickness = first_thickness * pow(rate,layer);
      total = total+thickness;
      xshift = thickness / tan(mach_angle_rad);
      if ( extrude_radially )
	{
	  RSS( ref_inflate_radially( ref_grid, faceids, thickness, xshift ), 
	       "inflate" );
	} 
      else 
	{
	  RSS( ref_inflate_face( ref_grid, faceids, thickness, xshift ), 
	       "inflate" );
	}
      printf("layer%5d of%5d : thickness %15.8e total %15.8e :%9d nodes\n",
	     layer+1,nlayers,thickness, total,
	     ref_node_n(ref_grid_node(ref_grid)));
    }

  RSS( ref_export_tec_surf( ref_grid, "inflated_boundary.tec" ), "tec" );
  RSS( ref_export_by_extension( ref_grid, "inflated.b8.ugrid" ), "b8" );

  RSS(ref_dict_free( faceids ), "free" );
  RSS(ref_grid_free(ref_grid),"free");

  return 0;
}
