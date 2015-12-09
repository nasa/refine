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
#include "ref_node.h"
#include "ref_sort.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_dict.h"
#include "ref_list.h"

#include "ref_edge.h"

#include "ref_args.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_DICT faceids;
  REF_INT arg, faceid, nlayers;
  REF_DBL first_thickness, total_thickness, mach;
  REF_DBL rate, total;
  REF_INT layer;
  REF_DBL thickness, xshift, mach_angle_rad;
  REF_INT aoa_pos;
  REF_DBL alpha_deg = 0;
  REF_DBL alpha_rad = 0;
  REF_INT origin_pos;
  REF_INT rotate_pos;
  REF_DBL rotate_deg = 0;
  REF_DBL rotate_rad = 0;
  REF_INT node;
  REF_DBL x, z;
  REF_BOOL extrude_radially = REF_FALSE;
  REF_DBL origin[3];
  REF_INT last_face_arg;

  if ( 7 > argc )
    {
      printf("usage: \n %s input.grid nlayers first_thickness total_thickness mach faceid [faceid...] [--aoa angle_of_attack_in_degrees] [--rotate angle_in_degrees]\n",
	     argv[0]);
      printf("  when first_thickness <= 0, it is set to a uniform grid,\n");
      printf("    first_thickness = total_thickness/nlayers\n");
      printf("  when nlayers < 0, extrude radially\n");
      printf("    (--aoa option only available for radial extrusion)\n");
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

  last_face_arg = argc;

  aoa_pos = REF_EMPTY;
  RXS( ref_args_find( argc, argv, "--aoa", &aoa_pos ),
       REF_NOT_FOUND, "aoa search" );

  if ( REF_EMPTY != aoa_pos )
    {
      if (aoa_pos >= argc)
	THROW("--aoa requires a value");
      if ( !extrude_radially )
	THROW("--aoa requires radial extrusion, nlayers < 0");
      alpha_deg = atof(argv[aoa_pos+1]);
      alpha_rad = ref_math_in_radians(alpha_deg);
      printf(" --aoa %f deg\n",alpha_deg);
      last_face_arg = MIN( last_face_arg, aoa_pos );
    }

  origin_pos = REF_EMPTY;
  RXS( ref_args_find( argc, argv, "--origin", &origin_pos ),
       REF_NOT_FOUND, "origin search" );

  if ( REF_EMPTY != origin_pos )
    {
      if (origin_pos >= argc-3)
	THROW("--origin requires a value");
      origin[0] = atof(argv[origin_pos+1]);
      origin[1] = atof(argv[origin_pos+2]);
      origin[2] = atof(argv[origin_pos+3]);
      printf(" --origin %f %f %f\n",origin[0],origin[1],origin[2]);
      last_face_arg = MIN( last_face_arg, origin_pos );
    }


  rotate_pos = REF_EMPTY;
  RXS( ref_args_find( argc, argv, "--rotate", &rotate_pos ),
       REF_NOT_FOUND, "rotate search" );

  if ( REF_EMPTY != rotate_pos )
    {
      if (rotate_pos >= argc)
	THROW("--rotate requires a value");
      rotate_deg = atof(argv[rotate_pos+1]);
      rotate_rad = ref_math_in_radians(rotate_deg);
      printf(" --rotate %f deg (%f rad)\n",rotate_deg,rotate_rad);
      last_face_arg = MIN( last_face_arg, rotate_pos );

      ref_node = ref_grid_node(ref_grid);
      each_ref_node_valid_node( ref_node, node )
	{
	  x = ref_node_xyz(ref_node,0,node);
	  z = ref_node_xyz(ref_node,2,node);
	  ref_node_xyz(ref_node,0,node) = x*cos(rotate_rad) - z*sin(rotate_rad);
	  ref_node_xyz(ref_node,2,node) = x*sin(rotate_rad) + z*cos(rotate_rad);
	}
    }

  printf("faceids");
  RSS( ref_dict_create( &faceids ), "create" );
  for( arg=6;arg<last_face_arg;arg++)
    {
      faceid = atoi( argv[arg] );
      RSS( ref_dict_store( faceids, faceid, REF_EMPTY ), "store" );
      printf(" %d",faceid);
    }
  printf("\n");

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

  if ( REF_EMPTY == origin_pos )
    RSS( ref_inflate_origin( ref_grid, faceids, origin ), "orig" );

  total = 0.0;
  for( layer=0;layer<nlayers;layer++)
    {
      thickness = first_thickness * pow(rate,layer);
      total = total+thickness;
      xshift = thickness / tan(mach_angle_rad);
      if ( extrude_radially )
	{
	  RSS( ref_inflate_radially( ref_grid, faceids, 
				     origin, thickness, 
				     mach_angle_rad, alpha_rad ), 
	       "inflate" );
	} 
      else 
	{
	  RSS( ref_inflate_face( ref_grid, faceids, 
				 origin, thickness, xshift ), 
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
