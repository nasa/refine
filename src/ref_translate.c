
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

#include "ref_export.h"
#include "ref_import.h"
#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_list.h"

#include "ref_math.h"
#include "ref_malloc.h"

static int print_usage(const char *name)
{
  printf("usage:\n");
  printf("  %s input_grid.extension output_grid.extension\n",name);
  printf("     [--shift dx dy dz]\n");
  printf("     [--scale s]\n");
  printf("     [--rotate degrees]\n");
  printf("     [--drop-face faceid]\n");
  printf("     [--compact-faceids]\n");
  return 0;
}

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_INT node;
  REF_DBL dx, dy, dz, ds;
  REF_DBL x, z, rotate_deg, rotate_rad;
  REF_INT faceid, ndrop;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  char *endptr;
  REF_INT pos;
  REF_INT min_faceid, max_faceid, nfaceid;
  REF_INT *new_faceid;

  if ( 3 > argc )
    return(print_usage(argv[0]));

  RSS( ref_import_by_extension( &ref_grid, argv[1] ), "import" );

  pos = 3;
  while( pos < argc ) {
    if( strcmp(argv[pos],"--shift") == 0 ) {
      printf("%d: --shift\n",pos);
      if ( pos+4 > argc )
	return(print_usage(argv[0]));
      ref_node = ref_grid_node(ref_grid);
      pos++; dx = strtod( argv[pos], &endptr);RAS(argv[pos]!=endptr,"parse dx");
      pos++; dy = strtod( argv[pos], &endptr);RAS(argv[pos]!=endptr,"parse dy");
      pos++; dz = strtod( argv[pos], &endptr);RAS(argv[pos]!=endptr,"parse dz");
      printf("%f %f %f\n",dx,dy,dz);
      each_ref_node_valid_node( ref_node, node )
	{
	  ref_node_xyz(ref_node,0,node) += dx;
	  ref_node_xyz(ref_node,1,node) += dy;
	  ref_node_xyz(ref_node,2,node) += dz;
	}
    }
    if( strcmp(argv[pos],"--scale") == 0 ) {
      printf("%d: --scale\n",pos);
      if ( pos+2 > argc )
	return(print_usage(argv[0]));
      ref_node = ref_grid_node(ref_grid);
      pos++; ds = strtod( argv[pos], &endptr);RAS(argv[pos]!=endptr,"parse ds");
      printf("%f\n",ds);
      each_ref_node_valid_node( ref_node, node )
	{
	  ref_node_xyz(ref_node,0,node) *= ds;
	  ref_node_xyz(ref_node,1,node) *= ds;
	  ref_node_xyz(ref_node,2,node) *= ds;
	}
    }
    if( strcmp(argv[pos],"--rotate") == 0 ) {
      printf("%d: --rotate\n",pos);
      if ( pos+2 > argc )
	return(print_usage(argv[0]));
      ref_node = ref_grid_node(ref_grid);
      pos++;
      rotate_deg = strtod( argv[pos], &endptr);
      rotate_rad = ref_math_in_radians(rotate_deg);
      RAS(argv[pos]!=endptr,"parse degree");
      printf("%f deg %f radian\n",rotate_deg,rotate_rad);
      each_ref_node_valid_node( ref_node, node )
	{
	  x = ref_node_xyz(ref_node,0,node);
	  z = ref_node_xyz(ref_node,2,node);
	  ref_node_xyz(ref_node,0,node) = x*cos(rotate_rad) - z*sin(rotate_rad);
	  ref_node_xyz(ref_node,2,node) = x*sin(rotate_rad) + z*cos(rotate_rad);
	}
    }
    if( strcmp(argv[pos],"--drop-face") == 0 ) {
      printf("%d: --drop-face\n",pos);
      if ( pos+2 > argc )
	return(print_usage(argv[0]));
      ref_node = ref_grid_node(ref_grid);
      pos++;
      faceid = strtol( argv[pos], &endptr, 10);
      RAS(argv[pos]!=endptr,"parse faceid to drop");
      printf(" dropping faceid %d\n",faceid);
      ref_cell = ref_grid_tri( ref_grid );
      ndrop = 0;
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  if ( faceid == nodes[ref_cell_node_per(ref_cell)] )
	    {
	      ndrop++;
	      RSS( ref_cell_remove( ref_cell, cell ), "drop" );
	    }
	}
      printf("dropped %d triangles from face %d\n",ndrop,faceid);
      ref_cell = ref_grid_qua( ref_grid );
      ndrop = 0;
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  if ( faceid == nodes[ref_cell_node_per(ref_cell)] )
	    {
	      ndrop++;
	      RSS( ref_cell_remove( ref_cell, cell ), "drop" );
	    }
	}
      printf("dropped %d quadrilaterals from face %d\n",ndrop,faceid);
    }
    if( strcmp(argv[pos],"--compact-faceids") == 0 ) {
      printf("%d: --compact-faceids\n",pos);
      RSS( ref_export_faceid_range( ref_grid, 
				    &min_faceid, &max_faceid ),
	   "min max faceid");
      printf(" faceid range %d %d\n", min_faceid, max_faceid);
      ref_malloc_init( new_faceid, max_faceid-min_faceid+1, REF_INT, REF_EMPTY);
      ref_cell = ref_grid_tri( ref_grid );
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  faceid = nodes[ref_cell_node_per(ref_cell)];
	  new_faceid[faceid-min_faceid] = 1;
	}
      ref_cell = ref_grid_qua( ref_grid );
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  faceid = nodes[ref_cell_node_per(ref_cell)];
	  new_faceid[faceid-min_faceid] = 1;
	}
      nfaceid=0;
      for ( faceid = 0; faceid < (max_faceid-min_faceid+1); faceid++)
	{
	  if ( new_faceid[faceid] > 0 )
	    {
	      nfaceid++;
	      new_faceid[faceid] = nfaceid;
	    }
	  printf("faceid %d is now %d\n",
		 faceid+min_faceid, new_faceid[faceid]);
	}
      printf("%d unique faceids detected\n",nfaceid);
      ref_cell = ref_grid_tri( ref_grid );
      each_ref_cell_valid_cell( ref_cell, cell )
	{
	  faceid = ref_cell_c2n( ref_cell, ref_cell_node_per(ref_cell), cell );
	  ref_cell_c2n( ref_cell, ref_cell_node_per(ref_cell), cell ) =
	    new_faceid[faceid - min_faceid];
	}
      ref_cell = ref_grid_qua( ref_grid );
      each_ref_cell_valid_cell( ref_cell, cell )
	{
	  faceid = ref_cell_c2n( ref_cell, ref_cell_node_per(ref_cell), cell );
	  ref_cell_c2n( ref_cell, ref_cell_node_per(ref_cell), cell ) =
	    new_faceid[faceid - min_faceid];
	}
    }
    pos++; 
  }
  RSS( ref_export_by_extension( ref_grid, argv[2] ), "export" );

  RSS(ref_grid_free(ref_grid),"free");

  return 0;
}
