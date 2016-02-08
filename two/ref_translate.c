#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"
#include "ref_import.h"
#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_list.h"

static int print_usage(char *name)
{
  printf("usage:\n");
  printf("  %s input_grid.extension output_grid.extension\n",name);
  printf("     [--shift dx dy dz]\n");
  printf("     [--scale s]\n");
  return 0;
}

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_INT node;
  REF_DBL dx, dy, dz, ds;
  REF_NODE ref_node;
  char *endptr;
  REF_INT pos;

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
      printf("%f %f %f\n",dx,dy,dy);
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
    pos++; 
  }
  RSS( ref_export_by_extension( ref_grid, argv[2] ), "export" );

  RSS(ref_grid_free(ref_grid),"free");

  return 0;
}
