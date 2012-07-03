#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



#include "ref_grid.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_matrix.h"

#include "ref_cell.h"
#include "ref_sort.h"

#include "ref_stitch.h"

#include "ref_import.c"
#include "ref_export.c"

#include "ref_dict.h"
#include "ref_mpi.h"
#include "ref_edge.h"

int main( int argc, char *argv[] )
{

  if (argc == 4)
    {      
      REF_GRID ref_grid;
      REF_INT tri_boundary;
      REF_INT qua_boundary;

      printf("importing %s\n",argv[1]);
      RSS(ref_import_by_extension( &ref_grid, argv[1] ),"imp");
      printf("complete.\n");

      RSS(ref_grid_inspect( ref_grid ), "inspection");

      tri_boundary = atoi(argv[2]);
      qua_boundary = atoi(argv[3]);

      RSS(ref_stitch_together( ref_grid, tri_boundary, qua_boundary ),"stitch");

      RSS(ref_grid_inspect( ref_grid ), "inspection");

      RSS(ref_export_by_extension( ref_grid, "ref_stitch_test.b8.ugrid" ),"out");
      RSS(ref_export_by_extension( ref_grid, "ref_stitch_test.tec" ),"out");

      RSS( ref_grid_free( ref_grid ), "free" );
      printf("done.\n");

      return 0;
    }

#define add_next_node( ref_grid, x, y, z )				\
  {REF_INT global,local;						\
    REF_NODE ref_node = ref_grid_node(ref_grid);			\
    RSS( ref_node_next_global(ref_node,&global), "get glob");		\
    RSS(ref_node_add(ref_node,global,&local),"add node");		\
    ref_node_xyz(ref_node,0,local) = (x);				\
    ref_node_xyz(ref_node,1,local) = (y);				\
    ref_node_xyz(ref_node,2,local) = (z);				\
  }

  {
    REF_GRID ref_grid;
    RSS(ref_grid_create(&ref_grid),"create");


    /*
                                 inode12
                                   /|\


                               inode11-----------inode10
                                 / \               /
                                /   \             /
                               /     \           /
                              /       \         /
                             /         \       /
                            /           \     /
                           /             \   /
                          /               \ /
                       inode8------------inode9
    */

    /*
                               inode7------------inode6
                                 /.                /|
                                / .               / |
                               /  .              /  |
                              /   .             /   |
                             /    .            /    |
                            /     .           /     |
                           /      .          /      |
                          /       .         /       |
                       inode4------------inode5     |
                         |      inode3.....|......inode2
                         |       .         |       /
                         |      .          |      /
     Z                   |     .           |     /
     ^   Y               |    .            |    /
     |  /                |   .             |   /
     | /                 |  .              |  /
     |/                  | .               | /
     +-----> X           |.                |/
                       inode0------------inode1

*/

    add_next_node( ref_grid, 0.0, 0.0, 0.0 );
    add_next_node( ref_grid, 0.0, 1.0, 0.0 );
    add_next_node( ref_grid, 1.0, 1.0, 0.0 );
    add_next_node( ref_grid, 0.0, 1.0, 0.0 );

    add_next_node( ref_grid, 0.0, 0.0, 1.0 );
    add_next_node( ref_grid, 0.0, 1.0, 1.0 );
    add_next_node( ref_grid, 1.0, 1.0, 1.0 );
    add_next_node( ref_grid, 0.0, 1.0, 1.0 );

    add_next_node( ref_grid, 0.0, 0.0, 1.0 );
    add_next_node( ref_grid, 0.0, 1.0, 1.0 );
    add_next_node( ref_grid, 1.0, 1.0, 1.0 );
    add_next_node( ref_grid, 0.0, 1.0, 1.0 );

    add_next_node( ref_grid, 0.5, 0.5, 2.0 );

    RSS(ref_grid_free(ref_grid),"create");
  }

  return 0;
}
