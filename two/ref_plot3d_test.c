#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_plot3d.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_dict.h"
#include "ref_node.h"
#include "ref_cell.h"
#include "ref_adj.h"
#include "ref_mpi.h"
#include "ref_matrix.h"
#include "ref_sort.h"
#include "ref_list.h"

#include "ref_math.h"

#include "ref_malloc.h"

static REF_STATUS set_up_2x2( REF_PATCH *ref_patch_ptr )
{
  REF_PATCH ref_patch;
  ref_malloc( *ref_patch_ptr, 1, REF_PATCH_STRUCT );
  ref_patch = *ref_patch_ptr;
  ref_patch->idim = 2;
  ref_patch->jdim = 2;
  ref_patch->kdim = 1;
  ref_malloc( ref_patch->xyz, 3*ref_patch->idim*ref_patch->jdim, REF_DBL );

  ref_patch_xyz(ref_patch,0,0,0) = 0.0;
  ref_patch_xyz(ref_patch,1,0,0) = 0.0;
  ref_patch_xyz(ref_patch,2,0,0) = 0.0;

  ref_patch_xyz(ref_patch,0,1,0) = 1.0;
  ref_patch_xyz(ref_patch,1,1,0) = 0.0;
  ref_patch_xyz(ref_patch,2,1,0) = 0.0;

  ref_patch_xyz(ref_patch,0,0,1) = 0.0;
  ref_patch_xyz(ref_patch,1,0,1) = 1.0;
  ref_patch_xyz(ref_patch,2,0,1) = 0.0;

  ref_patch_xyz(ref_patch,0,1,1) = 1.0;
  ref_patch_xyz(ref_patch,1,1,1) = 1.0;
  ref_patch_xyz(ref_patch,2,1,1) = 0.0;

  return REF_SUCCESS;
}

int main( int argc, char *argv[] )
{

  if (2 == argc) 
    {
      REF_PLOT3D ref_plot3d;
      RSS( ref_plot3d_from_file( &ref_plot3d, argv[1] ), "from file" );
      RSS( ref_plot3d_tec( ref_plot3d, "ref_plot3d_test.tec" ), "tec" );
      RSS( ref_plot3d_free( ref_plot3d ), "free" );
      return 0;
    }

  /*

~/refine/strict/two/ref_plot3d_test \
 ~/cases/dpw5/dpw5mbgrids_rev01/L1.T.rev01.p3d \
 ~/cases/dpw5/mavriplis/hybrid/L1/L1.T.rev01.p3d.hybrid.b8.ugrid

  */

  if (3 == argc) 
    {
      REF_PLOT3D ref_plot3d;
      REF_GRID ref_grid;

      printf("plot3d %s\n",argv[1]);
      RSS( ref_plot3d_from_file( &ref_plot3d, argv[1] ), "from file" );

      printf("unstruct %s\n",argv[2]);
      RSS( ref_import_by_extension( &ref_grid, argv[2] ), "by ext" );

      printf("mate\n");
      RSS( ref_plot3d_mate( ref_plot3d, ref_grid ), "mate" );

      RSS( ref_plot3d_free( ref_plot3d ), "free" );
      RSS( ref_grid_free( ref_grid ), "free" );

      return 0;
    }

  {
    REF_PATCH ref_patch;
    REF_DBL xyz[3], uv[2];
    RSS( set_up_2x2( &ref_patch ), "2x2" );

    uv[0] = 0.0; uv[1] = 0.0;
    RSS( ref_patch_xyz_at( ref_patch, uv, xyz ), "at");
    RWDS( 0.0, xyz[0], -1.0, "xyz[0]");
    RWDS( 0.0, xyz[1], -1.0, "xyz[1]");
    RWDS( 0.0, xyz[2], -1.0, "xyz[2]");

    uv[0] = 1.0; uv[1] = 0.0;
    RSS( ref_patch_xyz_at( ref_patch, uv, xyz ), "at");
    RWDS( 1.0, xyz[0], -1.0, "xyz[0]");
    RWDS( 0.0, xyz[1], -1.0, "xyz[1]");
    RWDS( 0.0, xyz[2], -1.0, "xyz[2]");

    uv[0] = 0.0; uv[1] = 1.0;
    RSS( ref_patch_xyz_at( ref_patch, uv, xyz ), "at");
    RWDS( 0.0, xyz[0], -1.0, "xyz[0]");
    RWDS( 1.0, xyz[1], -1.0, "xyz[1]");
    RWDS( 0.0, xyz[2], -1.0, "xyz[2]");

    uv[0] = 1.0; uv[1] = 1.0;
    RSS( ref_patch_xyz_at( ref_patch, uv, xyz ), "at");
    RWDS( 1.0, xyz[0], -1.0, "xyz[0]");
    RWDS( 1.0, xyz[1], -1.0, "xyz[1]");
    RWDS( 0.0, xyz[2], -1.0, "xyz[2]");

    uv[0] = 0.5; uv[1] = 0.0;
    RSS( ref_patch_xyz_at( ref_patch, uv, xyz ), "at");
    RWDS( 0.5, xyz[0], -1.0, "xyz[0]");
    RWDS( 0.0, xyz[1], -1.0, "xyz[1]");
    RWDS( 0.0, xyz[2], -1.0, "xyz[2]");

    RSS( ref_patch_free( ref_patch ), "free" );
  }

  {
    REF_PATCH ref_patch;
    REF_DBL xyz[3], uv[2];
    RSS( set_up_2x2( &ref_patch ), "2x2" );

    xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
    RSS( ref_patch_locate( ref_patch, xyz, uv ), "loc");
    RWDS( 0.0, uv[0], -1.0, "uv[0]");
    RWDS( 0.0, uv[1], -1.0, "uv[1]");

    xyz[0] = 1.0; xyz[1] = 0.0; xyz[2] = 0.0;
    RSS( ref_patch_locate( ref_patch, xyz, uv ), "loc");
    RWDS( 1.0, uv[0], -1.0, "uv[0]");
    RWDS( 0.0, uv[1], -1.0, "uv[1]");

    xyz[0] = 0.0; xyz[1] = 1.0; xyz[2] = 0.0;
    RSS( ref_patch_locate( ref_patch, xyz, uv ), "loc");
    RWDS( 0.0, uv[0], -1.0, "uv[0]");
    RWDS( 1.0, uv[1], -1.0, "uv[1]");

    xyz[0] = 0.3; xyz[1] = 0.4; xyz[2] = 0.0;
    RSS( ref_patch_locate( ref_patch, xyz, uv ), "loc");
    RWDS( 0.3, uv[0], -1.0, "uv[0]");
    RWDS( 0.4, uv[1], -1.0, "uv[1]");

    RSS( ref_patch_free( ref_patch ), "free" );
  }

  return 0;
}
