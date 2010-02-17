
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:mike.park@nasa.gov
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "unravel_fortran.h"

#include "grid.h"

static Grid *grid = NULL;

void unravel_start_( void )
{
  if (NULL != grid) gridFree(grid);

  grid = gridCreate(200, 200, 0, 0);

}

void unravel_tet_( int *c2n, double *x, double *y, double *z )
{
  int cell_node;
  int local_nodes[4];

  for ( cell_node = 0 ; cell_node < 4 ; cell_node++ )
    {
      local_nodes[cell_node] = gridGlobal2Local(grid, c2n[cell_node]-1 );
      if ( EMPTY == local_nodes[cell_node] )
	{
	  local_nodes[cell_node] = 
	    gridAddNodeWithGlobal(grid, 
				  x[cell_node], y[cell_node], z[cell_node], 
				  c2n[cell_node]-1 );
	  gridFreezeNode( grid, local_nodes[cell_node] );
	}
    }
  gridAddCell(grid, 
	      local_nodes[0], local_nodes[1], 
	      local_nodes[2], local_nodes[3] );
}

void unravel_thaw_( int *nodeid, int *status )
{
  *status = 0;
  if ( grid != gridThawNode( grid, (*nodeid)-1 ) ) *status = 1;
}

void unravel_it_( int *status )
{
  *status = -1;
}

void unravel_xyz_( int *nodeid, double *x, double *y, double *z, int *status )
{
  int local_node;
  double *xyz;

  *status = 0;
  
  local_node = gridGlobal2Local(grid, (*nodeid)-1 );
  if ( grid != gridNodeXYZ( grid, local_node, xyz ) )  
    {
      *status = -1;
      return;
    }
  *x = xyz[0];
  *y = xyz[1];
  *z = xyz[2];
}

void unravel_cleanup_( void )
{
 gridFree(grid);
 grid = NULL;
}

