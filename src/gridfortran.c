
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "gridfortran.h"
#include "grid.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "gridcad.h"
#include "gridinsert.h"

static Grid *grid;

int gridcreate_( int *nnode, double *x, double *y, double *z ,
		 int *ncell, int *maxcell, int *c2n )
{
  int node, cell;
  int nodes[4];
  double xyz[3];
  grid = gridCreate( *nnode, *ncell, 5000, 0);
  for ( node=0; node<*nnode; node++) gridAddNode(grid,x[node],y[node],z[node]);
  printf("populated grid object with %d nodes\n",gridNNode(grid));
  for ( cell=0; cell<*ncell; cell++) gridAddCell( grid,
						  c2n[cell+0*(*maxcell)] - 1,
						  c2n[cell+1*(*maxcell)] - 1,
						  c2n[cell+2*(*maxcell)] - 1,
						  c2n[cell+3*(*maxcell)] - 1 );
  printf("populated grid object with %d cells\n",gridNCell(grid));
  printf(" min AR %17.15f\n",gridMinAR(grid));
}

int gridfree_( )
{
  gridFree(grid);
}

int gridinsertboundary_( int *faceId, int *nnode, int *inode, 
			 int *nface, int *ndim, int *f2n )
{
  int face;
  int node0, node1, node2;
  for(face=0;face<*nface;face++){
    node0 = f2n[face+0*(*nface)] - 1;
    node1 = f2n[face+1*(*nface)] - 1;
    node2 = f2n[face+2*(*nface)] - 1;
    node0 = inode[node0] - 1;
    node1 = inode[node1] - 1;
    node2 = inode[node2] - 1;
    gridAddFace(grid, node0, node1, node2, *faceId);
  }
  printf( " boundary %4d has %6d nodes and %6d faces\n",
	  *faceId,*nnode,*nface);
  printf( " %d total faces with min MR of %18.15f\n",
	  gridNFace(grid),gridMinFaceMR(grid));
}

int gridsetmap_( int *nnode, double* map )
{
  int node;
  for ( node=0; node<*nnode; node++) 
    gridSetMap( grid, node,
		map[0+6*node], map[1+6*node], map[2+6*node],
		map[3+6*node], map[4+6*node], map[5+6*node] );
}

int gridswap_( )
{
  gridSwap(grid);
  printf(" post swap minn AR %17.15f\n",gridMinAR(grid));
}

int gridsmoothvolume_( )
{
  gridSmoothVolume(grid);
  printf(" post smooth min AR %17.15f\n",gridMinAR(grid));
}

int gridadaptwithoutcad_( double *minLength, double *maxLength )
{
  gridAdaptWithOutCAD(grid,*minLength, *maxLength);
  printf(" post adapt min AR %17.15f\n",gridMinAR(grid));
}

int gridwritetecplotsurfacezone_( )
{
  gridWriteTecplotSurfaceZone(grid);
}
