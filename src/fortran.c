
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "fortran.h"
#include "grid.h"

static Grid *grid;

int gridcreate( int *maxnode, int *maxcell, int *maxface, int *maxedge )
{
  printf("called gridcreate\n")
  grid = gridCreate( *maxnode, *maxcell, *maxface, *maxedge);
}

int gridcreate_( int *maxnode, int *maxcell, int *maxface, int *maxedge )
{
  printf("called gridcreate_\n")
  grid = gridCreate( *maxnode, *maxcell, *maxface, *maxedge);
}

int gridcreate__( int *maxnode, int *maxcell, int *maxface, int *maxedge )
{
  printf("called gridcreate__\n")
  grid = gridCreate( *maxnode, *maxcell, *maxface, *maxedge);
}

int GRIDCREATE( int *maxnode, int *maxcell, int *maxface, int *maxedge )
{
  printf("called GRIDCREATE\n")
  grid = gridCreate( *maxnode, *maxcell, *maxface, *maxedge);
}

