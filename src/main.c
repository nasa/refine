
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include <CADGeom/CADGeom.h>

Grid *gridLoadPart( char *project);

int main( int argc, char *argv[] )
{
  Grid *grid;
  char *filename;
  char *project;


  filename = "../test/om6_inv08.fgrid";
  printf("running default filename %s\n",filename);
  project = "../test/om6";
  printf("running default project %s\n",project);
  
  grid = gridImportFAST( filename );
  
  printf("orig size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));
  printf("minimum Aspect Ratio %12f\n",gridMinAR(grid));
  printf("minimum Volume %12.8e\n",gridMinVolume(grid));

  printf("edge swapping grid...\n");
  gridSwap(grid);

  printf("new size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));
  printf("minimum Aspect Ratio %12f\n",gridMinAR(grid));
  printf("minimum Volume %12.8e\n",gridMinVolume(grid));

  gridFree(grid);

  grid = gridLoadPart( project );

  printf("Done.\n");

  return;
}

Grid *gridLoadPart( char *project)
{
  Grid *grid;
  int vol=1;
  UGridPtr ugrid;
  int gridDimensions[3];
  int nnode, nface, ncell;
  int maxnode, maxface, maxcell;
  
  printf("calling CADGeom_Start ... \n");
  if ( ! CADGeom_Start( ) ){
    printf("ERROR: CADGeom_Start broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  

  printf("calling CADGeom_Load for project <%s> ... \n",project);
  if ( ! CADGeom_LoadPart( project ) ){
    printf("ERROR: CADGeom_LoadPart broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  if (NULL == (ugrid = CADGeom_VolumeGrid(vol)) ) {
    printf("ERROR: Can not find grid in restart. \n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  if( !UGrid_GetDims(ugrid,gridDimensions) ) {
    printf("ERROR: Could not get grid dimensions. \n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  nnode = gridDimensions[0];
  nface = gridDimensions[1];
  ncell = gridDimensions[2];

  printf("ugrid size: %d nodes %d faces %d cells.\n",nnode,nface,ncell);

  maxnode = nnode * 100;
  maxface = nface * 100;
  maxcell = ncell * 100;

  printf("max grid size: %d nodes %d faces %d cells.\n",
	 maxnode,maxface,maxcell);

  printf("calling CADGeom_Stop ... \n");
  if ( ! CADGeom_Stop( ) ){
    printf("ERROR: CADGeom_Stop broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  

  return grid;
}

 
