
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

int main( int argc, char *argv[] )
{
  Grid *grid;
  char *filename;

  if (argc == 2) {
    filename = argv[1];
  }else{
    printf("\nUsage: refine project.fgrid\n\n");
    filename = "../test/om6_inv08.fgrid";
    printf("running default filename %s\n",filename);
  }
  
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

  printf("calling CADGeom_Start ... \n");
  if ( ! CADGeom_Start( ) ){
    printf("Yo! it broke.\n");
  }  

  printf("calling CADGeom_Load ... \n");
  if ( ! CADGeom_LoadPart( "../test/om6ff25" ) ){
    printf("Yo! it broke.\n");
  }  

  printf("calling CADGeom_Stop ... \n");
  if ( ! CADGeom_Stop( ) ){
    printf("Yo! it broke.\n");
  }  

  printf("Done. \n");

  return;
}
