
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "grid.h"
#include "gridmetric.h"

int main( int argc, char *argv[] )
{
  Grid *grid;
  char *file_name;
  int end_of_string;
  double min_vol;
  GridBool valid_boundary;

  printf("refine %s\n",VERSION);

  if ( 3 != argc )
    {
      printf("usage : %s input.{ngp|fgrid|grd} output.{ngp|fgrid|dat} \n", argv[0] );  
      return 1;
    }

  file_name = argv[1];
  end_of_string = strlen(file_name);
 
  grid = NULL;
  if( strcmp(&file_name[end_of_string-4],".ngp") == 0 ) {
    printf("ngp input file %s\n", file_name);
    grid = gridImportNGP( file_name );
  } else if( strcmp(&file_name[end_of_string-6],".fgrid") == 0 ) {
    printf("fast input file %s\n", file_name);
    grid = gridImportFAST( file_name );
  } else if( strcmp(&file_name[end_of_string-4],".grd") == 0 ) {
    printf("fieldview input file %s\n", file_name);
    grid = gridImportFV( file_name );
  } else {
    printf("input file name extension unknown %s\n", file_name);
  }

  if ( NULL == grid )
    {
      printf("grid import failed\n");
      return 1;
    }

  min_vol = gridMinVolume(grid);
  valid_boundary = gridRightHandedBoundary(grid);
  printf("min volume %e with %s boundary faces\n",
	 min_vol,(valid_boundary?"valid":"invalid"));

  file_name = argv[2];
  end_of_string = strlen(file_name);
 
  if( strcmp(&file_name[end_of_string-4],".ngp") == 0 ) {
    printf("ngp output file %s\n", file_name);
    gridExportNGP( grid, file_name );
  } else if( strcmp(&file_name[end_of_string-6],".fgrid") == 0 ) {
    printf("fast output file %s\n", file_name);
    gridExportFAST( grid, file_name );
  } else if( strcmp(&file_name[end_of_string-4],".dat") == 0 ) {
    printf("tecplot output file %s\n", file_name);
    gridWriteTecplotSurfaceGeom( grid, file_name );
  } else {
    printf("output file name extension unknown %s\n", file_name);
    return 1;
  }

  printf("Done.\n");

  return 0;
}

