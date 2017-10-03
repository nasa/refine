
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
      printf("usage : %s input.{ngp|fgrid|grd|gri} output.{ngp|fgrid|dat} \n", argv[0] );  
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
  } else if( strcmp(&file_name[end_of_string-4],".gri") == 0 ) {
    printf("gri input file %s\n", file_name);
    grid = gridImportGRI( file_name );
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

