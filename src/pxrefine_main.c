
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include "grid.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "gridinsert.h"
#include "gridcad.h"

#define PRINT_STATUS {double l0,l1;gridEdgeRatioRange(grid,&l0,&l1);printf("Len %12.5e %12.5e AR %8.6f MR %8.6f Vol %10.6e\n", l0,l1, gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid)); fflush(stdout);}

#define STATUS { \
  PRINT_STATUS; \
}


#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  Grid *grid;
  char modeler[256];
  char gri_input[256];
  char metric_input[256];
  char gri_output[256];

  int i;
  int iview=0;

  sprintf( modeler,    "Unknown" );
  sprintf( gri_input,    "" );
  sprintf( gri_output,   "" );
  sprintf( metric_input, "" );

  i = 1;
  while( i < argc ) {
    if( strcmp(argv[i],"-p") == 0 ) {
      i++; sprintf( gri_input, "%s", argv[i] );
      printf("-p argument %d: %s\n",i, gri_input);
    } else if( strcmp(argv[i],"-m") == 0 ) {
      i++; sprintf( metric_input, "%s", argv[i] );
      printf("-m argument %d: %s\n",i, metric_input);
    } else if( strcmp(argv[i],"-o") == 0 ) {
      i++; sprintf( gri_output, "%s", argv[i] );
      printf("-o argument %d: %s\n",i, gri_output);
    } else if( strcmp(argv[i],"-h") == 0 ) {
      printf("Usage: flag value pairs:\n");
      printf(" -p input .gri name\n");
      printf(" -m input .metric name\n");
      printf(" -o output .gri name\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if ( (strcmp(gri_input,"")==0) ||
       (strcmp(gri_output,"")==0) ) {
    printf("no input of output specified.\n");
    printf("Done.\n");  
    return 1;
  }

  grid = gridImportGRI( gri_input );

  printf("grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  gridWriteTecplotSurfaceGeom(grid,NULL);

  printf("Spacing reset.\n");
  gridResetSpacing(grid);

  gridSetCostConstraint(grid,
			gridCOST_CNST_VOLUME );


  STATUS;

  gridExportFAST( grid, "grid_orig.fgrid" );

  if (!gridRightHandedBoundary(grid)) {
    printf("ERROR: loaded part does not have right handed boundaries\n");
    return 1;
  }

  printf("Done.\n");
  
  return 0;
}

