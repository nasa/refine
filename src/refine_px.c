
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


static void usage( char *executable );
static void usage( char *executable )
{
  printf("Usage: %s -g input.gri -g projectx.metric -o output.gri\n",
	 executable );
  printf(" -g input .gri name\n");
  printf(" -m input .metric name\n");
  printf(" -o output .gri name\n");
}

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  Grid *grid;
  char gri_input[256] = "";
  char metric_input[256] = "";
  char gri_output[256] = "";

  int i;

  int active_edges, out_of_tolerence_edges;

  printf("refine %s\n",VERSION);

  i = 1;
  while( i < argc ) {
    if( strcmp(argv[i],"-g") == 0 ) {
      i++; sprintf( gri_input, "%s", argv[i] );
      printf("-g argument %d: %s\n",i, gri_input);
    } else if( strcmp(argv[i],"-m") == 0 ) {
      i++; sprintf( metric_input, "%s", argv[i] );
      printf("-m argument %d: %s\n",i, metric_input);
    } else if( strcmp(argv[i],"-o") == 0 ) {
      i++; sprintf( gri_output, "%s", argv[i] );
      printf("-o argument %d: %s\n",i, gri_output);
    } else if( strcmp(argv[i],"-h") == 0 ) {
      usage( argv[0] );
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if ( (strcmp(gri_input,"")==0) ||
       (strcmp(gri_output,"")==0) ) {
    printf("no input or output file names specified.\n");
    usage( argv[0] );
    printf("Done.\n");  
    return 1;
  }

  grid = gridImportGRI( gri_input );

  if ( NULL == grid ) 
    {
      printf("read of %s failed. stop\n",gri_input);
      return(1);
    }

  printf("grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  gridRobustProject(grid);

  gridSetCostConstraint(grid,
			gridCOST_CNST_VOLUME );

  STATUS;

  if ( !(strcmp(metric_input,"")==0) )
    {
      printf("reading metric file %s\n",metric_input);
      gridImportAdapt(grid,metric_input);
      STATUS;
    } else {
      printf("Spacing reset.\n");
      gridResetSpacing(grid);
    }

  STATUS;

  gridWriteTecplotSurfaceGeom(grid,NULL);
  {
    int iteration;
    int iterations = 10;
    double ratioSplit, ratioCollapse;

    ratioCollapse = 0.3;
    ratioSplit    = 1.0;      

    STATUS;

    printf("edge swapping grid...\n");gridSwap(grid,0.9);
    STATUS;

    gridAdapt(grid, ratioCollapse, ratioSplit);
    STATUS;
    
    for ( iteration=0; (iteration<iterations) ; iteration++){
      
      for (i=0;i<1;i++){
	printf("edge swapping grid...\n");gridSwap(grid,0.9);
	STATUS;
	printf("node smoothin grid...\n");gridSmooth(grid,0.9,0.5);
	STATUS;
      }

      gridEdgeRatioTolerence(grid, ratioSplit, ratioCollapse,
			     &active_edges, &out_of_tolerence_edges );

      printf("edges %d of %d (%6.2f%%) out of tol\n",
	     out_of_tolerence_edges, active_edges, 
	     100.0*(double)out_of_tolerence_edges/(double)active_edges);

      if ( 0.001 > (double)out_of_tolerence_edges/(double)active_edges)
	break;

      gridAdapt(grid, ratioCollapse, ratioSplit);
      STATUS;
    
      gridWriteTecplotSurfaceGeom(grid,NULL);
    }
  

  }


  gridExportFAST( grid, "grid_orig.fgrid" );
  if ( !(strcmp(gri_output,"")==0) ) 
    gridExportGRI( grid, gri_output );


  printf("Done.\n");
  
  return 0;
}

