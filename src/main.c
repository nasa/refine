
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
#include "gridcad.h"
#include "gridfiller.h"
#include "layer.h"
#include <CADGeom/CADGeom.h>

#define PRINT_STATUS printf("minimum Thawed Aspect Ratio %8.6f Mean Ratio %8.6f Volume %10.6e\n", gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid));

#define DUMP_TEC iview++;printf("Frame %d\n",iview);gridWriteTecplotSurfaceZone(grid);

#define STATUS DUMP_TEC PRINT_STATUS

int main( int argc, char *argv[] )
{
  Grid *grid;
  Layer *layer;
  int jmax;
  double height;
  char project[256];
  char adaptfile[256], outputProject[256], outputFAST[256];
  int i, j, oldSize, newSize;
  int wiggleSteps, wiggle;
  double ratio=0.3;
  double ratioRefine, ratioCollapse;
  bool projected;
  bool boundaryLayerGrid = FALSE;
  bool debugInsert = FALSE;
  int iview = 0;

  sprintf( project,       "" );
  sprintf( outputProject, "" );
  sprintf( adaptfile,     "" );    
  sprintf( outputFAST,    "" );

  i = 1;
  while( i < argc ) {
    if( strcmp(argv[i],"-p") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-p argument %d: %s\n",i, project);
    } else if( strcmp(argv[i],"-o") == 0 ) {
      i++; sprintf( outputProject, "%s", argv[i] );
      printf("-o argument %d: %s\n",i, outputProject);
    } else if( strcmp(argv[i],"-a") == 0 ) {
      i++; sprintf( adaptfile,"%s",argv[i]  );
      printf("-a argument %d: %s\n",i, adaptfile);
    } else if( strcmp(argv[i],"-l") == 0 ) {
      boundaryLayerGrid = TRUE;
      printf("-l argument %d, ignoring -a \n",i);
    } else if( strcmp(argv[i],"-r") == 0 ) {
      i++; ratio = atof(argv[i]);
      printf("-r argument %d: %f\n",i, ratio);
    } else if( strcmp(argv[i],"-i") == 0 ) {
      debugInsert = TRUE;
      printf("-i argument %d\n",i);
    } else if( strcmp(argv[i],"-h") == 0 ) {
      printf("Usage: flag value pairs:\n");
      printf(" -p input project name\n");
      printf(" -o output project name\n");
      printf(" -a party project_adapt_hess file name\n");
      printf(" -l make a boundary layer grid -a ignored\n");
      printf(" -r initial edge length ratio for adapt\n");
      printf(" -i insert final advancing layer (debug)\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(strcmp(project,"")==0)       sprintf(project,"../test/om6" );
  if(strcmp(outputProject,"")==0) sprintf(outputProject,"%s_out", project );
  if(strcmp(adaptfile,"")==0)     sprintf(adaptfile,"none");
  if(strcmp(adaptfile,"")==0)     sprintf(adaptfile,"%s_adapt_hess",project);
  if(strcmp(outputFAST,"")==0)    sprintf(outputFAST,"%s.fgrid",outputProject);

  if(boundaryLayerGrid || debugInsert ) sprintf(adaptfile,"none");

  printf("running project %s\n",project);
  grid = gridLoadPart( project, 500000 );

  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: loaded part does not have right handed boundaries\n");

  printf("restart grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  if(strcmp(adaptfile,"none")==0) {
    printf("adapt parameter >none< selected. Spacing reset.\n");
    gridResetSpacing(grid);
    if (boundaryLayerGrid) {
      layer = formAdvancingFront(grid,project);
    }else{
      if (debugInsert) {
	layer = formAdvancingFront(grid,project);
	printf("Inserting Phantom front.\n");
	layerInsertPhantomFront( layer, 0.22 );
	ratio=0.8;
      }else{
	printf("Scaling spacing to refine a sphere.\n");
	gridScaleSpacingSphere(grid, 0.0, 0.0, 0.0, 1.0, 0.7 );
      }
    }
  }else{
    printf("reading adapt parameter from file %s ...\n",adaptfile);
    gridImportAdapt(grid, adaptfile); // Do not sort nodes before this call.
  }
  STATUS;

  for (i=0;i<3;i++){
    projected = ( grid == gridRobustProject(grid));
    if (projected) {
      printf("edge swapping grid...\n");gridSwap(grid);
      printf("node smoothing grid...\n");gridSmooth(grid);
    }else{
      printf("node smoothing volume grid...\n");gridSmoothVolume(grid);
    }
  }
  STATUS;

  oldSize = 1;
  newSize = gridNNode(grid);
  jmax = 40;
  for ( j=0; (j<jmax) && (
	(ratio < 0.99) || 
	  (((double)ABS(newSize-oldSize)/(double)oldSize)>0.001) ||
	  !projected );
	j++){

    if (boundaryLayerGrid) {
      height = 0.0016*pow(1.2,j);
      layerTerminateNormalWithSpacing(layer,height*2.5);
      if (layerNActiveNormal(layer) == 0 ) jmax=0;
      printf("insert layer height = %f\n",height);
      wiggleSteps = MIN(4,(int)(height/0.001)+1);
      height = height / (double)wiggleSteps;
      layerVisibleNormals(layer);
      layerAdvance(layer,height);
      for (i=1;i<wiggleSteps;i++) {
	printf("edge swapping grid...\n");gridSwap(grid);
	layerSmoothLayerNeighbors(layer );
	printf("node smoothing grid...\n");gridSmooth(grid);
	gridAdapt(grid,0.6,1.4);
	printf("edge swapping grid...\n");gridSwap(grid);
	layerSmoothLayerNeighbors(layer );
	printf("node smoothing grid...\n");gridSmooth(grid);
	printf("wiggle step %d of %d, minAR %8.5f\n",i+1,wiggleSteps,gridMinThawedAR(grid));
	layerWiggle(layer,height);
	//printf("minimum Volume %12.8e\n", gridMinVolume(grid));
      }
      printf("edge swapping grid...\n");gridSwap(grid);
	layerSmoothLayerNeighbors(layer );
      printf("node smoothing grid...\n");gridSmooth(grid);
      STATUS;
    }
    if (ratio<0.01) ratio = 0.01;
    if (ratio>1.0) ratio = 1.0;
    ratioCollapse = 0.4*ratio;
    ratioRefine   = 1.5/ratio;
    if (boundaryLayerGrid) {
      printf("adapt, ratio %4.2f, collapse limit %8.5f, refine limit %10.5f\n",
	     ratio, 0.4, 1.5 );
      gridAdapt(grid,0.4,1.5);
    }else{
      printf("adapt, ratio %4.2f, collapse limit %8.5f, refine limit %10.5f\n",
             ratio, ratioCollapse, ratioRefine );
      gridAdapt(grid,ratioCollapse,ratioRefine);
    }
    oldSize = newSize;
    newSize = gridNNode(grid) ;
    printf("%02d new size: %d nodes %d faces %d cells %d edge elements.\n",
	   j, gridNNode(grid),gridNFace(grid),gridNCell(grid),gridNEdge(grid));
    STATUS;
    if (debugInsert) 
      { layerVerifyPhantomEdges( layer ); layerVerifyPhantomFaces( layer ); }
        
    for (i=0;i<2;i++){
      projected = ( grid == gridRobustProject(grid));
      if (projected) {
	printf("edge swapping grid...\n");gridSwap(grid);
	printf("node smoothing grid...\n");gridSmooth(grid);
	if (((double)ABS(newSize-oldSize)/(double)oldSize)<0.3)
	  ratio = ratio + 0.025;
      }else{
	printf("node smoothing volume grid...\n");gridSmoothVolume(grid);
	ratio = ratio - 0.05;
      }
    }
    STATUS;
    if (boundaryLayerGrid) {
    }else{
      if (!debugInsert) gridFreezeGoodNodes(grid,0.6,0.4,1.5);
      printf("nodes frozen %d\n",gridNFrozen(grid));
    }
  }

  if (debugInsert) 
    { layerVerifyPhantomEdges( layer ); layerVerifyPhantomFaces( layer ); }


  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: modifed grid does not have right handed boundaries\n");

  printf("writing output project %s\n",outputProject);
  gridSavePart( grid, outputProject );

  printf("writing output FAST file %s\n",outputFAST);
  gridExportFAST( grid, outputFAST );

  printf("Done.\n");

  return;
}

