
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
#include "gridmove.h"
#include "gridfiller.h"
#include "CADGeom/CADGeom.h"

#define PRINT_STATUS printf("minimum Thawed Aspect Ratio %8.6f Mean Ratio %8.6f Volume %10.6e\n", gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid)); fflush(stdout);

#define DUMP_TEC if (tecplotOutput) {iview++;printf("Frame %d\n",iview);gridWriteTecplotSurfaceZone(grid,NULL);}

#define STATUS DUMP_TEC PRINT_STATUS

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  Grid *grid;
  int jmax;
  char project[256];
  char adaptfile[256];
  char linesfile[256];
  char outputProject[256];
  char outputFAST[256];
  char outputlines[256];
  int i, j, oldSize, newSize;
  double ratio=0.6;
  double spacing =1.0/3.0;
  double minAR=-1.0;
  double ratioSplit, ratioCollapse;
  GridBool projected;
  GridBool GridMoveProjection = FALSE;
  GridBool tecplotOutput = FALSE;
  int iview = 0;
  int maxnode = 50000;
  char modeler[81];

  sprintf( modeler,       "" );
  sprintf( project,       "" );
  sprintf( outputProject, "" );
  sprintf( adaptfile,     "" );    
  sprintf( linesfile,     "" );    
  sprintf( outputFAST,    "" );

  i = 1;
  while( i < argc ) {
    if( strcmp(argv[i],"-p") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-p argument %d: %s\n",i, project);
      sprintf( modeler, "Unknown" );
    } else if( strcmp(argv[i],"-felisa") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-felisa argument %d: %s\n",i, project);
      sprintf( modeler, "FELISA" );
    } else if( strcmp(argv[i],"-open") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-open argument %d: %s\n",i, project);
      sprintf( modeler, "OpenCASCADE" );
    } else if( strcmp(argv[i],"-proe") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-proe argument %d: %s\n",i, project);
      sprintf( modeler, "Pro/ENGINEER" );
    } else if( strcmp(argv[i],"-parasolid") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-parasolid argument %d: %s\n",i, project);
      sprintf( modeler, "Parasolid" );
    } else if( strcmp(argv[i],"-catia") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-catia argument %d: %s\n",i, project);
      sprintf( modeler, "CatiaV5" );
    } else if( strcmp(argv[i],"-ug") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-ug argument %d: %s\n",i, project);
      sprintf( modeler, "UniGraphics" );
    } else if( strcmp(argv[i],"-sw") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-sw argument %d: %s\n",i, project);
      sprintf( modeler, "SolidWorks" );
    } else if( strcmp(argv[i],"-o") == 0 ) {
      i++; sprintf( outputProject, "%s", argv[i] );
      printf("-o argument %d: %s\n",i, outputProject);
    } else if( strcmp(argv[i],"-r") == 0 ) {
      i++; ratio = atof(argv[i]);
      printf("-r argument %d: %f\n",i, ratio);
    } else if( strcmp(argv[i],"-s") == 0 ) {
      i++; spacing = atof(argv[i]);
      printf("-s argument %d: %f\n",i, spacing);
    } else if( strcmp(argv[i],"-v") == 0 ) {
      i++; minAR = atof(argv[i]);
      printf("-v argument %d: %f\n",i, minAR);
    } else if( strcmp(argv[i],"-f") == 0 ) {
      i++; sprintf( linesfile, "%s", argv[i] );
      printf("-l argument %d: %s\n",i, linesfile);
    } else if( strcmp(argv[i],"-n") == 0 ) {
      i++; maxnode = atoi(argv[i]);
      printf("-n argument %d: %d\n",i, maxnode);
     } else if( strcmp(argv[i],"-t") == 0 ) {
      tecplotOutput = TRUE;
      printf("-t argument %d\n",i);
   } else if( strcmp(argv[i],"-h") == 0 ) {
      printf("Usage: flag value pairs:\n");
#ifdef HAVE_CAPRI2
      printf(" -felisa input FELISA project name\n");
      printf(" -open input OpenCASCADE project name\n");
      printf(" -proe input Pro/E project name\n");
      printf(" -parasolid input Parasolid project name\n");
      printf(" -catia input Catia V5 project name\n");
      printf(" -ug input Unigraphics project name\n");
      printf(" -sw input SolidWorks project name\n");
#else
      printf(" -p input project name\n");
#endif
      printf(" -o output project name\n");
      printf(" -r initial edge length ratio for adapt\n");
      printf(" -s uniform grid size\n");
      printf(" -v freeze cells with small aspect ratio (viscous)\n");
      printf(" -f freeze nodes in this .lines file\n");
      printf(" -n max number of nodes in grid\n");
      printf(" -t write tecplot zones durring adaptation\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(strcmp(modeler,"")==0)       sprintf(modeler,"FELISA" );
  if(strcmp(project,"")==0)       sprintf(project,"../test/box1" );
  if(strcmp(outputProject,"")==0) sprintf(outputProject,"%s_out", project );
  if(strcmp(outputFAST,"")==0)    sprintf(outputFAST,"%s.fgrid",outputProject);

  printf("running project %s\n",project);
  grid = gridLoadPart( modeler, project, maxnode );

  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: loaded part does not have right handed boundaries\n");

  printf("restart grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  if(strcmp(linesfile,"")!=0) {
    printf("loading lines from file %s\n",linesfile);
    linesLoad(gridLines(grid), linesfile);
    printf("freezing line nodes...\n");
    gridFreezeLinesNodes(grid);
    printf("nodes frozen %d\n",gridNFrozen(grid));
  }

  if (minAR > 0) {
    printf("freezing cells with AR smaller than %f\n",minAR);
    gridFreezeSmallARCells(grid, minAR);
  }

  printf("Spacing reset.\n");
  gridResetSpacing(grid);
  printf("spacing set to constant.\n");
  gridSetGlobalMap(grid, 1.0/spacing/spacing, 0.0, 0.0, 1.0/spacing/spacing, 0.0, 1.0/spacing/spacing);
  STATUS;

  for (i=0;i<3;i++){
    projected = ( grid == gridRobustProject(grid));
    if (projected) {
      printf("edge swapping grid...\n");gridSwap(grid);
      printf("node smoothing grid...\n");gridSmooth(grid,1.0,1.0);
    }else{
      printf("could not project grid. stop.\n");
      return 1;
    }
  }
  STATUS;

  oldSize = 1;
  newSize = gridNNode(grid);
  jmax = 40;
  for ( j=0; 
	(j<jmax) && 
	  (((double)ABS(newSize-oldSize)/(double)oldSize)>0.001) || 
	  !projected;
	j++){

    ratioCollapse = 1.0*ratio;
    ratioSplit   = 1.0/sqrt(ratio);
    printf("adapt, ratio %4.2f, collapse limit %8.5f, refine limit %10.5f\n",
	   ratio, ratioCollapse, ratioSplit );
    gridAdapt(grid,ratioCollapse,ratioSplit, !GridMoveProjection);

    oldSize = newSize;
    newSize = gridNNode(grid) ;
    printf("%02d new size: %d nodes %d faces %d cells %d edge elements.\n",
	   j, gridNNode(grid),gridNFace(grid),gridNCell(grid),gridNEdge(grid));
    STATUS;

    for (i=0;i<2;i++){
      projected = ( grid == gridRobustProject(grid));
      if (projected) {
	printf("edge swapping grid...\n");gridSwap(grid);
	printf("node smoothing grid...\n");gridSmooth(grid,1.0,1.0);
	printf("node smoothing surface grid...\n");gridSmoothFaceMR(grid,1.0);
      }else{
	printf("node smoothing volume grid...\n");gridSmoothVolume(grid);
      }
    }
    STATUS;
  }

  for (i=0;i<5;i++){
    projected = ( grid == gridRobustProject(grid));
    if (projected) {
      printf("edge swapping grid...\n");gridSwap(grid);
      printf("node smoothing grid...\n");gridSmooth(grid,1.0,1.0);
      printf("node smoothing surface grid...\n");gridSmoothFaceMR(grid,1.0);
    }else{
      printf("node smoothing volume grid...\n");gridSmoothVolume(grid);
    }
    STATUS;
  }

  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: modifed grid does not have right handed boundaries\n");

  printf("writing output project %s\n",outputProject);
  gridSavePart( grid, outputProject );

  printf("writing output FAST file %s\n",outputFAST);
  gridExportFAST( grid, outputFAST );

  if(strcmp(linesfile,"")!=0) {
    sprintf(outputlines,"%s.lines",outputProject);
    printf("saving lines restart %s\n",outputlines);
    linesSave(gridLines(grid),outputlines);
  }

  printf("Done.\n");

  return 0;
}

