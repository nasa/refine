
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
#include "gridmpi.h"
#include "gridfortran.h"
#include "gridfiller.h"
#include "CADGeom/CADGeom.h"

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  Grid *grid;
  Queue *queue;

  char modeler[81];
  char project[256];
  char outputProject[256];
  char outputFAST[256];
  char filename[256];

  int i;

  int maxnode;

  sprintf( modeler,       "" );
  sprintf( project,       "" );
  sprintf( outputProject, "" );
  sprintf( outputFAST,    "" );

  maxnode = 50000;

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
    } else if( strcmp(argv[i],"-m") == 0 ) {
      i++; maxnode = atoi(argv[i]);
      printf("-m argument %d: %d\n",i, maxnode);
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
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(strcmp(modeler,"")==0)       sprintf(modeler,"FELISA" );
  if(strcmp(project,"")==0)       sprintf(project,"../test/om6" );
  if(strcmp(outputProject,"")==0) sprintf(outputProject,"%s_out", project );
  if(strcmp(outputFAST,"")==0)    sprintf(outputFAST,"%s.fgrid",outputProject);

  printf("running project %s\n",project);
  grid = gridLoadPart( modeler, project, maxnode );

  gridSetPartId(grid, 0 );
  queue = queueCreate( 9 );

  printf("initial grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  printf("Spacing reset.\n");
  gridResetSpacing(grid);
  printf("spacing set to high ar in x.\n");

  gridScaleSpacingSphereDirection( grid, 0.0, 0.0, 0.0, DBL_MAX,
				   10.0, 1.0, 1.0 );

  if (grid!=gridRobustProject(grid)) {
    printf("could not project grid. stop.\n");
    return 1;
  }

  printf("Swapping.\n");
  gridSwap(grid,-1);
  printf("Swapping.\n");
  gridSwap(grid,-1);
  printf("Swapping.\n");
  gridSwap(grid,-1);

  printf("writing output project %s\n",outputProject);
  gridSavePart( grid, outputProject );

  printf("Done.\n");

  return 0;
}

