
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
#include "gridshape.h"
#include "gridmove.h"
#include "gridmpi.h"
#include "gridfiller.h"
#include "CADGeom/CADGeom.h"

#define PRINT_STATUS printf("minimum Thawed Aspect Ratio %8.6f Mean Ratio %8.6f Volume %10.6e\n", gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid)); fflush(stdout);

#define DUMP_TEC if (tecplotOutput) { \
 iview++;printf("Frame %d\n",iview);\
 gridWriteTecplotSurfaceGeom(grid,NULL); \
  {int cell, nodes[4]; \
    for (cell=0;cell<gridMaxCell(grid);cell++) \
      if (grid==gridCell(grid, cell, nodes)) { \
        if ( -0.5 > gridAR(grid,nodes) ) \
          gridWriteTecplotCellGeom(grid,nodes,NULL,NULL); \
   }}; }

#define STATUS DUMP_TEC PRINT_STATUS

/* new format
 * order
 * number of edges
 * edge_n0 edge_n1 x y z x y z ...
 * number of triangle points
 * n0 n1 n2 b0 b1 b2 x y z
 */

Grid *gridDumpBentEdgesForPX(Grid *grid, int order, char *filename)
{
  int conn, total;
  int midnode;
  int nodes[2];
  double xyz[3];
  FILE *file;

  gridCreateConn(grid);
  total = 0;
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    if (0 != gridParentGeometry(grid, nodes[0], nodes[1]) ) {
      total++;
    }
  }
  total = total * (order-1);

  printf("%d edges bent of %d total edges.",total,gridNConn(grid));
  file = fopen(filename,"w");
  fprintf(file,"%10d edges in 1-base numbering,\n",total);
  
  total = 0;
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    if (0 != gridParentGeometry(grid, nodes[0], nodes[1]) ) {
      for (midnode=0;midnode<(order-1);midnode++) {
	total++;
	gridCurvedEdgeMidpoint(grid,nodes[0], nodes[1], xyz);
	fprintf(file,"%10d%10d%24.15e%24.15e%24.15e\n",
		nodes[0]+1,nodes[1]+1,xyz[0],xyz[1],xyz[2]);
      }
    }
  }
  fclose(file);

  gridEraseConn(grid);
  return grid;
}

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  Grid *grid;
  char project[256];
  char outputProject[256];
  char filename[256];
  double ratioSplit, ratioCollapse;
  int maxnode = 50000;
  int i;
  int cycle, invalid;
  int iview;
  GridBool tecplotOutput = FALSE;

  char modeler[81];

  sprintf( modeler,       "" );
  sprintf( project,       "" );
  sprintf( outputProject, "" );

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
      printf(" -n max number of nodes in grid\n");
      printf(" -t write tecplot zones durring adaptation\n");
      printf(" -h help\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(strcmp(modeler,"")==0)       sprintf(modeler,"FELISA" );

  if(strcmp(project,"")==0)       sprintf(project,"default_project" );

  if(strcmp(outputProject,"")==0) sprintf(outputProject,"%s_px", project );

  printf("running project %s\n",project);
  grid = gridLoadPart( modeler, project, maxnode );

  printf("restart grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  if (!gridRightHandedBoundary(grid)) {
    printf("ERROR: loaded part does not have right handed boundaries\n");
    return 1;
  }
  if (!gridRightHandedBoundaryUV(grid)) {
    int faceId;
    printf("ERROR: loaded part does not have right handed UV parameters\n");
    for(faceId=1;faceId<=gridNGeomFace(grid);faceId++)
      gridWriteTecplotGeomFaceUV(grid,"faceParameters.t",faceId);
    gridCloseTecplotGeomFile(grid);
    gridRobustProject(grid); 
    gridUntangleBadFaceParameters(grid);
    gridRobustProject(grid); 
    // return 1;
  }

  gridSetCostConstraint(grid,
			gridCOST_CNST_VALID  | 
			gridCOST_CNST_VOLUME | 
                        gridCOST_CNST_AREAUV );

  printf("Spacing reset.\n");
  gridResetSpacing(grid);

  printf("cost const %d.\n",gridCostConstraint(grid));
  invalid=gridNumberOfInvalidCells(grid);printf("invalid %d.\n",invalid);
  for (cycle=0;cycle<3&&invalid>0;cycle++){
    printf("edge swapping grid...\n");gridSwap(grid,-1.0);
    invalid=gridNumberOfInvalidCells(grid);printf("invalid %d.\n",invalid);
  }
  for (cycle=0;cycle<5&&invalid>0;cycle++){
    STATUS;
    printf("edge collapse...\n");gridCollapseInvalidCells(grid);
    invalid=gridNumberOfInvalidCells(grid);printf("invalid %d.\n",invalid);
    STATUS;
    printf("edge swapping grid...\n");gridSwap(grid,-1.0);
    invalid=gridNumberOfInvalidCells(grid);printf("invalid %d.\n",invalid);
    STATUS;
    printf("simplex node...\n");gridSmoothInvalidCellNodes(grid);
    invalid=gridNumberOfInvalidCells(grid);printf("invalid %d.\n",invalid);
  }
  if (invalid>0) {
    gridWriteTecplotCurvedGeom(grid,"invalid.t");
    gridWriteTecplotInvalid(grid,"invalid.t");
  }else{
    // gridJacVolRatio(grid);
    STATUS;
    printf("edge swapping grid...\n");gridSwap(grid,-1.0);
    STATUS;
    printf("node smoothing grid...\n");gridSmooth(grid,-1.0,-1.0);
    STATUS;
    printf("edge swapping grid...\n");gridSwap(grid,-1.0);
    STATUS;
    printf("node smoothing grid...\n");gridSmooth(grid,-1.0,-1.0);
    STATUS;
    // gridJacVolRatio(grid);

    sprintf(filename, "%s_midnodes.t", outputProject );
    gridWriteTecplotCurvedGeom(grid,filename);

    printf("writing output project %s\n",outputProject);
    gridSavePart( grid, outputProject );

    sprintf(filename,"%s.bent", outputProject );
    printf("dumping curved Tetrahedral sides to %s\n",filename);
    gridDumpBentEdgesForPX(grid,2,filename);
  }
  printf("Done.\n");
  return 0;
}

