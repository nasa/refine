
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

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  Grid *grid;
  int iteration;
  int iterations = 8;
  char project[256];
  char outputProject[256];
  char filename[256];
  double ratioSplit, ratioCollapse;
  GridBool tecplotOutput = FALSE;
  double LeadingEdgeScale = 1.0;
  int iview = 0;
  int maxnode = 50000;
  char modeler[81];

  int i;

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
    } else if( strcmp(argv[i],"-le") == 0 ) {
      i++; LeadingEdgeScale = atof(argv[i]);
      printf("-le argument %d: %f\n",i,LeadingEdgeScale);
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
      printf(" -le scale leading edge background grid\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(strcmp(modeler,"")==0)       sprintf(modeler,"FELISA" );

  if(strcmp(project,"")==0)       sprintf(project,"../test/le" );

  if(strcmp(outputProject,"")==0) sprintf(outputProject,"%s_out", project );

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
  }

  gridSetCostConstraint(grid,
			gridCOST_CNST_VOLUME | 
                        gridCOST_CNST_AREAUV );

  gridSetMinInsertCost( grid, 0.1 );
  gridSetMinSurfaceSmoothCost( grid, 0.1 );

  // gridConstrainSurfaceNode(grid);

  printf("Spacing reset.\n");
  gridResetSpacing(grid);

  if (TRUE) {
    int node;
    double xyz[3];
    double centerX;
    double radius, theta;
    double rSpace, tSpace, ySpace;
    double normal[3], tangent[3], third[3];
    printf("spacing set to Leading Edge.\n");
    for (node=0;node<gridMaxNode(grid);node++) {
      if (grid==gridNodeXYZ(grid,node,xyz)) {
	centerX = 0.05;
	radius = sqrt((xyz[0]-centerX)*(xyz[0]-centerX)+xyz[2]*xyz[2]);
	theta = atan2(xyz[2],xyz[0]-centerX);
	rSpace = 0.2*pow(MIN(4.0*radius,1.0),LeadingEdgeScale);
	tSpace = 0.05+0.1*radius;
	ySpace = 0.25;
	normal[0]=normal[1]=normal[2]=0;
	tangent[0]=tangent[1]=tangent[2]=0;
	third[0]=third[1]=third[2]=0;
	normal[0]=cos(theta);
	normal[2]=sin(theta);
	tangent[0]=-sin(theta);
	tangent[2]=cos(theta);
	third[1]=1.0;
	gridSetMapWithSpacingVectors(grid, node,
				     normal, tangent, third, 
				     rSpace, tSpace, ySpace);
      }
    }    
  }
  
  for (i=0;i<1;i++){
    printf("edge swapping grid...\n");gridSwap(grid,-1.0);
    STATUS;
    printf("node smoothing grid...\n");gridSmooth(grid,-1.0,-1.0);
  }
  STATUS;

  for ( iteration=0; (iteration<iterations) ; iteration++){

    ratioCollapse = 0.4;
    ratioSplit    = 1.0;

    gridAdapt(grid, ratioCollapse, ratioSplit);

    printf("%02d new size: %d nodes %d faces %d cells %d edge elements.\n",
	   iteration,
	   gridNNode(grid),gridNFace(grid),gridNCell(grid),gridNEdge(grid));
    STATUS;

    if (!gridSurfaceNodeConstrained(grid)){
      GridMove *gm;
      double minArea, minVolume;
      int untangling_steps;
      printf("Calling GridMove to project nodes...\n");
      gm = gridmoveCreate(grid);
      gridmoveProjectionDisplacements(gm);
      gridmoveRelaxation(gm,gridmoveELASTIC_SCHEME,1,2000);
      gridmoveApplyDisplacements(gm);
      gridmoveFree(gm);
      minArea = gridMinGridFaceAreaUV(grid); untangling_steps = 0;
      while (minArea < 1.0e-12) { // bump this up?
	printf("min face UV area %e\n",minArea);
	printf("relax neg faces...\n");
	gridParallelRelaxNegativeFaceAreaUV(grid,FALSE);
	minArea = gridMinGridFaceAreaUV(grid); untangling_steps++;
	if (untangling_steps >3) return 1;
      }
      STATUS; minVolume = gridMinVolume(grid); untangling_steps = 0;
      while (0.0>=minVolume) {
	printf("relax neg cells...\n");gridRelaxNegativeCells(grid,TRUE);
	printf("edge swapping grid...\n");gridSwap(grid,1.0);
	STATUS; minVolume = gridMinVolume(grid); untangling_steps++;
	if (untangling_steps >3) return 1;
      }
    }

    for (i=0;i<1;i++){
      printf("edge swapping grid...\n");gridSwap(grid,-1.0);
      STATUS;
      printf("node smoothing grid...\n");gridSmooth(grid,-1.0,-1.0);
    }
    STATUS;
  }

  for (i=0;i<2;i++){
    printf("edge swapping grid...\n");gridSwap(grid,-1.0);
    STATUS;
    printf("node smoothing grid...\n");gridSmooth(grid,-1.0,-1.0);
    STATUS;
  }

  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: modifed grid does not have right handed boundaries\n");

  printf("writing output project %s\n",outputProject);
  gridSavePart( grid, outputProject );

  printf("Done.\n");

  return 0;
}

