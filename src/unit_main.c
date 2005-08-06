
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

int gridNumberOfInvalidCells(Grid *grid)
{
  int cell, nodes[4];
  int invalid;
  invalid = 0;

  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, nodes)) {
      if ( -0.5 > gridAR(grid,nodes) ) {
	invalid++;
      }
    }
  }
  return invalid;
}

Grid *gridSmoothInvalidCellNodes(Grid *grid)
{
  int cell, nodes[4], i, node;

  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, nodes)) {
      if ( -0.5 > gridAR(grid,nodes) ) {
	for (i=0;i<4;i++) {
	  node = nodes[i];
	  if (!gridGeometryFace(grid,node)) {
	    gridSmartVolumeLaplacian( grid, node );	    
	    gridSmoothNodeVolumeSimplex(grid, node);
	  }
	}
      }
    }
  }
  return grid;
}

static int side2node0[] = {0, 0, 0, 1, 1, 2};
static int side2node1[] = {1, 2, 3, 2, 3, 3};

#define NEIGHBORDEG (500)

Grid *gridCollapseInvalidCells(Grid *grid)
{
  int cell, invalidnodes[4], nodes[4];
  int i, nlist, look, nodelist[NEIGHBORDEG];
  int corner, pivot, side, node0, node1;
  AdjIterator it;
  GridBool fixed, looking;

  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, invalidnodes)) {
      if ( -0.5 > gridAR(grid,invalidnodes) ) {
	fixed = FALSE;
	for (corner=0;corner<4 && !fixed;corner++) {

	  pivot = invalidnodes[corner];

	  nlist =0;
	  for ( it = adjFirst(gridCellAdj(grid),pivot); 
		adjValid(it); 
		it = adjNext(it) ){
	    gridCell(grid, adjItem(it), nodes);
	    for (i=0;i<4;i++) {
	      if (pivot != nodes[i]) {
		looking = (nlist<=NEIGHBORDEG);
		look = 0;
		for (look=0;look<nlist && looking ; look++){
		  looking = (nodelist[look] != nodes[i]);
		}
		if (looking && nlist<=NEIGHBORDEG){
		  nodelist[nlist] = nodes[i];
		  nlist++;
		}
	      }
	    }
	  }    

	  for (side=0;side<nlist && !fixed;side++) {
	    node0 = pivot;
	    node1 = nodelist[side];
	    if( (grid == gridCollapseEdge(grid, NULL, node0, node1, 0.00)) ||
		(grid == gridCollapseEdge(grid, NULL, node0, node1, 1.00)) ||
		(grid == gridCollapseEdge(grid, NULL, node0, node1, 0.50)) ) {
	      fixed = TRUE;
	    }
	  }

	}
      }
    }
  }

  return grid;
}

Grid *gridJacVolRatio(Grid *grid)
{
  int cell, nodes[4];
  double volume, jacobian, ratio;

  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, nodes)) {
      volume = gridVolume(grid, nodes );
      jacobian = gridMinCellJacDet2(grid, nodes );
      ratio = jacobian / volume / 6.0;
      if (ratio<0.9) printf("%10d %20.12e %10.8f\n",cell,volume,ratio);
    }
  }
  return grid;
}

Grid *gridUntangleBadFaceParameters(Grid *grid)
{
  int face, nodes[3], faceId;
  int node, *hits;
  
  hits = (int *)malloc( gridMaxNode(grid) * sizeof(int) );
  for (node=0;node<gridMaxNode(grid); node++) hits[node]=0;

  for (face=0;face<gridMaxFace(grid);face++) {
    if (grid == gridFace(grid,face,nodes,&faceId) ) {
      if ( 1.0e-14 > gridFaceAreaUV(grid, face) ) {
	hits[nodes[0]]++; hits[nodes[1]]++; hits[nodes[2]]++;
      }
    }
  }

  for (node=0;node<gridMaxNode(grid); node++) {
    if (hits[node] > 2) {
      printf("untangling node%10d with hits%3d\n",node,hits[node]);
      gridSmoothNodeFaceAreaUV(grid, node );
    }
  }  
  return grid;
}


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
  char filename[256];
  int i, j, oldSize, newSize;
  double ratio=0.6;
  double spacing = 1.0/3.0;
  double Zcommand = -1.0;
  double cyl = -1.0;
  int node;
  double Zspacing;
  double xyz[3];
  double minAR=-1.0;
  double ratioSplit, ratioCollapse;
  int EdgeBasedCycles = EMPTY;
  GridBool validate = FALSE;
  GridBool tecplotOutput = FALSE;
  GridBool LeadingEdgeBG = FALSE;
  double LeadingEdgeScale = 1.0;
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
    } else if( strcmp(argv[i],"-z") == 0 ) {
      i++; Zcommand = atof(argv[i]);
      printf("-z argument %d: %f\n",i, Zcommand);
    } else if( strcmp(argv[i],"-c") == 0 ) {
      i++; cyl = atof(argv[i]);
      printf("-c argument %d: %f\n",i, cyl);
    } else if( strcmp(argv[i],"-e") == 0 ) {
      i++; EdgeBasedCycles = atoi(argv[i]);
      printf("-e argument %d: %d\n",i, EdgeBasedCycles);
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
    } else if( strcmp(argv[i],"--validate") == 0 ) {
      validate = TRUE;
      printf("--validate argument %d\n",i);
    } else if( strcmp(argv[i],"-le") == 0 ) {
      LeadingEdgeBG = TRUE;
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
      printf(" -r initial edge length ratio for adapt\n");
      printf(" -s uniform grid size\n");
      printf(" -z linearly vary spacing to this in z dir\n");
      printf(" -c cylinder spacing\n");
      printf(" -e Number of Edge Based Operators Adaptation Cycles\n");
      printf(" -v freeze cells with small aspect ratio (viscous)\n");
      printf(" -f freeze nodes in this .lines file\n");
      printf(" -n max number of nodes in grid\n");
      printf(" -t write tecplot zones durring adaptation\n");
      printf(" -le scale leading edge background grid\n");
      printf(" --validate give grid valid cost constraints\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(strcmp(modeler,"")==0)       sprintf(modeler,"FELISA" );

  if(strcmp(project,"")==0)       sprintf(project,"../test/le" );
  LeadingEdgeBG = TRUE;

  if (validate) {
    if(strcmp(outputProject,"")==0) sprintf(outputProject,"%s_val", project );
  }else{
    if(strcmp(outputProject,"")==0) sprintf(outputProject,"%s_out", project );
  }
  if(strcmp(outputFAST,"")==0)    sprintf(outputFAST,"%s.fgrid",outputProject);

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

  if (EdgeBasedCycles!=EMPTY)gridSetCostFunction(grid,gridCOST_FCN_EDGE_LENGTH);

  gridSetCostConstraint(grid,
			gridCOST_CNST_VOLUME | 
                        gridCOST_CNST_AREAUV );
  //			gridCOST_CNST_VALID  |

  gridConstrainSurfaceNode(grid);

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
  printf("spacing set to constant %f with %f z variation.\n",spacing,Zcommand);
  for (node=0;node<gridMaxNode(grid);node++)
    if (grid==gridNodeXYZ(grid,node,xyz)) {
      Zspacing = spacing;
      if (Zcommand>0.0) 
	Zspacing = ABS(2.0*xyz[2]-1.0)*spacing 
	  + (1.0-ABS(2.0*xyz[2]-1.0))*Zcommand ;
      gridSetMap(grid, node,
		 1.0/spacing/spacing, 0.0, 0.0, 
		 1.0/spacing/spacing, 0.0, 
		 1.0/Zspacing/Zspacing);
    }
  if (cyl>0.0) {
    double radius, theta;
    double rSpace, tSpace, ySpace;
    double normal[3], tangent[3], third[3];
    printf("spacing set to cylinder with %f in normal direction.\n",cyl);
    for (node=0;node<gridMaxNode(grid);node++) {
      if (grid==gridNodeXYZ(grid,node,xyz)) {
	radius = sqrt(xyz[0]*xyz[0]+xyz[2]*xyz[2]);
	if (radius<1.0) radius = 2.0-radius;
	theta = atan2(xyz[2],xyz[0]);
	rSpace = cyl*radius*radius*radius;
	tSpace = 0.25*radius;
	ySpace = 0.25;
	normal[0]=normal[1]=normal[2]=0;
	tangent[0]=tangent[1]=tangent[2]=0;
	third[0]=third[1]=third[2]=0;
	normal[0]=cos(theta);
	normal[2]=sin(theta);
	tangent[0]=-sin(theta);
	tangent[2]=cos(theta);
	third[1]=1.0;
#ifdef ECHO_SPACING
	printf("X%6.2f Y%6.2f Z%6.2f\n",xyz[0],xyz[1],xyz[2]);
	printf("N%6.2f N%6.2f N%6.2f\n",normal[0],normal[1],normal[2]);
	printf("T%6.2f T%6.2f T%6.2f\n",tangent[0],tangent[1],tangent[2]);
	printf("Y%6.2f Y%6.2f Y%6.2f\n",third[0],third[1],third[2]);
	printf("\n");
#endif
	gridSetMapWithSpacingVectors(grid, node,
				     normal, tangent, third, 
				     rSpace, tSpace, ySpace);
      }
    }    
  }
  if (LeadingEdgeBG) {
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
  
  /*
    if (grid==gridRobustProject(grid)) {
    printf("projected grid to test params.\n");
    }else{
    printf("could not project grid. stop.\n");
    return 1;
    }
    gridUntangleBadFaceParameters(grid);
  */
  
  if (validate) {
    int cycle, invalid;
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
      printf("writing output project %s\n",outputProject);
      gridSavePart( grid, outputProject );
      printf("writing output FAST file %s\n",outputFAST);
      gridExportFAST( grid, outputFAST );
    }
    printf("Done.\n");
    return 0;
  }

  if (EMPTY!=EdgeBasedCycles) {
    sprintf(filename,"%s_surface.t",project);
    gridWriteTecplotSurfaceGeom(grid,filename); gridCloseTecplotGeomFile(grid);
    STATUS;
    printf("edge swapping grid...\n");gridSwap(grid,0.9);
    STATUS;
    for (j=0;j<EdgeBasedCycles;j++){
      printf("start edge based cycle %d\n",j);
      gridAdaptBasedOnConnRankings(grid);
      STATUS;
    }
    printf("edge swapping grid...\n");gridSwap(grid,0.9);
    STATUS;
    sprintf(filename,"%s_surface.t",outputProject);
    gridWriteTecplotSurfaceGeom(grid,filename); gridCloseTecplotGeomFile(grid);
    printf("writing output project %s\n",outputProject);
    gridSavePart( grid, outputProject );
    printf("Done.\n");
    return 0;
  }

  for (i=0;i<1;i++){
    printf("edge swapping grid...\n");gridSwap(grid,-1.0);
    STATUS;
    printf("node smoothing grid...\n");gridSmooth(grid,-1.0,-1.0);
  }
  STATUS;

  oldSize = 1;
  newSize = gridNNode(grid);
  jmax = 8;
  for ( j=0; 
	(j<jmax) ;//&& (((double)ABS(newSize-oldSize)/(double)oldSize)>0.01);
	j++){

    ratioCollapse = 0.4;
    ratioSplit   = 1.0;
    printf("adapt, ratio %4.2f, collapse limit %8.5f, refine limit %10.5f\n",
	   ratio, ratioCollapse, ratioSplit );
    if (gridCostConstraint(grid)&gridCOST_CNST_VALID) {
      if (grid != gridAdaptLongShortCurved(grid,ratioCollapse,ratioSplit,TRUE)){
	gridWriteTecplotCurvedGeom(grid, NULL );
	gridWriteTecplotSurfaceGeom(grid,NULL);
	return 1;
      }
    }else{
      if (grid != gridAdaptLongShortLinear(grid,ratioCollapse,ratioSplit,TRUE)){
	gridWriteTecplotCurvedGeom(grid, NULL );
	gridWriteTecplotSurfaceGeom(grid,NULL);
	return 1;
      }
    }
    oldSize = newSize;
    newSize = gridNNode(grid) ;
    printf("%02d new size: %d nodes %d faces %d cells %d edge elements.\n",
	   j, gridNNode(grid),gridNFace(grid),gridNCell(grid),gridNEdge(grid));
    STATUS;

    if (!gridSurfaceNodeConstrained(grid)){
      GridMove *gm;
      double minVolume;
      printf("Calling GridMove to project nodes...\n");
      gm = gridmoveCreate(grid);
      gridmoveProjectionDisplacements(gm);
      gridmoveRelaxation(gm,gridmoveELASTIC_SCHEME,1,100);
      gridmoveApplyDisplacements(gm);
      gridmoveFree(gm);
      STATUS; minVolume = gridMinVolume(grid);
      while (0.0>=minVolume) {
	printf("relax neg faces...\n");
	gridParallelRelaxNegativeFaceAreaUV(grid,FALSE);
	printf("relax neg cells...\n");gridRelaxNegativeCells(grid,TRUE);
	printf("edge swapping grid...\n");gridSwap(grid,1.0);
	STATUS; minVolume = gridMinVolume(grid);
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

