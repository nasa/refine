
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
#include "gridedger.h"
#include "gridfacer.h"

#define PRINT_STATUS {double l0,l1;gridEdgeRatioRange(grid,&l0,&l1);printf("Len %12.5e %12.5e AR %8.6f MR %8.6f Vol %10.6e\n", l0,l1, gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid)); fflush(stdout);}

#define DUMP_TEC if (tecplotOutput) { \
 iview++;printf("Frame %d\n",iview);\
 gridWriteTecplotSurfaceGeom(grid,NULL); \
  {int cell, nodes[4]; \
    for (cell=0;cell<gridMaxCell(grid);cell++) \
      if (grid==gridCell(grid, cell, nodes)) { \
        if ( -0.5 > gridAR(grid,nodes) ) \
          gridWriteTecplotCellGeom(grid,nodes,NULL,NULL); \
   }}; }

#define STATUS { \
  if (LeadingEdgeBG && update_spacing_durring_status) \
    leading_edge_spacing(grid, LeadingEdgeScale); \
  DUMP_TEC; \
  PRINT_STATUS; \
}

void leading_edge_spacing(Grid *grid, double LeadingEdgeScale ) {
  int node;
  double xyz[3];
  double centerX, nose_radius;
  double radius, theta;
  double rate;
  double rSpace, tSpace, ySpace;
  double normal[3], tangent[3], third[3];
  for (node=0;node<gridMaxNode(grid);node++) {
    if (grid==gridNodeXYZ(grid,node,xyz)) {
      centerX = 0.03;
      nose_radius = 0.03;
      radius = sqrt((xyz[0]-centerX)*(xyz[0]-centerX)+xyz[2]*xyz[2]);
      theta = atan2(xyz[2],xyz[0]-centerX);
      rate = 1.2;
      tSpace = 0.10+0.1*radius;
      rSpace = MIN( tSpace,
		    LeadingEdgeScale * 
		    pow(rate,(radius-nose_radius)/LeadingEdgeScale) );
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
      
Grid *gridSplitProblemProjectionEdges(Grid *grid) {
  int conn, nodes[2];
  int parent, newnode;
  int split_edges;
  double min_insert_cost;
  min_insert_cost = gridMinInsertCost( grid );
  gridSetMinInsertCost( grid, -100.0 );
  gridCreateConn(grid);
  split_edges = 0;
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    parent = gridParentGeometry(grid, nodes[0], nodes[1] );
    if ( ( gridGeometryFace( grid, nodes[0] ) &&
	   gridGeometryFace( grid, nodes[1] ) &&
	   0 == parent  ) ||
	 ( gridGeometryEdge( grid, nodes[0] ) &&
	   gridGeometryEdge( grid, nodes[1] ) &&
	   0 < parent  ) ) {
      newnode = gridSplitEdgeRatio( grid, NULL, nodes[0], nodes[1], 0.5);
      split_edges++;
    }
  }
  printf("split %d edges\n",split_edges);
  gridEraseConn(grid);
  gridSetMinInsertCost( grid, min_insert_cost );
  return grid;
}

Grid *gridHistogram( Grid *grid, char *filename ) 
{
  int nodes[4];
  int edge, edgeId;
  int face, faceId;
  int cell;
  double length, area, volume;

  FILE *file;

  if (NULL == filename) {
    file = fopen( "histogram.m", "w" );
  }else{
    file = fopen( filename, "w" );
  }

  fprintf( file, "edge_length = [\n" );
  for ( edge = 0 ; edge < gridMaxEdge(grid) ; edge++ ) {
    if ( grid == gridEdge( grid, edge, nodes, &edgeId ) ) {
      length = gridEdgeRatio(grid, nodes[0], nodes[1]);
      fprintf( file, "%e\n", length );
    }
  }
  fprintf( file, "];\n" );

  fprintf( file, "triangle_uv_area = [\n" );
  for ( face = 0 ; face < gridMaxFace(grid) ; face++ ) {
    if ( grid == gridFace( grid, face, nodes, &faceId ) ) {
      if ( gridGeometryEdge( grid, nodes[0] ) ||
	   gridGeometryEdge( grid, nodes[1] ) ||
	   gridGeometryEdge( grid, nodes[2] ) ) {
	area = gridFaceAreaUV(grid, face);
	fprintf( file, "%e\n", area );
      }
    }
  }
  fprintf( file, "];\n" );

  fprintf( file, "tet_volume = [\n" );
  for ( cell = 0 ; cell < gridMaxCell(grid) ; cell++ ) {
    if ( grid == gridCell( grid, cell, nodes ) ) {
      if ( gridGeometryEdge( grid, nodes[0] ) ||
	   gridGeometryEdge( grid, nodes[1] ) ||
	   gridGeometryEdge( grid, nodes[2] ) ||
	   gridGeometryEdge( grid, nodes[3] ) ) {
	volume = gridVolume(grid, nodes);
	fprintf( file, "%e\n", volume );
      }
    }
  }
  fprintf( file, "];\n" );

  fclose(file);

  return grid;
}

Grid *gridPhase1(Grid *grid )
{
  int i;
  int edge, edgeId;
  GridEdger **ge;

  gridSetPhase(grid, 1);
  gridSetMinInsertCost( grid, -0.5 );
  gridConstrainSurfaceNode(grid);
  gridSplitProblemProjectionEdges(grid);
  gridUntangle(grid);
  
  ge = (GridEdger **)malloc( gridNGeomEdge(grid) * sizeof(GridEdger *) );
  for ( edge = 0 ; edge < gridNGeomEdge(grid) ; edge++ ) {
    edgeId = edge+1;
    ge[edge] = gridedgerCreate(grid,edgeId);
  }

#define FREE_ALL_GE(ge) { \
for ( edge = 0 ; edge < gridNGeomEdge(grid) ; edge++ ) { \
  gridedgerFree(ge[edge]); \
} \
free(ge); }

  for ( edge = 0 ; edge < gridNGeomEdge(grid) ; edge++ ) {
    edgeId = edge+1;
    if ( ge[edge] != gridedgerDiscretizeEvenly(ge[edge]) ) {
      printf("gridedgerDiscretizeEvenly failed for edge %d\n",edgeId);
      FREE_ALL_GE(ge);
      return NULL;
    }
    if ( ge[edge] != gridedgerInsert(ge[edge]) ) {
      printf("gridedgerInsert failed for edge %d\n",edgeId);
      FREE_ALL_GE(ge);
      return NULL;
    }
    if (grid != gridUntangle(grid) ) {
      printf("gridUntangle failed for edge %d\n",edgeId);
      FREE_ALL_GE(ge);
      return NULL;      
    }
  }

  for ( edge = 0 ; edge < gridNGeomEdge(grid) ; edge++ ) {
    edgeId = edge+1;
    if ( ge[edge] != gridedgerRemoveUnused(ge[edge]) ) {
      printf("gridedgerRemoveUnused failed for edge %d\n",edgeId);
      FREE_ALL_GE(ge);
      return NULL;
    }
    for ( i = 0 ; i < gridedgerUnusedNodes(ge[edge]) ; i++ ) {
      double xyz[3];
      gridNodeXYZ(grid, gridedgerUnusedNode(ge[edge],i), xyz);
      printf("edge%4d xyz %f %f %f\n",
	     edgeId, xyz[0], xyz[1], xyz[2] );
    }
  }
    
  FREE_ALL_GE(ge);

  return grid;
}

Grid *gridPhase2(Grid *grid )
{
  int i;
  GridFacer *gf;
  int faceId;
  double ratio0, ratio1;
    
  gridSplitProblemProjectionEdges(grid);
  gridSetPhase(grid, 2);
  gridSetMinInsertCost( grid, -0.5 );
  gridConstrainSurfaceNode(grid);

  for (faceId = 1; faceId <= gridNGeomFace(grid); faceId++) {
    gf = gridfacerCreate(grid,faceId);

    gridfacerTurnCameraOn(gf);

    for ( i = 0 ; i < 20 ; i++ ) {
      gridfacerSwap(gf);
      gridfacerSplitProblemProjectionEdges(gf);
      gridSetMinInsertCost( grid, -10.0 );
      if (gf != gridfacerSplit(gf)){
	printf("gridfacerSplit failed for face %d\n",faceId);
	gridfacerFree(gf);
	gridWriteTecplotSurfaceGeom(grid,NULL);
	{int cell, nodes[4];
	for (cell=0;cell<gridMaxCell(grid);cell++) 
	  if (grid==gridCell(grid, cell, nodes)) { 
	    if ( -0.5 > gridAR(grid,nodes) ) 
	      gridWriteTecplotCellGeom(grid,nodes,NULL,NULL);
	  }
	}
	return NULL;      
      }
      if (grid != gridUntangle(grid) ) {
	printf("gridUntangle failed for face %d\n",faceId);
	gridfacerFree(gf);
	return NULL;      
      }
      gridSetMinInsertCost( grid, -0.5 );
      gridfacerRatioRange(gf,&ratio0,&ratio1);
      printf("face%6d cycle%3d ratios%8.3f%8.3f\n",faceId,i,ratio0,ratio1);
      gridfacerCollapse(gf);
      gridfacerRatioRange(gf,&ratio0,&ratio1);
      printf("face%6d cycle%3d ratios%8.3f%8.3f\n",faceId,i,ratio0,ratio1);
      if (1.0 >= ratio0) break;
    }
      
    gridfacerFree(gf);

    if (1.0 < ratio0) {
      printf("Phase 2 failed for face %d\n",faceId);  
      return NULL;
    }
  }

  return grid;
}

Grid *gridPhase3(Grid *grid )
{
  int iteration;
  int iterations;
  double longest_edge, shortest_edge;

  gridSetPhase(grid, 3);

  gridSetMinInsertCost( grid, -0.5 );
  gridSetCostFunction( grid, gridCOST_FCN_EDGE_LENGTH );
  iterations = 20;
  longest_edge = 2.0;
  for ( iteration=0; 
	(iteration<iterations) && (longest_edge>1.000001); 
	iteration++){
    gridAdaptVolumeEdges(grid);
    gridEdgeRatioRange(grid,&longest_edge,&shortest_edge);
    printf("Length Range %12.5e %12.5e\n", longest_edge, shortest_edge);
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
  int iteration;
  int iterations = 8;
  char project[256];
  char outputProject[256];
  char filename[256];
  double ratioSplit, ratioCollapse;
  GridBool tecplotOutput = FALSE;
  GridBool edge_based = FALSE;
  int phase = 0;
  GridBool LeadingEdgeBG = TRUE;
  double LeadingEdgeScale = 1.0;
  double global_scale = 1.0;
  GridBool update_spacing_durring_status = FALSE;
  GridBool FAST_output = FALSE;
  int iview = 0;
  int maxnode = 50000;
  char modeler[81];

  int i;

  double minArea, minVolume;
  int untangling_steps;

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
   } else if( strcmp(argv[i],"-s") == 0 ) {
      i++; global_scale = atof(argv[i]);
      LeadingEdgeBG = FALSE;
      printf("-s argument %d: %f\n",i,global_scale);
   } else if( strcmp(argv[i],"-e") == 0 ) {
      edge_based = TRUE;
      printf("-e argument %d\n",i);
   } else if( strcmp(argv[i],"--phase") == 0 ) {
      i++; phase = atoi(argv[i]);
      printf("--phase argument %d: %d\n",i, phase);
    } else if( strcmp(argv[i],"--fast") == 0 ) {
      i++; FAST_output = TRUE;
      printf("--fast argument %d\n",i );
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
      printf(" -s global scale of current spacing (deactivates LE)\n");
      printf(" -e edge length only adaptation\n");
      printf(" --phase phase of adaptation to run\n");
      printf(" --fast output FAST ASCII file\n");
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

  printf("Spacing reset.\n");
  gridResetSpacing(grid);
  if (LeadingEdgeBG) leading_edge_spacing(grid, LeadingEdgeScale);
  for( i = 0 ; i < gridMaxNode(grid) ; i++ ) { 
    gridScaleSpacing( grid, i, global_scale );
  }

  gridSetCostConstraint(grid,
			gridCOST_CNST_VOLUME | 
                        gridCOST_CNST_AREAUV );

  if (4 == phase) {
    gridCacheCurrentGridAndMap(grid);
    if (grid!=gridPhase1(grid)) return 1;
    if (grid!=gridPhase2(grid)) return 1;
    if (grid!=gridPhase3(grid)) return 1;
    STATUS;
    printf("writing output project %s\n",outputProject);
    gridSavePart( grid, outputProject );
    sprintf(filename,"%s.fgrid", outputProject );
    printf("dumping FAST file to %s\n",filename);
    gridExportFAST(grid,filename);
    printf("Done.\n");
    return 0;
  } else {
    gridSetPhase(grid, phase);
  }

  ratioCollapse = 0.3;
  ratioSplit    = 1.0;
      
  if (edge_based) { /* if edge_based */
    gridSetMinInsertCost( grid, 1.0e-5 );
    gridSetMinSurfaceSmoothCost( grid, 1.0e-5 );
    
    gridSetCostFunction( grid, gridCOST_FCN_EDGE_LENGTH );

    for (i=0;i<1;i++){
      printf("edge swapping grid...\n");gridSwap(grid,-1.0);
      STATUS;
    }
    for ( iteration=0; (iteration<iterations) ; iteration++){
      
      gridAdapt(grid, ratioCollapse, ratioSplit);
      gridSplitProblemProjectionEdges( grid );
      
      printf("%02d new size: %d nodes %d faces %d cells %d edge elements.\n",
	     iteration,
	     gridNNode(grid),gridNFace(grid),gridNCell(grid),gridNEdge(grid));
      STATUS;
      
      if (!gridSurfaceNodeConstrained(grid)){
	GridMove *gm;
	printf("Calling GridMove to project nodes...\n");
	gm = gridmoveCreate(grid);
	gridmoveProjectionDisplacements(gm);
	gridmoveApplyDisplacements(gm);
	gridmoveFree(gm);
	if (grid != gridUntangle(grid)) return 1;
      }

      for (i=0;i<1;i++){
	printf("edge swapping grid...\n");gridSwap(grid,-1.0);
	STATUS;
      }
    }

  }else{ /* not edge_based */
    gridSetMinInsertCost( grid, -0.5 );
    gridSetMinSurfaceSmoothCost( grid, -1.5 );
    
    for ( iteration=0; (iteration<iterations) ; iteration++){
      
      gridAdapt(grid, ratioCollapse, ratioSplit);
      gridSplitProblemProjectionEdges( grid );
      
      printf("%02d new size: %d nodes %d faces %d cells %d edge elements.\n",
	     iteration,
	     gridNNode(grid),gridNFace(grid),gridNCell(grid),gridNEdge(grid));
      STATUS;

      if ( grid != gridWholesaleEvaluation(grid) ) {
	printf("evaluate points on surface...FAILED\n");
	tecplotOutput=TRUE;
	DUMP_TEC;
	return 1;
      }
      
      for (i=0;i<1;i++){
	if (gridEDGE_PHASE != gridPhase(grid)) {
	  printf("edge swapping grid...\n");gridSwap(grid,-1.0);
	  STATUS;
	}
	printf("node smoothing grid...\n");gridSmooth(grid,-1.0,-1.0);
	STATUS;
	if (grid != gridUntangle(grid)) return 1;
	STATUS;
      }
    }
    
  } /* end edge_based */
  
  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: modifed grid does not have right handed boundaries\n");
  
  printf("writing output project %s\n",outputProject);
  gridSavePart( grid, outputProject );
  
  if (FAST_output) {
    sprintf(filename,"%s.fgrid", outputProject );
    printf("dumping FAST file to %s\n",filename);
    gridExportFAST(grid,filename);
  }

  printf("Done.\n");
  
  return 0;
}

