
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
#include "gridmove.h"
#include "layer.h"

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
  DUMP_TEC; \
  PRINT_STATUS; \
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

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  Grid *grid;
  int iteration;
  int iterations = 1;
  char modeler[256];
  char project[256];
  char ref_input[256];
  char ref_output[256];
  double ratioSplit, ratioCollapse;
  GridBool tecplotOutput = TRUE;

  int i;
  int iview=0;

  sprintf( modeler,    "Unknown" );
  sprintf( project,    "" );
  sprintf( ref_input,  "" );
  sprintf( ref_output, "" );

  i = 1;
  while( i < argc ) {
    if( strcmp(argv[i],"-p") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-p argument %d: %s\n",i, project);
    } else if( strcmp(argv[i],"-r") == 0 ) {
      i++; sprintf( ref_input, "%s", argv[i] );
      printf("-r argument %d: %s\n",i, ref_input);
    } else if( strcmp(argv[i],"-o") == 0 ) {
      i++; sprintf( ref_output, "%s", argv[i] );
      printf("-o argument %d: %s\n",i, ref_output);
    } else if( strcmp(argv[i],"-h") == 0 ) {
      printf("Usage: flag value pairs:\n");
      printf(" -p input project name\n");
      printf(" -r input ref name\n");
      printf(" -o output ref name\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(strcmp(ref_output,"")==0){
    printf("no output specified.\n");
    printf("Done.\n");  
    return 0;
  }

  if(!(strcmp(ref_input,"")==0)) {
    printf("running ref %s\n",ref_input);
    grid = gridImportRef( ref_input );
    if(strcmp(project,"")==0){
      printf("no geom specified.\n");
      printf("Done.\n");  
      return 0;
    }
    printf("loading geometry %s\n",project);
    gridGeomStartOnly( grid, project );
  }else{
    if(!(strcmp(project,"")==0)) {
      printf("running project %s\n",project);
      grid = gridLoadPart( modeler, project, 50000 );
    }
  }

  printf("restart grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  if (!gridRightHandedBoundary(grid)) {
    printf("ERROR: loaded part does not have right handed boundaries\n");
    return 1;
  }

  printf("Spacing reset.\n");
  gridResetSpacing(grid);

  if(TRUE) {
    int node;
    double xyz[3];
    double dx[3] = {1.0,0.0,0.0};
    double dy[3] = {0.0,1.0,0.0};
    double dz[3] = {0.0,0.0,1.0};
    double hx,hy,hz;
    for(node=0;node<gridMaxNode(grid);node++){
      if (grid==gridNodeXYZ(grid,node,xyz)) {
	hx=0.1; hy=0.1; hz=0.1;
	hz = ABS(xyz[2]-0.5)/0.5;
	hz = 0.1*(0.9*hz+0.1);
	gridSetMapWithSpacingVectors(grid,node,
				     dx,dy,dz,hx,hy,hz);
      }
    }
      
  }

  gridSetCostConstraint(grid,
			gridCOST_CNST_VOLUME | 
                        gridCOST_CNST_AREAUV );

  STATUS;

  gridCacheCurrentGridAndMap(grid);

  ratioCollapse = 0.3;
  ratioSplit    = 1.5;
      
  gridSetMinInsertCost( grid, 1.0e-5 );
  gridSetMinSurfaceSmoothCost( grid, 1.0e-5 );
    
  for (i=0;i<1;i++){
    printf("edge swapping grid...\n");gridSwap(grid,-1.0);
    STATUS;
  }
  for ( iteration=0; (iteration<iterations) ; iteration++){
      
    gridAdapt(grid, ratioCollapse, ratioSplit);

    for (i=0;i<1;i++){
      printf("edge swapping grid...\n");gridSwap(grid,-1.0);
      STATUS;
    }
  }
  
  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: modifed grid does not have right handed boundaries\n");
  
  gridExportFAST( grid, NULL );

  printf("writing output ref %s\n",ref_output);
  gridExportRef( grid, ref_output );

  if (tecplotOutput) gridWriteTecplotSurfaceGeom(grid, NULL );

  printf("Done.\n");
  
  return 0;
}

