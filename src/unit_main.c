
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

void bl_metric_flat(Grid *grid, double h0) {
  int node;
  double xyz[3];
  double dx[3] = {1.0,0.0,0.0};
  double dy[3] = {0.0,1.0,0.0};
  double dz[3] = {0.0,0.0,1.0};
  double hx,hy,hz;
  double y0;
  for(node=0;node<gridMaxNode(grid);node++){
    if (grid==gridNodeXYZ(grid,node,xyz)) {
      hx=0.1; hy=0.1; hz=0.1;
      //      y0 = 0.5-ABS(0.5*(xyz[0]-0.5)*(xyz[2]-0.5));
      y0 = 0.5;
      hy = ABS(xyz[1]-y0)/0.5;
      hy = MIN(1.0,hy);
      hy = 0.1*((1.0-h0)*hy+h0);
      gridSetMapWithSpacingVectors(grid,node,
				   dx,dy,dz,hx,hy,hz);
    }
  }   
}

void bl_metric_sphere(Grid *grid, double h0) {
  int node;
  double xyz[3];
  double x,y,z;
  double x0=-0.5;
  double y0=-0.5;
  double z0=-0.5;
  double dx[3] = {1.0,0.0,0.0};
  double dy[3] = {0.0,1.0,0.0};
  double dz[3] = {0.0,0.0,1.0};
  double hx,hy,hz;
  double h1 = 0.1;
  double r0 = sqrt(3.0);
  double r, rp;
  for(node=0;node<gridMaxNode(grid);node++){
    if (grid==gridNodeXYZ(grid,node,xyz)) {
      x = xyz[0] - x0;
      y = xyz[1] - y0;
      z = xyz[2] - z0;
      r  = sqrt(x*x+y*y+z*z);
      dx[0] = x/r;
      dx[1] = y/r;
      dx[2] = z/r;
      rp = sqrt(x*x+y*y);
      dy[0] = -y/rp;
      dy[1] = x/rp;
      dy[2] = 0;
      gridCrossProduct(dx,dy,dz);
      gridVectorNormalize(dz);
      hx=0.1; hy=0.1; hz=0.1;
      r = ABS(r - r0)/0.5;
      r = MIN(1.0,r);
      hx = 0.1*((1.0-h0)*r+h0);
      gridSetMapWithSpacingVectors(grid,node,
				   dx,dy,dz,hx,hy,hz);
    }
  }   
}

void bl_metric(Grid *grid, double h0) {
  bl_metric_flat(grid, h0);
}

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
  int orig_cost; \
  bl_metric(grid,h0); \
  PRINT_STATUS; \
  orig_cost = gridCostFunction(grid); \
  gridSetCostFunction(grid, gridCOST_FCN_CONFORMITY ); \
   PRINT_STATUS; \
   gridSetCostFunction(grid, orig_cost ); \
}

static int hist_index = 0;

Grid *gridHistogram( Grid *grid, char *filename ) 
{
  int nodes[4];
  int edge, edgeId;
  int face, faceId;
  int cell;
  double length, volume;

  int bins[200];
  int bin;

  char fn[1024];
  FILE *file;
  
  if (NULL == filename) {
    hist_index++;
    sprintf( fn, "histogram%04d.m", hist_index );
    printf("dump %s\n",fn);
    file = fopen( fn, "w" );
  }else{
    file = fopen( filename, "w" );
  }

  for (bin=0;bin<200;bin++) bins[bin] = 0;
  gridCreateConn(grid);
  for ( edge = 0 ; edge < gridNConn(grid) ; edge++ ) {
    gridConn2Node( grid, edge, nodes );
    length = gridEdgeRatio(grid, nodes[0], nodes[1]);
    bin = (int)(log10(length)*100.0)+100;
    bin = MAX(0,bin);
    bin = MIN(bin,199);
    bins[bin]++;
  }
  gridEraseConn(grid);

  fprintf( file, "edge_length = [\n" );
  for (bin=0;bin<200;bin++) 
    fprintf( file, "%f %d\n", (((double)bin)-100.0)/100.0, bins[bin] );
  fprintf( file, "];\n" );

  for (bin=0;bin<200;bin++) bins[bin] = 0;
  for ( cell = 0 ; cell < gridMaxCell(grid) ; cell++ ) {
    if ( grid == gridCell( grid, cell, nodes ) ) {
      volume = gridVolume(grid, nodes);
      bin = (int)(-log10(volume)*10.0);
      bin = MAX(0,bin);
      bin = MIN(bin,199);
      bins[bin]++;
    }
  }
  fprintf( file, "tet_volume = [\n" );
  for (bin=199;bin>=0;bin--) 
    fprintf( file, "%f %d\n", -(((double)bin))/10.0, bins[bin] );
  fprintf( file, "];\n" );

  fprintf( file, "[x,y]=bar(edge_length(:,1), edge_length(:,2));plot(x,y,';Edge Length;');\n" );
  fprintf( file, "print -deps histogram_edge.eps\n" );
  fprintf( file, "closeplot\n" );
  fprintf( file, "[x,y]=bar(tet_volume(:,1), tet_volume(:,2));plot(x,y,';Tetrahedral Volume;');\n" );
  fprintf( file, "print -deps histogram_vol.eps\n" );
  fprintf( file, "closeplot\n" );

  fclose(file);

  return grid;
}


void relax_grid(Grid *grid, double h0) {
  int i;
  int iteration;
  int iterations = 10;
  double ratioSplit, ratioCollapse;
  double len0,len1;

  ratioCollapse = 0.3;
  ratioSplit    = 1.0;      

  STATUS;

  gridHistogram( grid, NULL );
  printf("edge swapping grid...\n");gridSwap(grid,0.9);
  STATUS;

  gridHistogram( grid, NULL );
  gridAdapt(grid, ratioCollapse, ratioSplit);
  STATUS;
    
  gridEdgeRatioRange(grid,&len0,&len1);
  for ( iteration=0; 
	(iteration<iterations) && len0 > 1.000001 ; 
	iteration++){
      
    for (i=0;i<1;i++){
      gridHistogram( grid, NULL );
      printf("edge swapping grid...\n");gridSwap(grid,0.9);
      STATUS;
      gridHistogram( grid, NULL );
      printf("node smoothin grid...\n");gridSmooth(grid,0.9,0.5);
      STATUS;
    }

    gridHistogram( grid, NULL );
    gridAdapt(grid, ratioCollapse, ratioSplit);
    STATUS;
    
    gridEdgeRatioRange(grid,&len0,&len1);
    
  }
  
  gridHistogram( grid, NULL );
}

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{
  Grid *grid;
  char modeler[256];
  char project[256];
  char ref_input[256];
  char ref_output[256];
  double h0 = 1.0;
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
    } else if( strcmp(argv[i],"-s") == 0 ) {
      i++; h0 = atof(argv[i]);
      printf("-s argument %d: %f\n",i, h0);
    } else if( strcmp(argv[i],"-h") == 0 ) {
      printf("Usage: flag value pairs:\n");
      printf(" -p input project name\n");
      printf(" -r input ref name\n");
      printf(" -o output ref name\n");
      printf(" -s bl height\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(strcmp(project,"")==0) sprintf( project, "box1" );
  if(strcmp(ref_output,"")==0) sprintf( ref_output, "box_out" );

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

  gridSetCostConstraint(grid,
			gridCOST_CNST_VOLUME | 
                        gridCOST_CNST_AREAUV );

  STATUS;

  gridSetMinInsertCost( grid, 1.0e-5 );
  gridSetMinSurfaceSmoothCost( grid, 1.0e-2 );

  DUMP_TEC;
  gridExportFAST( grid, "grid_orig.fgrid" );

  h0 = 1.0;
  relax_grid(grid,h0);
  DUMP_TEC;
  gridExportFAST( grid, "grid_h1000.fgrid" );

  h0 = 0.1;
  relax_grid(grid,h0);
  DUMP_TEC;
  gridExportFAST( grid, "grid_h0100.fgrid" );

  h0 = 0.01;
  relax_grid(grid,h0);
  DUMP_TEC;
  gridExportFAST( grid, "grid_h0010.fgrid" );

  h0 = 0.001;
  relax_grid(grid,h0);
  DUMP_TEC;
  gridExportFAST( grid, "grid_h0001.fgrid" );

  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: modifed grid does not have right handed boundaries\n");
  

  printf("writing output ref %s\n",ref_output);
  gridExportRef( grid, ref_output );

  printf("Done.\n");
  
  return 0;
}

