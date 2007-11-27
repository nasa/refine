
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
#include "interp.h"
#include "plan.h"

static  GridBool tecplotOutput = TRUE;
static  int iview=0;

void interp_metric(Grid *grid) {
  int node;
  double xyz[3];
  double m[6];
  Interp *interp = gridInterp(grid);
  for(node=0;node<gridMaxNode(grid);node++){
    if (grid==gridNodeXYZ(grid,node,xyz)) {
      interpMetric(interp,xyz,m);
      gridSetMap(grid,node,
		 m[0], m[1], m[2],
		 m[3], m[4],
		 m[5]);
    }
  }
}

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
  if (FALSE) {
    int node;
    double xyz[3], map[6];
    for(node=0;node<gridMaxNode(grid);node++){
      if (grid==gridNodeXYZ(grid,node,xyz)) {
	gridMap( grid, node, map );
	printf(" %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
	       map[0], map[1], map[2], map[3], map[4], map[5]);
      }
    }
  }
}

//#define PRINT_STATUS {double l0,l1;gridEdgeRatioRange(grid,&l0,&l1);printf("Len %9.2e %9.2e AR %12.5e err %12.5e Vol %10.6e\n", l0,l1, gridMinThawedAR(grid),interpTotalError(grid), gridMinVolume(grid)); fflush(stdout);}
#define PRINT_STATUS {printf("AR %12.5e err %12.5e Vol %10.6e\n", gridMinThawedAR(grid),interpTotalError(grid), gridMinVolume(grid)); fflush(stdout);}

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
  PRINT_STATUS; \
  printf("%10d %12.5e %12.5e %%oct\n",gridNNode(grid),sqrt(gridMinThawedAR(grid)),interpTotalError(grid)); \
}

Grid *gridHistogram( Grid *grid, char *filename ) 
{
  int nodes[4];
  int cell;
  double cost;

  FILE *file;

  if (NULL == filename) {
    file = fopen( "histogram.m", "w" );
  }else{
    file = fopen( filename, "w" );
  }

  fprintf( file, "cost = [\n" );
  for ( cell = 0 ; cell < gridMaxCell(grid) ; cell++ ) {
    if ( grid == gridCell( grid, cell, nodes ) ) {
      cost = gridAR(grid, nodes);
      fprintf( file, "%e\n", 1.0/(1.0+sqrt(cost)) );
    }
  }
  fprintf( file, "];\n" );

  fclose(file);

  return grid;
}

void adapt3smooth(Grid *grid)
{
  Plan *plan;
  int node, ranking;
  double cost;
  int report;
  double min_cost;

  min_cost = 10;

  plan = planCreate( gridNNode(grid), MAX(gridNNode(grid)/10,1000) );
  for (node=0;node<gridMaxNode(grid);node++) {
    if (gridValidNode(grid, node)) {
      gridNodeAR(grid,node,&cost);
      if ( cost > min_cost ) planAddItemWithPriority( plan, node, cost );
    }
  }
  planDeriveRankingsFromPriorities( plan );
  report = 10; if (planSize(plan) > 100) report = planSize(plan)/10;
  for ( ranking=planSize(plan)-1; ranking>=0; ranking-- ) { 
    node = planItemWithThisRanking(plan,ranking);
    gridSmoothNode(grid, node, TRUE );
    if ( ranking/report*report==ranking ){
      gridNodeAR(grid,node,&cost);
      printf("rank %d cost %f\n",ranking,cost);
      fflush(stdout);
    }
  }
  planFree(plan);
}

void adapt3swap(Grid *grid)
{
  Plan *plan;
  int cell, nodes[4], ranking;
  double cost,cost1;
  int report;
  double min_cost;

  min_cost = 10;

  plan = planCreate( gridNCell(grid), MAX(gridNCell(grid)/10,1000) );
  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, nodes)) {
      cost = gridAR(grid,nodes);
      if ( cost > min_cost ) planAddItemWithPriority( plan, cell, cost );
    }
  }
  planDeriveRankingsFromPriorities( plan );
  report = 10; if (planSize(plan) > 100) report = planSize(plan)/10;
  for ( ranking=planSize(plan)-1; ranking>=0; ranking-- ) { 
    cell = planItemWithThisRanking(plan,ranking);
    if (grid==gridCell(grid, cell, nodes)) {
      cost = gridAR(grid,nodes);
      if ( cost > min_cost ) {
      if ( ( NULL != gridSwapFace( grid, NULL, nodes[1],nodes[2],nodes[3]) ) ||
	   ( NULL != gridSwapFace( grid, NULL, nodes[0],nodes[2],nodes[3]) ) ||
	   ( NULL != gridSwapFace( grid, NULL, nodes[0],nodes[1],nodes[3]) ) ||
	   ( NULL != gridSwapFace( grid, NULL, nodes[0],nodes[1],nodes[2]) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[1] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[2] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[0], nodes[3] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[1], nodes[2] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[1], nodes[3] ) ) ||
	   ( NULL != gridSwapEdge( grid, NULL, nodes[2], nodes[3] ) ) ) {
      }
      }
      if ( ranking/report*report==ranking ){
	printf("rank %d cost %f \n",ranking,cost);
	fflush(stdout);
      }
    }
  }
  planFree(plan);
}

void adapt3insert(Grid *grid)
{
  Plan *plan;
  int conn, nodes[2], ranking;
  double length,maxLength;
  int nnodeAdd, nnodeRemove;
  int report;
  double ratio;
  int newnode;
  double cost0, cost1;

  int i;

  maxLength = 1.5;

  gridCreateConn(grid);
  plan = planCreate( gridNConn(grid)/2, MAX(gridNConn(grid)/10,1000) );
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    length = gridEdgeRatio(grid,nodes[0],nodes[1]);
    if ( length >= maxLength ) planAddItemWithPriority( plan, conn, length );
  }
  planDeriveRankingsFromPriorities( plan );
  
  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (planSize(plan) > 100) report = planSize(plan)/10;

  for ( ranking=planSize(plan)-1; ranking>=0; ranking-- ) { 
    conn = planItemWithThisRanking(plan,ranking);
    if (ranking/report*report == ranking || ranking==planSize(plan)-1) {
      printf("adapt ranking%9d nnode%9d added%9d removed%9d err%6.2f\n",
	     ranking,gridNNode(grid),nnodeAdd,nnodeRemove,
	     planPriorityWithThisRanking(plan,ranking));
      fflush(stdout);
    }
    if (grid == gridConn2Node(grid,conn,nodes)){
      if ( gridCellEdge(grid, nodes[0], nodes[1]) &&
	   gridValidNode(grid, nodes[0]) && 
	   gridValidNode(grid, nodes[1]) && 
	   !gridNodeFrozen(grid, nodes[0]) &&
	   !gridNodeFrozen(grid, nodes[1]) ) {
	length = gridEdgeRatio(grid, nodes[0], nodes[1]);
	if (length >= maxLength) {
	  if (grid!=gridMakeGem(grid,nodes[0], nodes[1])) continue;
	  if (grid!=gridGemAR(grid,&cost0)) continue;
	  ratio = 0.5;
	  newnode = gridSplitEdgeRatio( grid, NULL,
					nodes[0], nodes[1], ratio );
	  if ( newnode != EMPTY ){
	    gridSmoothNode(grid, newnode, TRUE );
	    gridNodeAR(grid, newnode, &cost1 );
	    //printf("cost %15.7e %f %f \n",cost1-cost0,cost0,cost1);
	    if ( cost1 > cost0 ) {
	      gridCollapseEdge(grid, NULL, nodes[0], newnode, 0.0);
	    }else{
	      nnodeAdd++;
	    }
	  } 
	}
      }
    }
  }
  planFree(plan);
  gridEraseConn(grid);
}

GridBool splitit(Grid *grid, int en[2])
{
  double cost0,cost1;
  double ratio;
  int newnode;
  ratio = 0.5;
  if (grid!=gridMakeGem(grid,en[0], en[1])) return FALSE;
  if (grid!=gridGemAR(grid,&cost0)) return FALSE;
  newnode = gridSplitEdgeRatio( grid, NULL, en[0], en[1], ratio );
  if ( newnode != EMPTY ){
    gridSmoothNode(grid, newnode, TRUE );
    gridNodeAR(grid, newnode, &cost1 );
    if ( cost1 > cost0 ) {
      if (grid != gridCollapseEdge(grid, NULL, en[0], newnode, 0.0))
	printf("collapse failed\n");
      return FALSE;
    }else{
      return TRUE;
    }
  }
  return FALSE;
}

void adapt4insert(Grid *grid)
{
  Plan *plan;
  int cell, nodes[4], ranking, en[2];
  double cost;
  int report;
  double min_cost;
  int nnodeAdd, nnodeRemove;

  min_cost = 10;

  plan = planCreate( gridNCell(grid), MAX(gridNCell(grid)/10,1000) );
  for (cell=0;cell<gridMaxCell(grid);cell++) {
    if (grid==gridCell(grid, cell, nodes)) {
      cost = gridAR(grid,nodes);
      if ( cost > min_cost ) planAddItemWithPriority( plan, cell, cost );
    }
  }
  planDeriveRankingsFromPriorities( plan );
  report = 10; if (planSize(plan) > 100) report = planSize(plan)/10;
  for ( ranking=planSize(plan)-1; ranking>=0; ranking-- ) { 
    cell = planItemWithThisRanking(plan,ranking);
    if (grid==gridCell(grid, cell, nodes)) {
      cost = gridAR(grid,nodes);
      if ( ranking/report*report==ranking ){
	printf("rank %d cost %f add %d\n",ranking,cost,nnodeAdd);
	fflush(stdout);
      }
      if ( cost > min_cost ) {
	en[0]=nodes[0];	en[1]=nodes[1];
	if ( splitit(grid, en)){nnodeAdd++; continue;}
	en[0]=nodes[0];	en[1]=nodes[2];
	if ( splitit(grid, en)){nnodeAdd++; continue;}
	en[0]=nodes[0];	en[1]=nodes[3];
	if ( splitit(grid, en)){nnodeAdd++; continue;}
	en[0]=nodes[1];	en[1]=nodes[2];
	if ( splitit(grid, en)){nnodeAdd++; continue;}
	en[0]=nodes[1];	en[1]=nodes[3];
	if ( splitit(grid, en)){nnodeAdd++; continue;}
	en[0]=nodes[2];	en[1]=nodes[3];
	if ( splitit(grid, en)){nnodeAdd++; continue;}
      }
    }
  }
  planFree(plan);
}

void adapt3(Grid *grid)
{
  adapt4insert(grid);
  STATUS;
  DUMP_TEC;
  adapt3swap(grid);
  STATUS;
  DUMP_TEC;
  adapt3smooth(grid);
  STATUS;
  DUMP_TEC;
}

void relax_grid(Grid *grid)
{
  int node;
  if (TRUE) {
    STATUS;
    printf("adapt2\n");
    gridAdapt2( grid );
    STATUS;
  }
  if (TRUE) {
    STATUS;
    printf("swap\n");
    gridSwap(grid,1.0);
    STATUS;
  }
  if (TRUE) {
    STATUS;
    printf("smooth\n");
    for (node=0;node<gridMaxNode(grid);node++) {
      if ( gridValidNode(grid,node) && !gridNodeFrozen( grid, node ) ) {
	if ( gridGeometryNode( grid, node ) ) continue;
	if ( gridGeometryEdge( grid, node ) ) {
	  //gridLineSearchTForCost(grid, node );
	  continue;
	}
	if ( gridGeometryBetweenFace( grid, node ) ) continue;
	if ( gridGeometryFace( grid, node ) ) {
	  gridSmoothNodeARFace(grid, node );
	  continue;
	}
	gridSmartLaplacian(grid, node );
	gridSmoothNodeARSimplex(grid, node );
      }
    }
    STATUS;
  }
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

  int function_id=1;
  int order=1;

  Interp *temp_interp;

  int i;

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

  temp_interp = interpCreate( grid, function_id, order, order );
  interpTecplot( temp_interp, "p_fit.t" );

  gridInterp(grid) = interpContinuousReconstruction( temp_interp, 
						     order+1, order );
  interpTecplot( gridInterp(grid), "p_rec.t" );

  return 0;

  gridSetCostFunction(grid, gridCOST_FCN_INTERPOLATION );
  gridSetCostConstraint(grid, gridCOST_CNST_VOLUME );
  gridSetMinInsertCost(grid, 1.0e99 );

  h0 = 0.1;
  interp_metric(grid);
  gridCacheCurrentGridAndMap(grid);
  STATUS;

  gridHistogram(grid,"hist0.m");

  DUMP_TEC;

  for(i=0;i<1;i++) {
    adapt3(grid);
  }

  gridHistogram(grid,"hist1.m");

  gridExportFAST( grid, "grid_h1000.fgrid" );

  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: modifed grid does not have right handed boundaries\n");
  
  printf("writing output ref %s\n",ref_output);
  gridExportRef( grid, ref_output );

  printf("Done.\n");
  
  return 0;
}

