/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <values.h>

#include "queue.h"
#include "plan.h"
#include "gridmath.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridswap.h"
#include "gridinsert.h"
#include "gridfacer.h"

GridFacer *gridfacerCreate( Grid *grid, int faceId )
{
  int face;
  int nodes[3];
  int id;

  GridFacer *gf;
  gf = malloc(sizeof(GridFacer));
  gf->grid = grid;

  gf->faceId = faceId;

  gf->camera = FALSE;
  gf->tecplotFile = NULL;

  gridAttachPacker( grid, gridfacerPack, (void *)gf );
  gridAttachNodeSorter( grid, gridfacerSortNode, (void *)gf );
  gridAttachReallocator( grid, gridfacerReallocator, (void *)gf );
  gridAttachFreeNotifier( grid, gridfacerGridHasBeenFreed, (void *)gf );

  gf->nedge = 0;
  gf->maxedge = 0;
  gf->e2n = NULL;

  for ( face = 0 ; face < gridMaxFace(grid) ; face++ ) {
    if ( grid == gridFace( grid, face, nodes, &faceId ) ) {
      if ( gridfacerFaceId(gf) == faceId ) {
	gridfacerAddUniqueEdge( gf, nodes[0], nodes[1] );
	gridfacerAddUniqueEdge( gf, nodes[1], nodes[2] );
	gridfacerAddUniqueEdge( gf, nodes[2], nodes[0] );
      }
    }
  }

  return gf;
}

Grid *gridfacerGrid(GridFacer *gf)
{
  return gf->grid;
}

void gridfacerFree(GridFacer *gf)
{
  if (NULL != gf->tecplotFile) fclose(gf->tecplotFile);
  if (NULL != gf->e2n) free(gf->e2n);
  if (NULL != gf->grid) { 
    gridDetachPacker( gf->grid );
    gridDetachNodeSorter( gf->grid );
    gridDetachReallocator( gf->grid );
    gridDetachFreeNotifier( gf->grid );
  }
  free(gf);
}

void gridfacerPack(void *voidGridFacer, 
		   int nnode, int maxnode, int *nodeo2n,
		   int ncell, int maxcell, int *cello2n,
		   int nface, int maxface, int *faceo2n,
		   int nedge, int maxedge, int *edgeo2n)
{
  int i;
  GridFacer *gf = (GridFacer *)voidGridFacer;
  for ( i = 0 ; i < 2*gridfacerEdges(gf) ; i++ ) {
    gf->e2n[i] = nodeo2n[gf->e2n[i]];
  }
}

void gridfacerSortNode(void *voidGridFacer, int maxnode, int *o2n)
{
  int i;
  GridFacer *gf = (GridFacer *)voidGridFacer;
  for ( i = 0 ; i < 2*gridfacerEdges(gf) ; i++ ) {
    gf->e2n[i] = o2n[gf->e2n[i]];
  }
}

void gridfacerReallocator(void *voidGridFacer, int reallocType, 
			  int lastSize, int newSize)
{

}

void gridfacerGridHasBeenFreed(void *voidGridFacer )
{
  GridFacer *gf = (GridFacer *)voidGridFacer;
  gf->grid = NULL;
}

int gridfacerFaceId(GridFacer *gf)
{
  return gf->faceId;
}

int gridfacerEdges(GridFacer *gf)
{
  return gf->nedge;
}

GridFacer *gridfacerAddUniqueEdge(GridFacer *gf, int nodeA, int nodeB)
{
  int node0, node1;
  int edge;
  int chunk;

  node0 = MIN(nodeA,nodeB);
  node1 = MAX(nodeA,nodeB);

  for ( edge = 0 ; edge < gridfacerEdges(gf) ; edge++ ) {
    if ( (node0 == gf->e2n[0+2*edge]) && 
	 (node1 == gf->e2n[1+2*edge]) ) return gf; 
  }

  gf->nedge++;
  if ( gf->nedge >= gf->maxedge ) {
    chunk = 2000;
    gf->maxedge += chunk;
    if ( NULL == gf->e2n ) {
      gf->e2n = (int *)malloc(2 * gf->maxedge * sizeof(int));
    } else {
      gf->e2n = (int *)realloc( gf->e2n, 2 * gf->maxedge * sizeof(int));
    }
  }
  edge = gf->nedge - 1;

  gf->e2n[0+2*edge]=node0; 
  gf->e2n[1+2*edge]=node1;

  return gf;
}

GridFacer *gridfacerRemoveEdge(GridFacer *gf, int nodeA, int nodeB)
{
  int node0, node1;
  int edge;
  int found;

  if (gridfacerCameraActive(gf)) gridfacerTecplot(gf,NULL);

  node0 = MIN(nodeA,nodeB);
  node1 = MAX(nodeA,nodeB);

  found = EMPTY;
  for ( edge = 0 ; edge < gridfacerEdges(gf) ; edge++ ) {
    if ( (node0 == gf->e2n[0+2*edge]) && 
	 (node1 == gf->e2n[1+2*edge]) ) {
      found = edge;
      break;
    }
  }

  if ( EMPTY == found ) return NULL;

  gf->nedge--;
  for ( edge = found ; edge < gridfacerEdges(gf) ; edge++ ) {
    gf->e2n[0+2*edge] = gf->e2n[0+2*(edge+1)];
    gf->e2n[1+2*edge] = gf->e2n[1+2*(edge+1)];
  }

  return gf;
}

GridFacer *gridfacerExamine(GridFacer *gf)
{
  int edge;
  int node0, node1;
  double ratio;
  Grid *grid = gridfacerGrid( gf );

  for ( edge = 0 ; edge < gridfacerEdges(gf) ; edge++ ) {
    node0 = gf->e2n[0+2*edge];
    node1 = gf->e2n[1+2*edge];
    if ( 0 < gridParentGeometry(grid,node0,node1 ) ) {
      ratio = gridEdgeRatio(grid,node0,node1);
      printf("edge%8d ratio%8.4f\n",edge,ratio);
    }
  }
  return gf; 
}

GridBool gridfacerCameraActive(GridFacer *gf)
{
  return gf->camera;
}

GridFacer *gridfacerTurnCameraOn(GridFacer *gf)
{
  gf->camera = TRUE;
  return gf;
}

GridFacer *gridfacerTecplot(GridFacer *gf, char *filename)
{
  int nnode, nface;
  int *f2n, *o2n, *n2o;
  int node, face, nodes[3], faceId;
  double xyz[3], uv[2];
  Grid *grid = gridfacerGrid( gf );

  if (NULL == gf->tecplotFile) {
    if (NULL == filename) {
      char temp_filename[256];
      sprintf(temp_filename,"gridface%04d.t",gridfacerFaceId(gf));
      gf->tecplotFile = fopen(temp_filename,"w");
    }else{
      gf->tecplotFile = fopen(filename,"w");
    } 
    fprintf(gf->tecplotFile, "title=\"tecplot gridfacer\"\n");
    fprintf(gf->tecplotFile, "variables=\"X\",\"Y\",\"Z\",\"U\",\"V\"\n");
  }

  nnode=0;
  nface=0;

  f2n = (int *)malloc( 3 * gridNFace(grid) * sizeof(int) );
  for (face = 0; face < gridMaxFace(grid) ; face++ ) {
    if ( grid == gridFace(grid, face, nodes, &faceId) ) {
      if (gridfacerFaceId(gf) == faceId) {
	f2n[0+3*nface] = nodes[0];
	f2n[1+3*nface] = nodes[1];
	f2n[2+3*nface] = nodes[2];
	nface++;
      }
    }
  }

  o2n = (int *)malloc( gridMaxNode(grid) * sizeof(int) );
  for (node = 0 ; node < gridMaxNode(grid); node++ ) o2n[node] = EMPTY;
  for (face = 0; face < nface ; face++ ) {
    if ( EMPTY == o2n[f2n[0+3*face]] ) {
      o2n[f2n[0+3*face]] = nnode; nnode++;
    }
    f2n[0+3*face] = o2n[f2n[0+3*face]];
    if ( EMPTY == o2n[f2n[1+3*face]] ) {
      o2n[f2n[1+3*face]] = nnode; nnode++;
    }
    f2n[1+3*face] = o2n[f2n[1+3*face]];
    if ( EMPTY == o2n[f2n[2+3*face]] ) {
      o2n[f2n[2+3*face]] = nnode; nnode++;
    }
    f2n[2+3*face] = o2n[f2n[2+3*face]];
  }

  n2o = (int *)malloc( gridMaxNode(grid) * sizeof(int) );
  for (node = 0 ; node < gridMaxNode(grid); node++ ) {
    if (EMPTY != o2n[node]) n2o[o2n[node]] = node;
  }     
  free(o2n);
  
  fprintf(gf->tecplotFile,
	  "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	  nnode, nface);
  for (node = 0 ; node < nnode; node++ ) {
    gridNodeXYZ(grid,n2o[node],xyz);
    gridNodeUV(grid,n2o[node],gridfacerFaceId(gf),uv);
    fprintf(gf->tecplotFile, 
	    "%23.15e%23.15e%23.15e%23.15e%23.15e\n",
	    xyz[0],xyz[1],xyz[2],uv[0],uv[1]);
  }
  free(n2o);

  for ( face=0; face<nface ; face++ ){
    fprintf(gf->tecplotFile, " %9d %9d %9d\n",
	    f2n[0+3*face]+1,f2n[1+3*face]+1,f2n[2+3*face]+1);
  }
  free(f2n);

  fflush(gf->tecplotFile);

  return gf;
}

GridFacer *gridfacerRatioRange(GridFacer *gf, 
			       double *longest_ratio, double *shortest_ratio)
{
  int edge;
  int node0, node1;
  double ratio;
  Grid *grid = gridfacerGrid( gf );

  *longest_ratio  =   0.0;
  *shortest_ratio = 999.0;

  for ( edge = 0 ; edge < gridfacerEdges(gf) ; edge++ ) {
    node0 = gf->e2n[0+2*edge];
    node1 = gf->e2n[1+2*edge];
    if ( 0 < gridParentGeometry(grid,node0,node1 ) ) {
      ratio = gridEdgeRatio(grid,node0,node1);
      *longest_ratio  = MAX((*longest_ratio),ratio);
      *shortest_ratio = MIN((*shortest_ratio),ratio);
    }
  }
  return gf;
}


GridFacer *gridfacerSwap(GridFacer *gf)
{
  int edge;
  int node0, node1;
  int face0, face1;
  int nodes[3], faceId;
  int node2, node3;
  double current_ratio;
  double swapped_ratio;
  double cost, cost_improvement;

  Grid *grid = gridfacerGrid( gf );

  cost = gridMinSwapCost(grid);
  cost_improvement = gridMinSwapCostImprovement(grid);

  gridSetMinSwapCost(grid,-0.5);
  gridSetMinSwapCostImprovement(grid,-10.0);

  edge = 0;
  while ( edge < gridfacerEdges(gf) ) {
    node0 = gf->e2n[0+2*edge];
    node1 = gf->e2n[1+2*edge];
    if ( 0 > gridParentGeometry(grid,node0,node1 ) ) {
      edge++;
      continue;
    }
    face0 = gridFindFaceWithNodesUnless(grid, node0, node1, EMPTY);
    if (EMPTY == face0) {
      printf("%s: %d: face0 EMPTY.\n",__FILE__,__LINE__);
      gridSetMinSwapCost(grid,cost);
      gridSetMinSwapCostImprovement(grid,cost_improvement);
      return NULL;
    }
    face1 = gridFindFaceWithNodesUnless(grid, node0, node1, face0);
    if (EMPTY == face1) {
      printf("%s: %d: face1 EMPTY for interior face.\n",__FILE__,__LINE__);
      gridSetMinSwapCost(grid,cost);
      gridSetMinSwapCostImprovement(grid,cost_improvement);
      return NULL;
    }
    gridFace(grid,face0,nodes,&faceId);
    node2 = nodes[0] + nodes[1] + nodes[2] - node0 - node1;
    gridFace(grid,face1,nodes,&faceId);
    node3 = nodes[0] + nodes[1] + nodes[2] - node0 - node1;

    current_ratio = gridEdgeRatio(grid,node0,node1);
    swapped_ratio = gridEdgeRatio(grid,node2,node3);

    if ( ((swapped_ratio < current_ratio) && (current_ratio > 1.0)) ||
	 ((swapped_ratio > current_ratio) && (swapped_ratio < 1.0)) ){
      if ( grid == gridSwapEdge( grid, NULL, node0, node1 ) ) {
	gridfacerRemoveEdge(gf, node0, node1);
	gridfacerAddUniqueEdge(gf, node2, node3);
	continue;
      }
    }
    edge++;
  }

  gridSetMinSwapCost(grid,cost);
  gridSetMinSwapCostImprovement(grid,cost_improvement);
  return gf; 
}

GridFacer *gridfacerSplit(GridFacer *gf)
{
  int edge;
  int node0, node1;
  double ratio;
  int rank;
  int *local_e2n;
  int face0, face1;
  int nodes[3], faceId;
  int node2, node3;
  int newnode;

  Plan *plan;

  Grid *grid = gridfacerGrid( gf );

  plan = planCreate( gridfacerEdges(gf), 100 );

  for ( edge = 0 ; edge < gridfacerEdges(gf) ; edge++ ) {
    node0 = gf->e2n[0+2*edge];
    node1 = gf->e2n[1+2*edge];
    if ( 0 < gridParentGeometry(grid,node0,node1 ) ) {
      ratio = gridEdgeRatio(grid,node0,node1);
      if ( ratio > 1.0 ) {
	planAddItemWithPriority( plan, edge, ratio );
      }
    }
  }

  planDeriveRankingsFromPriorities( plan );

  local_e2n = (int *)malloc( 2 * planSize( plan ) * sizeof(int) );
  for ( rank = planSize( plan )-1 ; rank >=0  ; rank-- ) {
    edge = planItemWithThisRanking( plan, rank);
    local_e2n[0+2*rank] = gf->e2n[0+2*edge];
    local_e2n[1+2*rank] = gf->e2n[1+2*edge];
  }

  for ( rank = planSize( plan )-1 ; rank >=0  ; rank-- ) {
    node0 = local_e2n[0+2*rank];
    node1 = local_e2n[1+2*rank];
    face0 = gridFindFaceWithNodesUnless(grid, node0, node1, EMPTY);
    if (EMPTY == face0) {
      printf("%s: %d: face0 EMPTY.\n",__FILE__,__LINE__);
      return NULL;
    }
    face1 = gridFindFaceWithNodesUnless(grid, node0, node1, face0);
    if (EMPTY == face1) {
      printf("%s: %d: face1 EMPTY for interior face.\n",__FILE__,__LINE__);
      return NULL;
    }
    gridFace(grid,face0,nodes,&faceId);
    node2 = nodes[0] + nodes[1] + nodes[2] - node0 - node1;
    gridFace(grid,face1,nodes,&faceId);
    node3 = nodes[0] + nodes[1] + nodes[2] - node0 - node1;
    newnode = gridSplitEdgeRatio(grid, NULL, node0, node1, 0.5);
    if ( EMPTY != newnode ) {
      gridfacerRemoveEdge(gf, node0, node1);
      gridfacerAddUniqueEdge(gf, node0, newnode);
      gridfacerAddUniqueEdge(gf, node1, newnode);
      gridfacerAddUniqueEdge(gf, node2, newnode);
      gridfacerAddUniqueEdge(gf, node3, newnode);
      if (grid!=gridUntangle(grid)) {
	printf("%s: %d: gridUntangle NULL.\n",__FILE__,__LINE__);
	return NULL;
      }
    }
  }
  free(local_e2n);
  planFree(plan);
  return gf; 
}

GridFacer *gridfacerSplitProblemProjectionEdges(GridFacer *gf) {
  int conn, nodes[2];
  int parent, newnode;
  int split_edges;
  double min_insert_cost;

  Grid *grid = gridfacerGrid( gf );

  min_insert_cost = gridMinInsertCost( grid );
  gridSetMinInsertCost( grid, -100.0 );
  gridCreateConn(grid);
  split_edges = 0;
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    parent = gridParentGeometry(grid, nodes[0], nodes[1] );
    if ( ( gridGeometryFace( grid, nodes[0] ) &&
	   gridGeometryFace( grid, nodes[1] ) &&
	   0 == parent  ) ) {
      newnode = gridSplitEdgeRatio( grid, NULL, nodes[0], nodes[1], 0.5);
      split_edges++;
    }
  }
  printf("split %d problem edges\n",split_edges);
  gridEraseConn(grid);
  gridSetMinInsertCost( grid, min_insert_cost );
  return gf;
}


GridFacer *gridfacerCollapseEdge( GridFacer *gf, int node0, int node1 )
{
  AdjIterator it;
  int nodes[3], faceId;
  double ratio0, ratio1;

  Grid *grid = gridfacerGrid( gf );

  printf("%s: %d: gridfacerCollapseEdge not implemented.\n",
	 __FILE__,__LINE__);
  return NULL;
  
  for ( it = adjFirst(gridFaceAdj(grid),node1); 
	adjValid(it); 
	it = adjNext(it) ) {
    gridFace(grid, adjItem(it), nodes, &faceId);

    ratio0 = gridEdgeRatio(grid,nodes[0],nodes[1]);
    ratio1 = ratio0;
    if ( nodes[0] == node0 ) ratio1 = gridEdgeRatio(grid,node0,nodes[1]);
    if ( nodes[1] == node0 ) ratio1 = gridEdgeRatio(grid,nodes[0],node0);

  }

  return NULL;
}

GridFacer *gridfacerCollapse(GridFacer *gf)
{
  int edge;
  int node0, node1;
  double limit;
  double ratio;
  int rank;
  int *local_e2n;
  int face0, face1;
  int nodes[3], faceId;
  int node2, node3;
  int newnode;

  Plan *plan;

  Grid *grid = gridfacerGrid( gf );

  plan = planCreate( gridfacerEdges(gf), 100 );
  limit = 0.2;
  for ( edge = 0 ; edge < gridfacerEdges(gf) ; edge++ ) {
    node0 = gf->e2n[0+2*edge];
    node1 = gf->e2n[1+2*edge];
    if ( 0 < gridParentGeometry(grid,node0,node1 ) ) {
      ratio = gridEdgeRatio(grid,node0,node1);
      if ( ratio < limit ) {
	planAddItemWithPriority( plan, edge, ratio );
      }
    }
  }

  planDeriveRankingsFromPriorities( plan );

  local_e2n = (int *)malloc( 2 * planSize( plan ) * sizeof(int) );
  for ( rank = 0 ; rank < planSize( plan ) ; rank++ ) {
    edge = planItemWithThisRanking( plan, rank);
    local_e2n[0+2*rank] = gf->e2n[0+2*edge];
    local_e2n[1+2*rank] = gf->e2n[1+2*edge];
  }

  for ( rank = 0 ; rank < planSize( plan ) ; rank++ ) {
    node0 = local_e2n[0+2*rank];
    node1 = local_e2n[1+2*rank];
    ratio = gridEdgeRatio(grid,node0,node1);
    if (ratio >= limit) continue;
    printf("rank %d len%10.6f\n",rank,ratio);
    if (gf != gridfacerCollapseEdge( gf, node0, node1) ) 
      gridfacerCollapseEdge( gf, node1, node0);
  }
  free(local_e2n);
  planFree(plan);
  return gf; 
}

