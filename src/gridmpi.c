
/* Computes metrics from faces and tets 
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "plan.h"
#include "gridmetric.h"
#include "gridinsert.h"
#include "gridswap.h"
#include "gridcad.h"
#include "gridmpi.h"

Grid *gridIdentityNodeGlobal(Grid *grid, int offset )
{
  int node;
  for (node = 0; node < gridNNode(grid) ; node++ )
    if (grid != gridSetNodeGlobal(grid,node,node+offset)) return NULL;
  return grid;
}

Grid *gridIdentityCellGlobal(Grid *grid, int offset )
{
  int cell;
  for (cell = 0; cell < gridNCell(grid) ; cell++ )
    if (grid != gridSetCellGlobal(grid,cell,cell+offset)) return NULL;
  return grid;
}

Grid *gridSetAllLocal(Grid *grid )
{
  int node;
  for (node = 0; node < gridNNode(grid) ; node++ )
    if (grid != gridSetNodePart(grid,node,gridPartId(grid))) return NULL;
  return grid;
}

Grid *gridSetGhost(Grid *grid, int node )
{
  return gridSetNodePart(grid,node,EMPTY);
}

Grid *gridParallelAdapt(Grid *grid, Queue *queue,
			double minLength, double maxLength )
{
  int n0, n1, adaptnode, origNNode, newnode;
  int nnodeAdd, nnodeRemove;
  double ratio;
  double ar, arLimit;

  GridBool use_a_plan;
  int conn, ranking;
  int nodes[2];
  double length;
  Plan *plan;
  char *env;
  
  origNNode   = gridNNode(grid);
  adaptnode   = 0;
  nnodeAdd    = 0;
  nnodeRemove = 0;

  use_a_plan = FALSE;
  env = getenv("REFINE_WITH_PLAN");
  if (env != NULL) {
    if (strcasecmp(env,"NO") == 0) use_a_plan = FALSE;
    else                           use_a_plan = TRUE;
  }

  if (use_a_plan) {

    /* split long edges */
    gridCreateConn(grid);
    plan = planCreate( gridNConn(grid)/2, MAX(gridNConn(grid)/10,1000) );
    for(conn=0;conn<gridNConn(grid);conn++) {
      gridConn2Node(grid,conn,nodes);
      ratio = gridEdgeRatio(grid,nodes[0],nodes[1]);
      if ( ratio >= maxLength ) planAddItemWithPriority( plan, conn, ratio );
    }
    planDeriveRankingsFromPriorities( plan );

    for ( ranking=planSize(plan)-1; ranking>=0; ranking-- ) { 
      conn = planItemWithThisRanking(plan,ranking);

      if (grid == gridConn2Node(grid,conn,nodes)){
	if ( gridCellEdge(grid, nodes[0], nodes[1]) &&
	     gridValidNode(grid, nodes[0]) && 
	     gridValidNode(grid, nodes[1]) && 
	     !gridNodeFrozen(grid, nodes[0]) &&
	     !gridNodeFrozen(grid, nodes[1]) &&
	     ( gridNodeLocal(grid,nodes[0]) || 
	       gridNodeLocal(grid,nodes[1]) ) ) {
	  length = gridEdgeRatio(grid, nodes[0], nodes[1]);
	  if (length >= maxLength) {
	    ratio = 0.5;
	    newnode = gridParallelEdgeSplit( grid, queue, nodes[0], nodes[1] );
	    if ( newnode != EMPTY ){
	      //gridSwapNearNode( grid, newnode, 1.0 );
	    }
	  }
	}
      }
    }
    planFree(plan);
    gridEraseConn(grid);

    /* collapse sort edges */
    gridCreateConn(grid);
    plan = planCreate( gridNConn(grid)/2, MAX(gridNConn(grid)/10,1000) );
    for(conn=0;conn<gridNConn(grid);conn++) {
      gridConn2Node(grid,conn,nodes);
      ratio = gridEdgeRatio(grid,nodes[0],nodes[1]);
      if ( ratio <= minLength ) planAddItemWithPriority( plan, conn, ratio );
    }
    planDeriveRankingsFromPriorities( plan );
  
    for ( ranking=0; ranking<planSize(plan); ranking++ ) { 
      conn = planItemWithThisRanking(plan,ranking);
      if (grid == gridConn2Node(grid,conn,nodes)){
	if ( gridCellEdge(grid, nodes[0], nodes[1]) &&
	     gridValidNode(grid, nodes[0]) && 
	     gridValidNode(grid, nodes[1]) && 
	     !gridNodeFrozen(grid, nodes[0]) &&
	     !gridNodeFrozen(grid, nodes[1]) &&
	     gridNodeLocal(grid,nodes[0]) &&
	     gridNodeLocal(grid,nodes[1]) ) {
	  length = gridEdgeRatio(grid, nodes[0], nodes[1] );
	  if (length <= minLength) {
	    if ( grid == 
		 gridParallelEdgeCollapse(grid, queue, nodes[0], nodes[1]) ) {
	      // gridSwapNearNode( grid, nodes[0], 1.0 );
	    }
	  }
	}
      }
    }
    planFree(plan);
    gridEraseConn(grid);

  }else{
    arLimit = 0.05;

    for ( n0=0; adaptnode<origNNode && n0<gridMaxNode(grid); n0++ ) { 
      adaptnode++;
      if ( gridValidNode( grid, n0) && 
	   !gridNodeFrozen( grid, n0 ) && 
	   gridNodeLocal( grid, n0 ) ) {
	if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) )
	  return NULL;
	if ( !gridNodeFrozen( grid, n1 ) && ratio > maxLength ) {
	  gridNodeAR(grid,n0,&ar);
	  if (ar > arLimit) {
	    newnode = gridParallelEdgeSplit(grid, queue, n0, n1);
	    if ( newnode != EMPTY ){
	      nnodeAdd++;
	    }
	  }
	}else{
	  if ( NULL == gridSmallestRatioEdge( grid, n0, &n1, &ratio) ) 
	    return NULL;
	  if ( !gridNodeFrozen( grid, n1 ) && ratio < minLength ) { 
	    if ( grid == gridParallelEdgeCollapse(grid, queue, n0, n1) ) {
	      nnodeRemove++;
	    }
	  }
	}
      }else{
	adaptnode++;
      }
    }
  }

#ifdef PARALLEL_VERBOSE 
  if ( NULL == queue ) {
    printf("local added%9d remov%9d AR%14.10f\n",
	   nnodeAdd,nnodeRemove,gridMinAR(grid));
  } else {
    printf("ghost added%9d remov%9d AR%14.10f\n",
	   nnodeAdd,nnodeRemove,gridMinAR(grid));
  }
#endif
  return grid;
}

Grid *gridParallelPreProject(Grid *grid, Queue *queue )
{
  int cell, conn, newnode;
  int nodes[4];
  int parent;
  double min_insert_cost;

  /* remove all tets with two or three faces on the boundary */
  /* these two- or three- face tets cause projection problems */
  for (cell=0;cell<gridMaxCell(grid);cell++){
    if ( grid==gridCell( grid, cell, nodes) ) {
      if ( NULL == queue && gridCellHasGhostNode(grid,nodes) ) continue;
      if ( gridGeometryFace(grid, nodes[0]) ||
	   gridGeometryFace(grid, nodes[1]) ||
	   gridGeometryFace(grid, nodes[2]) ||
	   gridGeometryFace(grid, nodes[3]) ){
	gridRemoveTwoFaceCell(grid, queue, cell);
	gridRemoveThreeFaceCell(grid, queue, cell);
      }
    }
  }

  /* ensure mesh is topologically correct for projection... */
  /*   by spliting any interior edge that has both nodes on boundary. */
  gridCreateConn(grid);
  min_insert_cost = gridMinInsertCost( grid ); /* save orig cost */
  gridSetMinInsertCost( grid, -100.0 ); /* split at any cost (pun intended)*/
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    if ( gridNodeLocal(grid,nodes[0]) || gridNodeLocal(grid,nodes[1]) ) {
      parent = gridParentGeometry(grid, nodes[0], nodes[1] );
      if ( ( gridGeometryFace( grid, nodes[0] ) &&
	     gridGeometryFace( grid, nodes[1] ) &&
	     0 == parent  ) ||
	   ( gridGeometryEdge( grid, nodes[0] ) &&
	     gridGeometryEdge( grid, nodes[1] ) &&
	     0 < parent  ) ) {
	newnode = gridParallelEdgeSplit( grid, queue, nodes[0], nodes[1] );
	/* to allow second pass at interior */
	if ( EMPTY == newnode && NULL != queue ){
	  newnode = gridParallelEdgeSplit( grid, NULL, nodes[0], nodes[1] );
	}
      }
    }
  }
  gridSetMinInsertCost( grid, min_insert_cost );  /* reset to orig cost */
  gridEraseConn(grid);

  return grid;
}

int gridParallelEdgeSplit(Grid *grid, Queue *queue, int node0, int node1 )
{
  GridBool gemLocal;

  if ( gridNodeGhost(grid,node0) && gridNodeGhost(grid,node1) ) return EMPTY;

  gridMakeGem(grid, node0, node1 );
  gemLocal = gridGemIsAllLocal(grid);
  if ( NULL == queue && !gemLocal) return EMPTY;
  if ( NULL != queue && gemLocal) return EMPTY;

  if (NULL != queue) queueNewTransaction(queue);
  return gridSplitEdgeRatio( grid, queue, node0, node1, 0.5 );
}

Grid *gridParallelEdgeCollapse(Grid *grid, Queue *queue, int node0, int node1 )
{
  Grid *result;

  if ( gridNodeGhost(grid,node0) && gridNodeGhost(grid,node1) ) return NULL;

  if ( NULL != queue ) return NULL;
  if ( gridNodeNearGhost(grid,node0) ) return NULL;
  if ( gridNodeNearGhost(grid,node1) ) return NULL;
  if ( !gridGeometryEdge(grid,node0) && gridGeometryBetweenFace(grid,node0) ) 
    return NULL;
  if ( !gridGeometryEdge(grid,node1) && gridGeometryBetweenFace(grid,node1) ) 
    return NULL;

  result = gridCollapseEdge(grid, queue, node0, node1, 0.5 );

  if (grid==result) {
  }

  return result;
}

Grid *gridParallelSmooth( Grid *grid, GridBool localOnly,
			  double optimizationLimit, double laplacianLimit,
                          GridBool smoothOnSurface )
{
  int node;
  double ar;
  GridBool nearGhost;
  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) && 
	 !gridNodeFrozen( grid, node ) && 
	 gridNodeLocal(grid,node) ) {
      nearGhost = gridNodeNearGhost(grid, node);
      if ( localOnly != nearGhost ) {
	gridNodeAR(grid,node,&ar);
	if (ar < optimizationLimit) {
	  gridSmoothNode( grid, node, smoothOnSurface );
	}else{
	  if (ar < laplacianLimit && !gridGeometryFace( grid, node )) {
	    gridSmartLaplacian( grid, node ); 
	  }
	}
      }
    }
  }
  return grid;
}

#define USE_LINEAR_PROGRAMMING_FOR_UNTANGLING (FALSE)

Grid *gridParallelRelaxNegativeCells( Grid *grid, 
				      GridBool localOnly,
				      GridBool smoothOnSurface )
{
  int node;
  double nodeCostValid;
  GridBool nearGhost;
  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) && 
	 !gridNodeFrozen( grid, node ) && 
	 gridNodeLocal(grid,node) ) {
      nearGhost = gridNodeNearGhost(grid, node);
      if ( localOnly != nearGhost ) {
	gridNodeCostValid(grid,node,&nodeCostValid);
	if ( -0.5 > nodeCostValid ) {
	  if (USE_LINEAR_PROGRAMMING_FOR_UNTANGLING) {
	    gridUntangleVolume( grid, node, 3, !localOnly );
	  } else {
	    gridSmoothVolumeNearNode( grid, node, 
				      smoothOnSurface );
	  }
	}
      }
    }
  }
  return grid;
}

Grid *gridParallelRelaxNegativeFaceAreaUV( Grid *grid, 
					   GridBool localOnly )
{
  int node;
  double nodeCostValid;
  GridBool nearGhost;
  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) &&
	 gridGeometryFace( grid, node ) &&
	 !gridNodeFrozen( grid, node ) && 
	 gridNodeLocal(grid,node) ) {
      nearGhost = gridNodeNearGhost(grid, node);
      if ( localOnly != nearGhost ) {
	gridNodeCostValid(grid,node,&nodeCostValid);
	if ( -1.5 > nodeCostValid ) {
	  if (USE_LINEAR_PROGRAMMING_FOR_UNTANGLING) {
	    gridUntangleAreaUV( grid, node, 1, !localOnly );
	  } else {
	    gridSmoothNodeFaceAreaUV( grid, node );
	  }
	}
      }
    }
  }
  return grid;
}

Grid *gridParallelSwap(Grid *grid, Queue *queue, double ARlimit )
{
  int cell, maxcell;
  int nodes[4];

  maxcell = gridMaxCell(grid);

  for (cell=0;cell<maxcell;cell++){
    if ( grid==gridCell( grid, cell, nodes) ) {
      if ( NULL == queue && gridCellHasGhostNode(grid,nodes) ) continue;
      if ( gridAR(grid, nodes)<ARlimit ) {
	if ( NULL == queue ) {
	  if ( grid == gridSwapFace(grid, queue, nodes[1],nodes[2],nodes[3]) )
	    continue;
	  if ( grid == gridSwapFace(grid, queue, nodes[0],nodes[2],nodes[3]) )
	    continue;
	  if ( grid == gridSwapFace(grid, queue, nodes[0],nodes[1],nodes[3]) )
	    continue;
	  if ( grid == gridSwapFace(grid, queue, nodes[0],nodes[1],nodes[2]) )
	    continue;
	}
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[0], nodes[1] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[0], nodes[2] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[0], nodes[3] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[1], nodes[2] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[1], nodes[3] ) )
	  continue;
	if ( grid == gridParallelEdgeSwap(grid, queue, nodes[2], nodes[3] ) )
	  continue;
      }
    }
  }

#ifdef PARALLEL_VERBOSE 
  printf(" final AR%14.10f\n",gridMinAR(grid));
#endif
  return grid;
}

Grid *gridParallelEdgeSwap(Grid *grid, Queue *queue, int node0, int node1 )
{
  GridBool gemLocal;
  Grid *result;

  if ( gridNodeGhost(grid,node0) && gridNodeGhost(grid,node1) ) return NULL;

  if ( grid != gridMakeGem(grid, node0, node1 ) ) {
    printf("gridParallelEdgeSwap gem error\n");
    return NULL;
  }
  gemLocal = gridGemIsAllLocal(grid);
  if ( NULL == queue && !gemLocal) return NULL;
  if ( NULL != queue && gemLocal) return NULL;

  queueNewTransaction(queue);
  result = gridSwapEdge( grid, queue, node0, node1 );

  if (NULL != result) {
    if (0==gridCellDegree(grid,node0)) gridRemoveNodeWithOutGlobal(grid,node0);
    if (0==gridCellDegree(grid,node1)) gridRemoveNodeWithOutGlobal(grid,node1);
  }

  return result;
}

Grid *gridApplyQueue(Grid *grid, Queue *gq )
{
  int transaction;
  int removed, removedcell, removedface, removededge;
  int i, globalnodes[4], globalCellId, nodeParts[4], localnodes[4];
  int cell, face, edge, faceId, edgeId;
  int added, addedcell, addedface, addededge;
  double xyz[1000], uv[6], ts[2];
  int dim, aux;
  Queue *lq;

  dim = 3 + 6 + gridNAux(grid);
  if (dim>250) printf( "ERROR: %s: %d: undersized static xyz.\n", 
		       __FILE__, __LINE__);

  lq = queueCreate( 1 ); /* only used for queuing local removed nodes */

  removedcell = 0;
  removedface = 0;
  removededge = 0;
  addedcell = 0;
  addedface = 0;
  addededge = 0;
  for (transaction=0;transaction<queueTransactions(gq);transaction++){
    for (removed=0;removed<queueRemovedCells(gq,transaction);removed++) {
      queueRemovedCellNodeParts( gq, removedcell, nodeParts );      
      if ( gridPartId(grid) == nodeParts[0] ||
	   gridPartId(grid) == nodeParts[1] ||
	   gridPartId(grid) == nodeParts[2] ||
	   gridPartId(grid) == nodeParts[3] ) {
	queueRemovedCellNodes( gq, removedcell, globalnodes );
	for(i=0;i<4;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
	cell = gridFindCell(grid,localnodes);
	if (grid == gridRemoveCellWithOutGlobal(grid,cell)) 
	  queueRemoveCell(lq,localnodes,nodeParts);
      }
      removedcell++;
    }
    for (removed=0;removed<queueRemovedFaces(gq,transaction);removed++) {
      queueRemovedFaceNodeParts( gq, removedface, nodeParts );
      if ( gridPartId(grid) == nodeParts[0] ||
	   gridPartId(grid) == nodeParts[1] ||
	   gridPartId(grid) == nodeParts[2] ) {
	queueRemovedFaceNodes( gq, removedface, globalnodes );
	for(i=0;i<3;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
	face = gridFindFace(grid,localnodes[0],localnodes[1],localnodes[2]);
	gridRemoveFace(grid,face);
      }
      removedface++;
    }
    for (removed=0;removed<queueRemovedEdges(gq,transaction);removed++) {
      queueRemovedEdgeNodeParts( gq, removededge, nodeParts );
      if ( gridPartId(grid) == nodeParts[0] ||
	   gridPartId(grid) == nodeParts[1] ) {
	queueRemovedEdgeNodes( gq, removededge, globalnodes );
	for(i=0;i<2;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
	edge = gridFindEdge(grid,localnodes[0],localnodes[1]);
	gridRemoveEdge(grid,edge);
      }
      removededge++;
    }

    for(added=0;added<queueAddedCells(gq,transaction);added++) {
      queueAddedCellNodeParts( gq, addedcell, nodeParts );
      if ( gridPartId(grid) == nodeParts[0] ||
	   gridPartId(grid) == nodeParts[1] ||
	   gridPartId(grid) == nodeParts[2] ||
	   gridPartId(grid) == nodeParts[3] ) {
	queueAddedCellNodes( gq, addedcell, globalnodes );
	queueAddedCellId( gq, addedcell, &globalCellId );
	queueAddedCellXYZs( gq, addedcell, xyz );
	for(i=0;i<4;i++) {
	  localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
	  if ( EMPTY == localnodes[i] ) {
	    localnodes[i]=gridAddNodeWithGlobal( grid,
						 xyz[0+dim*i],
						 xyz[1+dim*i],
						 xyz[2+dim*i],
						 globalnodes[i]);
	    gridSetMap(grid,localnodes[i],
		       xyz[3+dim*i],xyz[4+dim*i],xyz[5+dim*i],
		       xyz[6+dim*i],xyz[7+dim*i],xyz[8+dim*i]);
	    for ( aux = 0 ; aux < gridNAux(grid) ; aux++ )     
	      gridSetAux(grid, localnodes[i], aux, xyz[aux+9+dim*i]);
	    gridSetNodePart(grid, localnodes[i], nodeParts[i]);
	  }
	}
	cell = gridAddCellWithGlobal(grid,
				     localnodes[0],localnodes[1],
				     localnodes[2],localnodes[3],
				     globalCellId);
      }
      addedcell++;
    }
    for(added=0;added<queueAddedFaces(gq,transaction);added++) {
      queueAddedFaceNodeParts( gq, addedface, nodeParts );
      if ( gridPartId(grid) == nodeParts[0] ||
	   gridPartId(grid) == nodeParts[1] ||
	   gridPartId(grid) == nodeParts[2] ) {
	queueAddedFaceNodes( gq, addedface, globalnodes );
	queueAddedFaceId( gq, addedface, &faceId );
	queueAddedFaceUVs( gq, addedface, uv );
	for(i=0;i<3;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
	face = gridAddFaceUV(grid,
			     localnodes[0],uv[0],uv[1],
			     localnodes[1],uv[2],uv[3],
			     localnodes[2],uv[4],uv[5],
			     faceId);
      }
      addedface++;
    }
    for(added=0;added<queueAddedEdges(gq,transaction);added++) {
      queueAddedEdgeNodeParts( gq, addededge, nodeParts );
      if ( gridPartId(grid) == nodeParts[0] ||
	   gridPartId(grid) == nodeParts[1] ) {
	queueAddedEdgeNodes( gq, addededge, globalnodes );
	queueAddedEdgeId( gq, addededge, &edgeId );
	queueAddedEdgeTs( gq, addededge, ts );
	for(i=0;i<2;i++)localnodes[i]=gridGlobal2Local(grid,globalnodes[i]);
	edge = gridAddEdge(grid,localnodes[0],localnodes[1],edgeId,ts[0],ts[1]);
      }
      addededge++;
    }
    queueNewTransaction(lq);
  }
  for (removed=0;removed<queueTotalRemovedCells(lq);removed++) {
    queueRemovedCellNodes( lq, removed, localnodes );
    for (i=0;i<4;i++) {
      if ( 0 == gridCellDegree(grid,localnodes[i]) ) {
	gridRemoveNodeWithOutGlobal(grid,localnodes[i]);
      }
    }
  }
 
#ifdef PARALLEL_VERBOSE 
  printf( " %6d queue applied                           nnode %8d AR%14.10f\n",
	  gridPartId(grid),gridNNode(grid),gridMinAR(grid) );
  fflush(stdout);
#endif

  queueFree(lq);
  return grid;
}

Grid *gridGhostDataCountByPartition(Grid *grid, int total_number_of_partitions, 
				    int *partition_data_count)
{
  int node, part, faces, edges;

  for(node=0;node<total_number_of_partitions;node++) 
    partition_data_count[node] = 0;

  for(node=0;node<gridMaxNode(grid);node++) {
    if (gridNodeGhost(grid,node)) { 
      part = gridNodePart(grid,node);
      if (part < 0 || part >= total_number_of_partitions) 
	printf("%s: %d: gridNodePart error, %d part, %d npart\n",
	       __FILE__, __LINE__, part, total_number_of_partitions );
      partition_data_count[part]++;
      faces = gridNodeFaceIdDegree(grid,node);
      if (faces>0) partition_data_count[part] += (faces+1);
      edges = gridNodeEdgeIdDegree(grid,node);
      if (edges>0) partition_data_count[part] += (edges+1);
      if (faces==0 && edges>0) 
	printf("%s: %d: gridGhostDataCountByPartition error, %d faces, %d edges\n",
	       __FILE__, __LINE__, faces, edges);
    }
  }
  return grid;
}
