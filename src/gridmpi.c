
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

  arLimit = 0.05;

  origNNode   = gridNNode(grid);
  adaptnode   = 0;
  nnodeAdd    = 0;
  nnodeRemove = 0;

  for ( n0=0; adaptnode<origNNode && n0<gridMaxNode(grid); n0++ ) { 
    adaptnode++;
    if ( gridValidNode( grid, n0) && 
	 !gridNodeFrozen( grid, n0 ) && 
	 gridNodeLocal( grid, n0 ) ) {
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) return NULL;
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
	if ( -0.5 > nodeCostValid) gridSmoothVolumeNearNode( grid, node, 
							     smoothOnSurface );
      }
    }
  }
  return grid;
}

Grid *gridParallelRelaxNegativeFaceAreaUV( Grid *grid, 
					   GridBool localOnly )
{
  int node;
  double nodeFaceAreaUV;
  GridBool nearGhost;
  for (node=0;node<gridMaxNode(grid);node++) {
    if ( gridValidNode( grid, node ) &&
	 gridGeometryFace( grid, node ) &&
	 !gridNodeFrozen( grid, node ) && 
	 gridNodeLocal(grid,node) ) {
      nearGhost = gridNodeNearGhost(grid, node);
      if ( localOnly != nearGhost ) {
	gridMinFaceAreaUV(grid,node,&nodeFaceAreaUV);
	if (1.0e-14 > nodeFaceAreaUV) gridSmoothNodeFaceAreaUV( grid, node );
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
      if ( grid == gridRemoveTwoFaceCell(grid, queue, cell) ) continue;
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
