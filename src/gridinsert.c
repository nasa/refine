
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
#include <values.h>
#include "plan.h"
#include "ring.h"
#include "gridshape.h"
#include "gridmath.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "gridinsert.h"
#include "gridcad.h"

Grid *gridThrash(Grid *grid)
{
  int ncell, nodelimit, cellId, nodes[4];
  ncell = gridNCell(grid);
  nodelimit = gridNNode(grid)*3/2;
  for (cellId=0;cellId<ncell && gridNNode(grid)<nodelimit;cellId++)
    if ( NULL != gridCell( grid, cellId, nodes) )
      gridSplitEdge( grid, nodes[0], nodes[1] );
  
  return grid;
}

Grid *gridRemoveAllNodes(Grid *grid )
{
  AdjIterator it;
  int n0, i, cell, nodes[4];
  GridBool nodeExists;
 
  for ( n0=0; n0<gridMaxNode(grid); n0++ ) { 
    if ( gridValidNode(grid, n0) && !gridNodeFrozen(grid, n0) ) {
      nodeExists = TRUE;
      for ( it = adjFirst(gridCellAdj(grid),n0); 
	    nodeExists && adjValid(it); 
	    it = adjNext(it) ){
	cell = adjItem(it);
	gridCell( grid, cell, nodes);
	for (i=0;nodeExists && i<4;i++){
	  if (n0 != nodes[i]) {
	    nodeExists = (grid == gridCollapseEdge(grid,NULL,nodes[i],n0,1.0));
	  }
	}
      }
    }
  }
  return grid;
}

Grid *gridAdapt(Grid *grid, double minLength, double maxLength )
{
  int n0, n1, adaptnode, origNNode, newnode;
  int report, nnodeAdd, nnodeRemove;
  double ratio;
  double ar;
  
  origNNode = gridNNode(grid);
  adaptnode =0;
  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (gridNNode(grid) > 100) report = gridNNode(grid)/10;

  for ( n0=0; 
	adaptnode<origNNode && n0<gridMaxNode(grid); 
	n0++ ) { 
    if (adaptnode > 100 &&adaptnode/report*report == adaptnode ) {
      printf("adapt node %8d nnode %8d added %8d removed %8d\n",
	     adaptnode,gridNNode(grid),nnodeAdd,nnodeRemove);
    }
    if ( gridValidNode( grid, n0) && !gridNodeFrozen( grid, n0 ) ) {
      adaptnode++;
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) return NULL;
      if ( !gridNodeFrozen( grid, n1 ) && ratio > maxLength ) {
	newnode = gridSplitEdge(grid, n0, n1);
	if ( newnode != EMPTY ){
	  nnodeAdd++;
	  gridNodeAR(grid,newnode,&ar); 
	  if (ar<0.0) {
	    printf("%s: %d: split to node %d has ar %f\n",
		   __FILE__,__LINE__,newnode,ar);
	    exit(1);
	  }   
	  gridSwapNearNode( grid, newnode, -1.0 );
	  gridNodeAR(grid,newnode,&ar); 
	  if (ar<0.0) {
	    printf("%s: %d: split swap to node %d has ar %f\n",
		   __FILE__,__LINE__,newnode,ar);
	    exit(1);
	  }   
	}
      }else{
	if ( NULL == gridSmallestRatioEdge( grid, n0, &n1, &ratio) ) 
	  return NULL;
	if ( !gridNodeFrozen( grid, n1 ) && ratio < minLength ) { 
	  if ( grid == gridCollapseEdge(grid, NULL, n0, n1, 0.5) ) {
	    nnodeRemove++;
	    gridNodeAR(grid,n0,&ar); 
	    if (ar<0.0) {
	      printf("%s: %d: collapse to node %d has ar %f\n",
		     __FILE__,__LINE__,n0,ar);
	      exit(1);
	    }   
	    gridSwapNearNode( grid, n0, -1.0 );
	    gridNodeAR(grid,n0,&ar); 
	    if (ar<0.0) {
	      printf("%s: %d: collapse swap to node %d has ar %f\n",
		     __FILE__,__LINE__,n0,ar);
	      exit(1);
	    }   
	  }
	}
      }
    }else{
      adaptnode++;
    }
  }
  return grid;
}

static Grid *gridCollapseEdgeToBestExistingConfiguration( Grid *grid,
							  Queue *queue, 
							  int node0, int node1 )
{
  double currentCost, node0Cost, node1Cost;
  double ratio;
  if (grid!=gridCollapseCost(grid, node0, node1, 
			     &currentCost, &node0Cost, &node1Cost))
    return NULL;
  if (node0Cost<currentCost && node1Cost<currentCost) return NULL;
  ratio = (node0Cost>node1Cost?0.0:1.0);
  if (grid==gridCollapseEdge(grid, queue, node0, node1, ratio)) return grid;
  if (MIN(node0Cost,node1Cost)>currentCost)
    return gridCollapseEdge(grid, queue, node0, node1, 1.0-ratio);
  return NULL;
}

Grid *gridAdaptBasedOnConnRankings(Grid *grid )
{
  int ranking, conn, nodes[2];
  int report, nnodeAdd, nnodeRemove;
  double ratios[3];
  double dist, ratio;
  int i, newnode;
  Plan *plan;

  gridCreateConn(grid);
  plan = planCreate( gridNConn(grid), MAX(gridNConn(grid)/10,1000) );
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    planAddItemWithPriority( plan,
			     conn, gridEdgeRatioError(grid,nodes[0],nodes[1]) );
  }
  planDeriveRankingsFromPriorities( plan );

  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (planSize(plan) > 100) report = planSize(plan)/10;

  for ( ranking=planSize(plan)-1; ranking>=0; ranking-- ) { 
    conn = planItemWithThisRanking(plan,ranking);
    if (ranking/report*report == ranking || ranking==gridNConn(grid)-1) {
      printf("adapt ranking%9d nnode%9d added%9d removed%9d err%6.2f\n",
	     ranking,gridNNode(grid),nnodeAdd,nnodeRemove,
	     planPriorityWithThisRanking(plan,ranking));
    }
    if (grid == gridConn2Node(grid,conn,nodes)){
      if ( gridCellEdge(grid, nodes[0], nodes[1]) &&
	   gridValidNode(grid, nodes[0]) && 
	   gridValidNode(grid, nodes[1]) && 
	   !gridNodeFrozen(grid, nodes[0]) &&
	   !gridNodeFrozen(grid, nodes[1]) ) {
	if (grid == gridEdgeRatio3(grid, nodes[0], nodes[1], ratios ) ) {
	  if ( ratios[2] > 1.55 ) {
	    if (ratios[0]<ratios[1]) {
	      dist = 1.0/ratios[0];
              ratio = dist;
	    }else{
	      dist = 1.0/ratios[1];
              ratio = (1.0-dist);
	    }
	    newnode = gridSplitEdgeRatio( grid, NULL,
                                          nodes[0], nodes[1], ratio );
	    if ( newnode != EMPTY ){
	      nnodeAdd++;
	      gridSwapNearNode( grid, newnode, 1.0 );
	    }
	  }else if (ratios[2] < 1.0/1.35) {
	    if ( grid == 
		 gridCollapseEdgeToBestExistingConfiguration(grid, NULL, 
							      nodes[0], 
							      nodes[1] ) ) {
	      nnodeRemove++;
  	      gridSwapNearNode( grid, nodes[0], 1.0 );
	    }
	  }
	}
      }
    }
  }
  planFree(plan);
  gridEraseConn(grid);

  return grid;
}

Grid *gridCollapseEdgeToAnything( Grid *grid,
				  Queue *queue, 
				  int node0, int node1 )
{
  double currentCost, node0Cost, node1Cost;
  double ratio;
  if (grid!=gridCollapseCost(grid, node0, node1, 
			     &currentCost, &node0Cost, &node1Cost))
    return NULL;
  ratio = (node0Cost>node1Cost?0.0:1.0);
  if (grid==gridCollapseEdge(grid, queue, node0, node1, ratio)) return grid;
  ratio = 1-ratio;
  if (grid==gridCollapseEdge(grid, queue, node0, node1, ratio)) return grid;
  ratio = 0.5;
  if (grid==gridCollapseEdge(grid, queue, node0, node1, ratio)) return grid;
  return NULL;
}

Grid *gridAdaptLongShortCurved(Grid *grid, double minLength, double maxLength,
			       GridBool debug_split )
{
  int ranking, conn, nodes[2];
  int report, nnodeAdd, nnodeRemove;
  double length, ratio;
  int i, newnode;
  Plan *plan;

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
    }
    if (grid == gridConn2Node(grid,conn,nodes)){
      if ( gridCellEdge(grid, nodes[0], nodes[1]) &&
	   gridValidNode(grid, nodes[0]) && 
	   gridValidNode(grid, nodes[1]) && 
	   !gridNodeFrozen(grid, nodes[0]) &&
	   !gridNodeFrozen(grid, nodes[1]) ) {
	length = gridEdgeRatio(grid, nodes[0], nodes[1]);
	if (length >= maxLength) {
	  ratio = 0.5;
	  newnode = gridSplitEdgeRatio( grid, NULL,
					nodes[0], nodes[1], ratio );
	  if ( newnode != EMPTY ){
	    nnodeAdd++;
	    gridSwapNearNode( grid, newnode, 1.0 );
	  } else {
	    newnode = gridSplitEdgeRepeat( grid, NULL, nodes[0], nodes[1],
					   debug_split );
	    if ( newnode != EMPTY ){
	      nnodeAdd++;
	      gridSwapNearNode( grid, newnode, 1.0 );
	    }else{
	      if (debug_split) {
		printf("Edge%10d%10d will not split face%2d%2d err%6.2f\n",
		       nodes[0],nodes[1],
		       gridGeometryFace(grid,nodes[0]),
		       gridGeometryFace(grid,nodes[1]),
		       length);
		if (0 != gridParentGeometry(grid, nodes[0], nodes[1]) ) {
		  int igem;
		  gridWriteTecplotEquator(grid, nodes[0], nodes[1],
					  "edge_split_equator.t");
		  for ( igem=0 ; igem<gridNGem(grid) ; igem++ ){
		    gridWriteTecplotCellJacDet(grid,gridGem(grid,igem),NULL);
		  }
		    return NULL;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  planFree(plan);
  gridEraseConn(grid);

  gridCreateConn(grid);
  plan = planCreate( gridNConn(grid)/2, MAX(gridNConn(grid)/10,1000) );
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    length = gridEdgeRatio(grid,nodes[0],nodes[1]);
    if ( length <= minLength ) planAddItemWithPriority( plan, conn, length );
  }
  planDeriveRankingsFromPriorities( plan );
  
  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (planSize(plan) > 100) report = planSize(plan)/10;

  for ( ranking=0; ranking<planSize(plan); ranking++ ) { 
    conn = planItemWithThisRanking(plan,ranking);
    if (ranking/report*report == ranking || ranking==gridNConn(grid)-1) {
      printf("adapt ranking%9d nnode%9d added%9d removed%9d err%6.2f\n",
	     ranking,gridNNode(grid),nnodeAdd,nnodeRemove,
	     planPriorityWithThisRanking(plan,ranking));
    }
    if (grid == gridConn2Node(grid,conn,nodes)){
      if ( gridCellEdge(grid, nodes[0], nodes[1]) &&
	   gridValidNode(grid, nodes[0]) && 
	   gridValidNode(grid, nodes[1]) && 
	   !gridNodeFrozen(grid, nodes[0]) &&
	   !gridNodeFrozen(grid, nodes[1]) ) {
	length = gridEdgeRatio(grid, nodes[0], nodes[1] );
	if (length <= minLength) {
	  if ( grid == 
	       gridCollapseEdgeToAnything(grid, NULL, 
					  nodes[0], 
					  nodes[1] ) ) {
	    nnodeRemove++;
	    gridSwapNearNode( grid, nodes[0], 1.0 );
	  }
	}
      }
    }
  } 
  planFree(plan);
  gridEraseConn(grid);

  return grid;
}
Grid *gridAdaptLongShortLinear(Grid *grid, double minLength, double maxLength,
			       GridBool debug_split )
{
  int ranking, conn, nodes[2];
  int report, nnodeAdd, nnodeRemove;
  double length, ratio;
  int i, newnode;
  Plan *plan;

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
    }
    if (grid == gridConn2Node(grid,conn,nodes)){
      if ( gridCellEdge(grid, nodes[0], nodes[1]) &&
	   gridValidNode(grid, nodes[0]) && 
	   gridValidNode(grid, nodes[1]) && 
	   !gridNodeFrozen(grid, nodes[0]) &&
	   !gridNodeFrozen(grid, nodes[1]) ) {
	length = gridEdgeRatio(grid, nodes[0], nodes[1]);
	if (length >= maxLength) {
	  ratio = 0.5;
	  newnode = gridSplitEdgeRatio( grid, NULL,
					nodes[0], nodes[1], ratio );
	  if ( newnode != EMPTY ){
	    nnodeAdd++;
	    gridSwapNearNode( grid, newnode, 1.0 );
	  } else {
	    newnode = gridReconstructSplitEdgeRatio( grid, NULL,
						     nodes[0], nodes[1],
						     ratio );
	    if (newnode != EMPTY) {
	      nnodeAdd++;
	      gridSwapNearNode( grid, newnode, 1.0 );
	    }else{
	      printf("Edge%10d%10d will not split face%2d%2d err%6.2f\n",
		     nodes[0],nodes[1],
		     gridGeometryFace(grid,nodes[0]),
		     gridGeometryFace(grid,nodes[1]),
		     length);

	      gridWriteTecplotEquator(grid, nodes[0], nodes[1],
				      "edge_split_equator.t");
	      
	      return NULL;
	    }
	  }
	}
      }
    }
  }
  planFree(plan);
  gridEraseConn(grid);

  gridCreateConn(grid);
  plan = planCreate( gridNConn(grid)/2, MAX(gridNConn(grid)/10,1000) );
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    length = gridEdgeRatio(grid,nodes[0],nodes[1]);
    if ( length <= minLength ) planAddItemWithPriority( plan, conn, length );
  }
  planDeriveRankingsFromPriorities( plan );
  
  nnodeAdd = 0;
  nnodeRemove = 0;

  report = 10; if (planSize(plan) > 100) report = planSize(plan)/10;

  for ( ranking=0; ranking<planSize(plan); ranking++ ) { 
    conn = planItemWithThisRanking(plan,ranking);
    if (ranking/report*report == ranking || ranking==gridNConn(grid)-1) {
      printf("adapt ranking%9d nnode%9d added%9d removed%9d err%6.2f\n",
	     ranking,gridNNode(grid),nnodeAdd,nnodeRemove,
	     planPriorityWithThisRanking(plan,ranking));
    }
    if (grid == gridConn2Node(grid,conn,nodes)){
      if ( gridCellEdge(grid, nodes[0], nodes[1]) &&
	   gridValidNode(grid, nodes[0]) && 
	   gridValidNode(grid, nodes[1]) && 
	   !gridNodeFrozen(grid, nodes[0]) &&
	   !gridNodeFrozen(grid, nodes[1]) ) {
	length = gridEdgeRatio(grid, nodes[0], nodes[1] );
	if (length <= minLength) {
	  if ( grid == 
	       gridCollapseEdgeToAnything(grid, NULL, 
					  nodes[0], 
					  nodes[1] ) ) {
	    nnodeRemove++;
	    gridSwapNearNode( grid, nodes[0], 1.0 );
	  }
	}
      }
    }
  } 
  planFree(plan);
  gridEraseConn(grid);

  return grid;
}

int gridFindEnclosingCell(Grid *grid, int starting_guess, double *target )
{
  int current_cell;
  int tries;
  AdjIterator it;

  it = adjFirst(gridCellAdj(grid),starting_guess);
  current_cell = adjItem(it);
  if (EMPTY == current_cell) return EMPTY;
  for (tries=1;tries<=1000;tries++) {
    int nodes[4];
    double newxyz[3], xyz0[3], xyz1[3], xyz2[3], xyz3[3];
    double edge0[3], edge1[3];
    double norm0[3], norm1[3], norm2[3], norm3[3];
    double dir0[3], dir1[3], dir2[3], dir3[3];
    double dot0, dot1, dot2, dot3;
    int other_cell;
    
    gridCell( grid, current_cell, nodes );
    // gridWriteTecplotCellGeom( grid, nodes, NULL, NULL);

    gridNodeXYZ(grid, nodes[0], xyz0);
    gridNodeXYZ(grid, nodes[1], xyz1);
    gridNodeXYZ(grid, nodes[2], xyz2);
    gridNodeXYZ(grid, nodes[3], xyz3);

    gridSubtractVector(xyz3, xyz1, edge0);
    gridSubtractVector(xyz2, xyz1, edge1);
    gridCrossProduct(edge0,edge1,norm0);
    gridSubtractVector(target, xyz1, dir0);
    dot0 = gridDotProduct(dir0,norm0);

    gridSubtractVector(xyz2, xyz0, edge0);
    gridSubtractVector(xyz3, xyz0, edge1);
    gridCrossProduct(edge0,edge1,norm1);
    gridSubtractVector(target, xyz0, dir1);
    dot1 = gridDotProduct(dir1,norm1);

    gridSubtractVector(xyz3, xyz0, edge0);
    gridSubtractVector(xyz1, xyz0, edge1);
    gridCrossProduct(edge0,edge1,norm2);
    gridSubtractVector(target, xyz0, dir2);
    dot2 = gridDotProduct(dir2,norm2);

    gridSubtractVector(xyz1, xyz0, edge0);
    gridSubtractVector(xyz2, xyz0, edge1);
    gridCrossProduct(edge0,edge1,norm3);
    gridSubtractVector(target, xyz0, dir3);
    dot3 = gridDotProduct(dir3,norm3);

    // printf("cell%11d dots%23.15e%23.15e%23.15e%23.15e\n",
    //   current_cell, dot0,dot1,dot2,dot3);

    if ( dot0 >= 0.0 && dot1 >= 0.0 &&  dot2 >= 0.0 && dot3 >= 0.0 ) {
      return current_cell;
    }

    if  ( dot0 < 0.0 ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[1], nodes[2], nodes[3],
					       current_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( dot1 < 0.0 ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[2], nodes[3],
					       current_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( dot2 < 0.0 ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[1], nodes[3],
					       current_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    if  ( dot3 < 0.0 ) {
      other_cell = gridFindOtherCellWith3Nodes(grid,
					       nodes[0], nodes[1], nodes[2],
					       current_cell );
      if (EMPTY != other_cell ) {
	current_cell = other_cell;
	continue;
      }
    }

    it = adjNext(it);
    current_cell = adjItem(it);
    if (EMPTY == current_cell) {
      printf("%s: %d: gridFindEnclosingCell round-off: no more next \n");
      return EMPTY;
    }
  }

  printf("%s: %d: gridFindEnclosingCell exhausted tries.\n",
	 __FILE__,__LINE__);

  return EMPTY;
}

Grid *gridThreadCurveThroughVolume(Grid *grid, int parent, 
				   int n0, double *tuv0,
				   int n1, double *tuv1,
				   int *nextnode, double *nexttuv )
{
  AdjIterator it;
  int cell;
  int cell_nodes[4];
  int face_nodes[3];
  double xyz[3], bary[3];
  double origxyz[3];
  int node_index;

  (*nextnode) = EMPTY; nexttuv[0] = DBL_MAX; nexttuv[1] = DBL_MAX;

  if ( !gridValidNode(grid, n0) || !gridValidNode(grid, n1) ) return NULL; 

  if (gridCellEdge( grid, n0, n1)) {
    printf("gridThreadCurveThroughVolume already have %d %d.\n",n0,n1);
    (*nextnode) = n1; nexttuv[0] = tuv1[0]; nexttuv[1] = tuv1[1];
    return grid;
  }

  if (parent>0){
    printf("thread face curve %d %d.\n",n0,n1);
  }else{
    printf("thread edge curve %d %d.\n",n0,n1);
  }

  it = adjFirst(gridCellAdj(grid),n0);
  while (adjValid(it)){
    cell = adjItem(it);
    gridCell(grid, cell, cell_nodes);
    gridFaceOppositeCellNode(grid, cell_nodes, n0, face_nodes);
    
    //printf("cell %d.\n",cell);

    if ( grid == gridCurveIntersectsFace(grid, face_nodes,
					 parent, tuv0, tuv1,
					 nexttuv, xyz, bary ) ) {
      printf("intersect at tuv %f %f.\n",nexttuv[0],nexttuv[1]);
      printf("intersect at xyz %f %f %f.\n",xyz[0],xyz[1],xyz[2]);

      for (node_index=0;EMPTY==(*nextnode) && node_index<3;node_index++) {
	if (EMPTY == (*nextnode) && bary[node_index] > 0.95) {
	  double minAR;
	  gridNodeXYZ(grid, face_nodes[node_index], origxyz);
	  gridSetNodeXYZ(grid, face_nodes[node_index], xyz);
	  gridNodeAR(grid, face_nodes[node_index], &minAR );
	  if (minAR < gridMinInsertCost(grid) ) {
	    gridSetNodeXYZ(grid, face_nodes[node_index], origxyz);
	  } else {
	    (*nextnode) = face_nodes[node_index];
	    printf("move node %d of opposite face.\n",(*nextnode));
	    return grid;
	  }
	}
	if (EMPTY == (*nextnode) && bary[node_index] < 0.05) {
	  int edge0, edge1;
	  edge0 = node_index + 1; if (edge0>2) edge0 -= 3;
	  edge1 = node_index + 2; if (edge1>2) edge1 -= 3;
	  (*nextnode) = gridSplitEdgeRatio(grid, NULL,
				       face_nodes[edge0], face_nodes[edge1],
				       bary[edge1]);
	  if (EMPTY != (*nextnode)) {
	    printf("split edge %d of opposite face.\n",(*nextnode));
	    return grid;
	  }
	}
      }
      if (EMPTY == (*nextnode)) {
	printf("try to split opposite face.\n");
	(*nextnode) = gridSplitFaceAt(grid, face_nodes, xyz);
	if (EMPTY != (*nextnode)) {
	  printf("split opposite face with %d.\n",(*nextnode));
	  return grid;
	}
      }
    }
    
    it = adjNext(it);
  }

  printf("%s: %d: gridThreadCurveThroughVolume couldnt find a candidate face\n",
	 __FILE__,__LINE__);

  for ( it = adjFirst(gridCellAdj(grid),n0); adjValid(it); it = adjNext(it) ) {
    cell = adjItem(it);
    gridCell(grid, cell, cell_nodes);
    gridWriteTecplotCellGeom(grid, cell_nodes, NULL, "edge_split_equator.t");
  }

  return NULL;
}

Grid *gridRemoveCellsOutsideOfFaces( Grid *grid, int n0, int n1, int n2 )
{
  int cell0, cell1;
  int original_nodes[4], nodes[4];

  cell0 = gridFindOtherCellWith3Nodes( grid, n0, n1, n2, EMPTY );
  if (EMPTY == cell0) return grid;

  if (grid != gridCell(grid, cell0, original_nodes) ) return NULL;
  nodes[0] = n0; nodes[1] = n1;

  if (grid != gridOrient(grid, original_nodes, nodes ) ) return NULL;

  if ( nodes[2] == n2 ) {
    gridRemoveCell(grid, cell0);
  }else{
    cell1 = gridFindOtherCellWith3Nodes( grid, n0, n1, n2, cell0 );
    if (EMPTY == cell1) return grid;

    if (grid != gridCell(grid, cell1, original_nodes) ) return NULL;
    nodes[0] = n0; nodes[1] = n1;

    if (grid != gridOrient(grid, original_nodes, nodes ) ) return NULL;
    gridRemoveCell(grid, cell1);
  }

  if (EMPTY == gridFindFace(grid, nodes[1], nodes[2], nodes[3] ) ) {
    if (grid != gridRemoveCellsOutsideOfFaces(grid,
					      nodes[1],
					      nodes[2],
					      nodes[3])) return NULL;
  }
  
  if (EMPTY == gridFindFace(grid, nodes[0], nodes[3], nodes[2] ) ) {
    if (grid != gridRemoveCellsOutsideOfFaces(grid,
					      nodes[0],
					      nodes[3],
					      nodes[2])) return NULL;
  }
  
  if (EMPTY == gridFindFace(grid, nodes[0], nodes[1], nodes[3] ) ) {
    if (grid != gridRemoveCellsOutsideOfFaces(grid,
					      nodes[0],
					      nodes[1],
					      nodes[3])) return NULL;
  }
  
  if ( 0 == gridCellDegree(grid, nodes[0]) ) gridRemoveNode(grid, nodes[0]);
  if ( 0 == gridCellDegree(grid, nodes[1]) ) gridRemoveNode(grid, nodes[1]);
  if ( 0 == gridCellDegree(grid, nodes[2]) ) gridRemoveNode(grid, nodes[2]);
  if ( 0 == gridCellDegree(grid, nodes[3]) ) gridRemoveNode(grid, nodes[3]);
  
  return grid;

}

Grid *gridFillRingWithFaces( Grid *grid, Ring *ring, int faceId )
{
  int equator;
  int node0, node1, node2;
  double uv0[2], uv1[2], uv2[2];
  double area;
  GridBool could_not_complete_triangle;

  while ( ringSegments(ring) > 0 ) {
    ringInspect( ring );
    if (ring != ringSegment( ring, 0, // always grab first segment off stack 
			     &node0, &node1,
			     uv0, uv1 ) ) {
      printf("%s: %d: gridFillRingWithFaces: ringSegment NULL\n",
	     __FILE__,__LINE__);    
      return NULL;
    }

    if (grid != gridEquator(grid, node0, node1)) {
      printf("%s: %d: gridFillRingWithFaces: gridEquator NULL\n",
	     __FILE__,__LINE__);    
      return NULL;
    }
    
    could_not_complete_triangle = TRUE;
    for( equator = 0 ; 
	 could_not_complete_triangle && equator < gridNEqu(grid) ;
	 equator++ ) {

      node2 = gridEqu( grid, equator );
      printf("close nodes %6d%6d%6d\n",
	     node0,node1,node2);
      if ( ringSegmentsContainNode( ring, node2, uv2) ) {
	area = gridFaceAreaUVDirect(grid,uv0,uv1,uv2,faceId);
	printf("area %e\n",area);
	if (area>1.0e-12 ) {
	  if (ring == ringAddTriangle( ring, node0, node1, node2, uv2 ) ) {
	    could_not_complete_triangle = FALSE;
	  }else{
	    printf("%s: %d: gridFillRingWithFaces: ringAddTriangle NULL\n",
		   __FILE__,__LINE__);    
	    return NULL;
	  }
	}
      }
    }
    if (could_not_complete_triangle){
      printf("%s: %d: gridFillRingWithFaces: could_not_complete_triangle\n",
	     __FILE__,__LINE__);    
      return NULL;
    }
  }
  printf("remaining segments %d\n",ringSegments(ring));
  return grid;
}

int gridReconstructSplitEdgeRatio(Grid *grid, Queue *queue,
				  int n0, int n1, double ratio )
{
  int gap0, gap1, face0, face1, nodes[3], faceId0, faceId1;
  double uvgap0[2], uvgap1[2];
  int parent;
  double xyz[3];
  double tuv0[2], tuv1[2], newtuv[2];
  int newnode;
  int enclosing_cell;
  int node, nnode;
  double tuvs[2*MAXDEG];
  int curve[MAXDEG];
  Ring *ring0[MAXDEG], *ring1[MAXDEG];
  int triangle;
  int tri0, tri1, tri2;
  double uv0[2], uv1[2], uv2[2];

  if ( !gridValidNode(grid, n0) || !gridValidNode(grid, n1) ) return EMPTY; 
  if ( NULL == gridEquator( grid, n0, n1) ) return EMPTY;

  gap0 = gridEqu(grid,0);
  gap1 = gridEqu(grid,gridNGem(grid));
  if ( gridContinuousEquator(grid) ) {
    printf("%s: %d: gridReconstructSplitEdgeRatio has ContinuousEquator?\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }
  face0 = gridFindFace(grid, n0, n1, gap0 );
  face1 = gridFindFace(grid, n0, n1, gap1 );
  if ( face0 == EMPTY || face1 == EMPTY ) {
    printf("%s: %d: gridReconstructSplitEdgeRatio has empty face %d %d\n",
	   __FILE__,__LINE__,face0,face1);
    return EMPTY;
  }
  gridFace(grid, face0, nodes, &faceId0 );
  gridFace(grid, face1, nodes, &faceId1 );
  gridNodeUV(grid,gap0,faceId0,uvgap0);
  gridNodeUV(grid,gap1,faceId1,uvgap1);

  parent = gridParentGeometry(grid, n0, n1);
  if (parent == 0) return EMPTY;

  tuv0[0] = tuv0[1] = tuv1[0] = tuv1[1] = newtuv[0] = newtuv[1] = DBL_MAX;
  if (parent > 0) {
    gridNodeUV(grid,n0,parent,tuv0);
    gridNodeUV(grid,n1,parent,tuv1);
    newtuv[0] = (1-ratio)*tuv0[0]+ratio*tuv1[0];
    newtuv[1] = (1-ratio)*tuv0[1]+ratio*tuv1[1];
    gridEvaluateOnFace(grid, parent, newtuv, xyz );
  }else{
    gridNodeT(grid,n0,-parent,tuv0);
    gridNodeT(grid,n1,-parent,tuv1);
    newtuv[0] = (1-ratio)*tuv0[0]+ratio*tuv1[0];
    gridEvaluateOnEdge(grid, -parent, newtuv[0], xyz );
  }

  printf("new node at%23.15e%23.15e%23.15e\n",xyz[0],xyz[1],xyz[2]);

  enclosing_cell = gridFindEnclosingCell(grid, n0, xyz );

  if ( EMPTY == enclosing_cell ) {
    printf("%s: %d: gridReconstructSplitEdgeRatio could not find cell.\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }

  newnode = gridSplitCellAt(grid, enclosing_cell, xyz );
  if ( newnode == EMPTY ) return EMPTY;

  printf("new node %d inserted into cell %d\n",newnode,enclosing_cell);
  
  nnode=1;
  tuvs[0+2*(nnode-1)] = tuv0[0];
  tuvs[1+2*(nnode-1)] = tuv0[1];
  curve[(nnode-1)] = n0;
  
  while ( newnode != curve[nnode-1] ) {
    if (grid != gridThreadCurveThroughVolume(grid, parent, 
					     curve[(nnode-1)],
					     &(tuvs[2*(nnode-1)]),
					     newnode,   newtuv,
					     &(curve[nnode]),
					     &(tuvs[2*nnode]) ) ) {
      printf("%s: %d: gridThreadCurveThroughVolume failed n0\n",
	     __FILE__,__LINE__);
      for (node=0;node<(nnode-1);node++) {
	gridWriteTecplotEquator(grid, curve[node], curve[node+1],
				"edge_split_equator.t");
      }
      return EMPTY;
    }
    nnode++;
    if (nnode >= MAXDEG ) {
      printf("%s: %d: gridReconstructSplitEdgeRatio nnode exhasted MAXDEG\n",
	     __FILE__,__LINE__);
      return EMPTY;
    }
  }
  while ( n1 != curve[nnode-1] ) {
    if (grid != gridThreadCurveThroughVolume(grid, parent, 
					     curve[(nnode-1)],
					     &(tuvs[2*(nnode-1)]),
					     n1,   tuv1,
					     &(curve[nnode]),
					     &(tuvs[2*nnode]) ) ) {
      printf("%s: %d: gridThreadCurveThroughVolume failed n1\n",
	     __FILE__,__LINE__);
      for (node=0;node<(nnode-1);node++) {
	gridWriteTecplotEquator(grid, curve[node], curve[node+1],
				"edge_split_equator.t");
      }
      return EMPTY;
    }
    nnode++;
    if (nnode >= MAXDEG ) {
      printf("%s: %d: gridReconstructSplitEdgeRatio nnode exhasted MAXDEG\n",
	     __FILE__,__LINE__);
      return EMPTY;
    }
  }

  ring0[0] = ringCreate(  );
  for (node=0;node<(nnode-1);node++) {
    ringAddSegment( ring0[0], curve[node], curve[node+1], 
		    &(tuvs[2*node]), &(tuvs[2*(node+1)]) );
  }
  ringAddSegment( ring0[0], curve[nnode-1], gap0, 
		  &(tuvs[2*(nnode-1)]), uvgap0 );
  ringAddSegment( ring0[0], gap0, curve[0], 
		  uvgap0, tuvs );
  
  if (grid != gridFillRingWithFaces(grid, ring0[0], faceId0 ) ){
    printf("%s: %d: gridFillRingWithFaces ring0[0] returned NULL.\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }

  ring1[0] = ringCreate(  );
  for (node=(nnode-2);node >= 0;node--) {
    ringAddSegment( ring1[0], curve[node+1], curve[node], 
		    &(tuvs[2*(node+1)]), &(tuvs[2*node]) );
  }
  ringAddSegment( ring1[0], curve[0], gap1, 
		  tuvs, uvgap1 );
  ringAddSegment( ring1[0], gap1, curve[nnode-1], 
		  uvgap1, &(tuvs[2*(nnode-1)]) );
  
  if (grid != gridFillRingWithFaces(grid, ring1[0], faceId1 ) ) {
    printf("%s: %d: gridFillRingWithFaces ring1[0] returned NULL.\n",
	   __FILE__,__LINE__);
    return EMPTY;
  }
  
  for ( triangle=0 ; triangle < ringTriangles(ring0[0]) ; triangle++ ) {
    if ( ring0[0] != ringTriangle( ring0[0], triangle,
				&tri0, &tri1, &tri2,
				uv0, uv1, uv2 ) ) return EMPTY;
    gridAddFaceUV( grid, 
		   tri0, uv0[0], uv0[1],
		   tri1, uv1[0], uv1[1],
		   tri2, uv2[0], uv2[1],
		   faceId0 );
  }
  gridRemoveFace(grid,face0);

  for ( triangle=0 ; triangle < ringTriangles(ring1[0]) ; triangle++ ) {
    if ( ring1[0] != ringTriangle( ring1[0], triangle,
				&tri0, &tri1, &tri2,
				uv0, uv1, uv2 ) ) return EMPTY;
    gridAddFaceUV( grid, 
		   tri0, uv0[0], uv0[1],
		   tri1, uv1[0], uv1[1],
		   tri2, uv2[0], uv2[1],
		   faceId1 );
  }
  gridRemoveFace(grid,face1);

  printf("n0 %d n1 %d\n",n0,n1);
  if ( 29==n0 && 207==n1) {
    return EMPTY;
  }

  for ( triangle=0 ; triangle < ringTriangles(ring0[0]) ; triangle++ ) {
    if ( ring0[0] != ringTriangle( ring0[0], triangle,
				&tri0, &tri1, &tri2,
				uv0, uv1, uv2 ) ) return EMPTY;
    gridRemoveCellsOutsideOfFaces( grid, tri1, tri0, tri2 );
  }
  for ( triangle=0 ; triangle < ringTriangles(ring1[0]) ; triangle++ ) {
    if ( ring1[0] != ringTriangle( ring1[0], triangle,
				&tri0, &tri1, &tri2,
				uv0, uv1, uv2 ) ) return EMPTY;
    gridRemoveCellsOutsideOfFaces( grid, tri1, tri0, tri2 );
  }

  ringFree( ring0[0] );
  ringFree( ring1[0] );

  return newnode;
}

int gridSplitEdge(Grid *grid, int n0, int n1)
{
  return gridSplitEdgeRatio(grid, NULL, n0, n1, 0.5 );
}

int gridSplitEdgeRatio(Grid *grid, Queue *queue, int n0, int n1, double ratio )
{
  int igem, cell, nodes[4], inode, node;
  double xyz0[3], xyz1[3], xyz[3], dummy_xyz[3];
  int newnode, newnodes0[4], newnodes1[4];
  int gap0, gap1, face0, face1, faceNodes0[3], faceNodes1[3], faceId0, faceId1;
  int newface_gap0n0, newface_gap0n1, newface_gap1n0, newface_gap1n1;
  int edge, edgeId;
  int newedge0, newedge1;
  double t0, t1, newT;
  double minAR;

  if ( !gridValidNode(grid, n0) || !gridValidNode(grid, n1) ) return EMPTY; 
  if ( NULL == gridEquator( grid, n0, n1) ) return EMPTY;

  /* If the equator has a gap find the faces to be split or return */
  gap0 = gridEqu(grid,0);
  gap1 = gridEqu(grid,gridNGem(grid));
  face0 = face1 = EMPTY;
  if ( !gridContinuousEquator(grid) ){
    face0 = gridFindFace(grid, n0, n1, gap0 );
    face1 = gridFindFace(grid, n0, n1, gap1 );
    if ( face0 == EMPTY || face1 == EMPTY ) return EMPTY;
  }

  /* create new node and initialize */
  if (grid != gridNodeXYZ(grid, n0, xyz0) ) return EMPTY;
  if (grid != gridNodeXYZ(grid, n1, xyz1) ) return EMPTY;
  for (inode = 0 ; inode < 3 ; inode++) 
    xyz[inode] = (1-ratio)*xyz0[inode] + ratio*xyz1[inode]; 
  newnode = gridAddNode(grid, xyz[0], xyz[1], xyz[2] );
  if ( newnode == EMPTY ) return EMPTY;
  gridSetMapMatrixToAverageOfNodes2(grid, newnode, n0, n1 );
  gridSetAuxToAverageOfNodes2(grid, newnode, n0, n1 );

  /* insert new edges to use for projection and validity check */
  newedge0 = newedge1 = EMPTY;
  edge = gridFindEdge(grid,n0,n1);
  if ( edge != EMPTY ) {
    edgeId = gridEdgeId(grid,n0,n1);
    gridNodeT(grid,n0,edgeId,&t0);
    gridNodeT(grid,n1,edgeId,&t1);
    newT = (1-ratio)*t0+ratio*t1;

    if ( gridSurfaceNodeConstrained(grid) ) {
      gridEvaluateOnEdge(grid, edgeId, newT, xyz );
      gridSetNodeXYZ(grid, newnode, xyz);
    }

    newedge0 = gridAddEdgeAndQueue(grid,queue,n0,newnode,edgeId,t0,newT);
    newedge1 = gridAddEdgeAndQueue(grid,queue,n1,newnode,edgeId,t1,newT);
  }

  /* insert new faces to use for projection and validity check */
  newface_gap0n0 = newface_gap0n1 = EMPTY;
  newface_gap1n0 = newface_gap1n1 = EMPTY;
  if ( !gridContinuousEquator(grid) ){
    double n0Id0uv[2], n1Id0uv[2], n0Id1uv[2], n1Id1uv[2];
    double gap0uv[2], gap1uv[2], newId0uv[2], newId1uv[2]; 
    gridFace(grid,face0,faceNodes0,&faceId0);
    gridFace(grid,face1,faceNodes1,&faceId1);
    gridNodeUV(grid,n0,faceId0,n0Id0uv);
    gridNodeUV(grid,n1,faceId0,n1Id0uv);
    gridNodeUV(grid,n0,faceId1,n0Id1uv);
    gridNodeUV(grid,n1,faceId1,n1Id1uv);
    gridNodeUV(grid,gap0,faceId0,gap0uv);
    gridNodeUV(grid,gap1,faceId1,gap1uv);
    newId0uv[0] = (1-ratio)*n0Id0uv[0] + ratio*n1Id0uv[0];
    newId0uv[1] = (1-ratio)*n0Id0uv[1] + ratio*n1Id0uv[1];
    newId1uv[0] = (1-ratio)*n0Id1uv[0] + ratio*n1Id1uv[0];
    newId1uv[1] = (1-ratio)*n0Id1uv[1] + ratio*n1Id1uv[1];

    if ( gridSurfaceNodeConstrained(grid) ) {
      if ( EMPTY != edge ) {
	/* update uv parameters only for faces next to geom edge */
	gridResolveOnFace(grid, faceId0, newId0uv, xyz, dummy_xyz);
	gridResolveOnFace(grid, faceId1, newId1uv, xyz, dummy_xyz);
      } else {
	/* assume id0==id1 or fake geometry */
	gridEvaluateOnFace(grid, faceId0, newId0uv, xyz);
	gridSetNodeXYZ(grid, newnode, xyz);
      }
    }

    newface_gap0n0 = gridAddFaceUVAndQueue(grid, queue,
					   n0, n0Id0uv[0], n0Id0uv[1],
					   newnode, newId0uv[0],newId0uv[1],
					   gap0, gap0uv[0], gap0uv[1],
					   faceId0 );
    newface_gap0n1 = gridAddFaceUVAndQueue(grid, queue,
					   n1, n1Id0uv[0], n1Id0uv[1], 
					   gap0, gap0uv[0], gap0uv[1], 
					   newnode, newId0uv[0], newId0uv[1], 
					   faceId0 );
    newface_gap1n0 = gridAddFaceUVAndQueue(grid, queue,
					   n0, n0Id1uv[0], n0Id1uv[1], 
					   gap1, gap1uv[0], gap1uv[1], 
					   newnode, newId1uv[0], newId1uv[1], 
					   faceId1 );
    newface_gap1n1 = gridAddFaceUVAndQueue(grid, queue,
					   n1, n1Id1uv[0], n1Id1uv[1], 
					   newnode, newId1uv[0], newId1uv[1], 
					   gap1, gap1uv[0], gap1uv[1], 
					   faceId1 );
  }

  /* find the worst cell */
  minAR = 2.0; 
  for ( igem=0 ; igem<gridNGem(grid) ; igem++ ){
    cell = gridGem(grid,igem);
    gridCell(grid, cell, nodes);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      node = nodes[inode];
      newnodes0[inode]=node;
      newnodes1[inode]=node;
      if ( node == n0 ) newnodes0[inode] = newnode;
      if ( node == n1 ) newnodes1[inode] = newnode;
    }
    minAR = MIN(minAR,gridAR(grid,newnodes0));
    minAR = MIN(minAR,gridAR(grid,newnodes1));
  }

  /* if the worst cell is not good enough then undo the split and return */
  if (minAR < gridMinInsertCost(grid) ) {
    /* the remove and queue methods test for EMPTY==target */
    gridRemoveFaceAndQueue(grid, queue, newface_gap0n0 );
    gridRemoveFaceAndQueue(grid, queue, newface_gap0n1 );
    gridRemoveFaceAndQueue(grid, queue, newface_gap1n0 );
    gridRemoveFaceAndQueue(grid, queue, newface_gap1n1 );
    gridRemoveEdgeAndQueue(grid, queue, newedge0 );
    gridRemoveEdgeAndQueue(grid, queue, newedge1 );
    queueResetCurrentTransaction(queue);

    gridRemoveNode(grid,newnode);
    return EMPTY;
  }

  /* update cell connectivity */
  for ( igem=0 ; igem<gridNGem(grid) ; igem++ ){
    cell = gridGem(grid,igem);
    gridCell(grid, cell, nodes);
    gridRemoveCellAndQueue(grid, queue, cell);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      node = nodes[inode];
      newnodes0[inode]=node;
      newnodes1[inode]=node;
      if ( node == n0 ) newnodes0[inode] = newnode;
      if ( node == n1 ) newnodes1[inode] = newnode;
    }
    gridAddCellAndQueue( grid, queue,
			 newnodes0[0],newnodes0[1],newnodes0[2],newnodes0[3]);
    gridAddCellAndQueue( grid, queue,
			 newnodes1[0],newnodes1[1],newnodes1[2],newnodes1[3]);
  }

  /* remove original faces and edge */
  gridRemoveFaceAndQueue(grid, queue, face0 );
  gridRemoveFaceAndQueue(grid, queue, face1 );
  gridRemoveEdgeAndQueue(grid,queue,edge);

  return newnode;
}

int gridSplitEdgeForce(Grid *grid, Queue *queue, int n0, int n1,
		       GridBool debug_split )
{
  int igem, cell, nodes[4], inode, node;
  double xyz0[3], xyz1[3], xyz[3], dummy_xyz[3];
  int newnode, newnodes0[4], newnodes1[4];
  int gap0, gap1, face0, face1, faceNodes0[3], faceNodes1[3], faceId0, faceId1;
  int newface_gap0n0, newface_gap0n1, newface_gap1n0, newface_gap1n1;
  int edge, edgeId;
  int newedge0, newedge1;
  double t0, t1, newT;
  double minAR, minJac;
  double ratio;
  AdjIterator it;

  ratio = 0.5;

  if ( !gridValidNode(grid, n0) || !gridValidNode(grid, n1) ) return EMPTY; 
  if ( NULL == gridEquator( grid, n0, n1) ) return EMPTY;

  /* If the equator has a gap find the faces to be split or return */
  gap0 = gridEqu(grid,0);
  gap1 = gridEqu(grid,gridNGem(grid));
  face0 = face1 = EMPTY;
  if ( !gridContinuousEquator(grid) ){
    face0 = gridFindFace(grid, n0, n1, gap0 );
    face1 = gridFindFace(grid, n0, n1, gap1 );
    if ( face0 == EMPTY || face1 == EMPTY ) return EMPTY;
  }

  /* create new node and initialize */
  if (grid != gridNodeXYZ(grid, n0, xyz0) ) return EMPTY;
  if (grid != gridNodeXYZ(grid, n1, xyz1) ) return EMPTY;
  for (inode = 0 ; inode < 3 ; inode++) 
    xyz[inode] = (1-ratio)*xyz0[inode] + ratio*xyz1[inode]; 
  newnode = gridAddNode(grid, xyz[0], xyz[1], xyz[2] );
  if ( newnode == EMPTY ) return EMPTY;
  gridSetMapMatrixToAverageOfNodes2(grid, newnode, n0, n1 );
  gridSetAuxToAverageOfNodes2(grid, newnode, n0, n1 );

  /* insert new edges to use for projection and validity check */
  newedge0 = newedge1 = EMPTY;
  edge = gridFindEdge(grid,n0,n1);
  if ( edge != EMPTY ) {
    edgeId = gridEdgeId(grid,n0,n1);
    gridNodeT(grid,n0,edgeId,&t0);
    gridNodeT(grid,n1,edgeId,&t1);
    newT = (1-ratio)*t0+ratio*t1;

    if ( gridSurfaceNodeConstrained(grid) ) {
      gridEvaluateOnEdge(grid, edgeId, newT, xyz );
      gridSetNodeXYZ(grid, newnode, xyz);
    }

    newedge0 = gridAddEdgeAndQueue(grid,queue,n0,newnode,edgeId,t0,newT);
    newedge1 = gridAddEdgeAndQueue(grid,queue,n1,newnode,edgeId,t1,newT);
  }

  /* insert new faces to use for projection and validity check */
  newface_gap0n0 = newface_gap0n1 = EMPTY;
  newface_gap1n0 = newface_gap1n1 = EMPTY;
  if ( !gridContinuousEquator(grid) ){
    double n0Id0uv[2], n1Id0uv[2], n0Id1uv[2], n1Id1uv[2];
    double gap0uv[2], gap1uv[2], newId0uv[2], newId1uv[2]; 
    gridFace(grid,face0,faceNodes0,&faceId0);
    gridFace(grid,face1,faceNodes1,&faceId1);
    gridNodeUV(grid,n0,faceId0,n0Id0uv);
    gridNodeUV(grid,n1,faceId0,n1Id0uv);
    gridNodeUV(grid,n0,faceId1,n0Id1uv);
    gridNodeUV(grid,n1,faceId1,n1Id1uv);
    gridNodeUV(grid,gap0,faceId0,gap0uv);
    gridNodeUV(grid,gap1,faceId1,gap1uv);
    newId0uv[0] = (1-ratio)*n0Id0uv[0] + ratio*n1Id0uv[0];
    newId0uv[1] = (1-ratio)*n0Id0uv[1] + ratio*n1Id0uv[1];
    newId1uv[0] = (1-ratio)*n0Id1uv[0] + ratio*n1Id1uv[0];
    newId1uv[1] = (1-ratio)*n0Id1uv[1] + ratio*n1Id1uv[1];

    if ( gridSurfaceNodeConstrained(grid) ) {
      if ( EMPTY != edge ) {
	/* update uv parameters only for faces next to geom edge */
	gridResolveOnFace(grid, faceId0, newId0uv, xyz, dummy_xyz);
	gridResolveOnFace(grid, faceId1, newId1uv, xyz, dummy_xyz);
      } else {
	/* assume id0==id1 or fake geometry */
	gridEvaluateOnFace(grid, faceId0, newId0uv, xyz);
	gridSetNodeXYZ(grid, newnode, xyz);
      }
    }

    newface_gap0n0 = gridAddFaceUVAndQueue(grid, queue,
					   n0, n0Id0uv[0], n0Id0uv[1],
					   newnode, newId0uv[0],newId0uv[1],
					   gap0, gap0uv[0], gap0uv[1],
					   faceId0 );
    newface_gap0n1 = gridAddFaceUVAndQueue(grid, queue,
					   n1, n1Id0uv[0], n1Id0uv[1], 
					   gap0, gap0uv[0], gap0uv[1], 
					   newnode, newId0uv[0], newId0uv[1], 
					   faceId0 );
    newface_gap1n0 = gridAddFaceUVAndQueue(grid, queue,
					   n0, n0Id1uv[0], n0Id1uv[1], 
					   gap1, gap1uv[0], gap1uv[1], 
					   newnode, newId1uv[0], newId1uv[1], 
					   faceId1 );
    newface_gap1n1 = gridAddFaceUVAndQueue(grid, queue,
					   n1, n1Id1uv[0], n1Id1uv[1], 
					   newnode, newId1uv[0], newId1uv[1], 
					   gap1, gap1uv[0], gap1uv[1], 
					   faceId1 );
  }

  /* insert cells */
  for ( igem=0 ; igem<gridNGem(grid) ; igem++ ){
    cell = gridGem(grid,igem);
    gridCell(grid, cell, nodes);
    for ( inode = 0 ; inode < 4 ; inode++ ){
      node = nodes[inode];
      newnodes0[inode]=node;
      newnodes1[inode]=node;
      if ( node == n0 ) newnodes0[inode] = newnode;
      if ( node == n1 ) newnodes1[inode] = newnode;
    }
    gridAddCellAndQueue( grid, queue,
			 newnodes0[0],newnodes0[1],newnodes0[2],newnodes0[3]);
    gridAddCellAndQueue( grid, queue,
			 newnodes1[0],newnodes1[1],newnodes1[2],newnodes1[3]);
  }

   if ( gridContinuousEquator(grid) ){
     gridSmoothNodeVolumeSimplex( grid, newnode );
   }else{
     /* will not work between faces */
     gridSmoothNodeVolumeUVSimplex( grid, newnode );
   }
  gridNodeMinCellJacDet2(grid, newnode, &minJac );
  gridNodeAR(grid, newnode, &minAR );
  //  printf("min AR%20.15f Jac%20.15f\n",minAR, minJac);

  if (minAR < gridMinInsertCost(grid) ) {

    if (debug_split) {
      for ( it = adjFirst(gridCellAdj(grid),newnode); 
	    adjValid(it); 
	    it = adjNext(it) ) {
	cell = adjItem(it);
	gridWriteTecplotCellJacDet(grid,cell,NULL);
      }
    }

    it = adjFirst(gridCellAdj(grid),newnode);
    while (adjValid(it)){
      cell = adjItem(it);
      gridRemoveCellAndQueue( grid, queue, cell );
      it = adjFirst(grid->cellAdj,newnode);
    }
    /* the remove and queue methods test for EMPTY==target */
    gridRemoveFaceAndQueue(grid, queue, newface_gap0n0 );
    gridRemoveFaceAndQueue(grid, queue, newface_gap0n1 );
    gridRemoveFaceAndQueue(grid, queue, newface_gap1n0 );
    gridRemoveFaceAndQueue(grid, queue, newface_gap1n1 );
    gridRemoveEdgeAndQueue(grid, queue, newedge0 );
    gridRemoveEdgeAndQueue(grid, queue, newedge1 );
    gridRemoveNode(grid,newnode);
    queueResetCurrentTransaction(queue);

    return EMPTY;
  }else{

    /* remove old cells */
    for ( igem=0 ; igem<gridNGem(grid) ; igem++ ){
      cell = gridGem(grid,igem);
      gridCell(grid, cell, nodes);
      gridRemoveCellAndQueue(grid, queue, cell);
    }

    /* remove original faces and edge */
    gridRemoveFaceAndQueue(grid, queue, face0 );
    gridRemoveFaceAndQueue(grid, queue, face1 );
    gridRemoveEdgeAndQueue(grid,queue,edge);

    return newnode;
  }
}

int gridSplitEdgeRepeat(Grid *grid, Queue *queue, int n0, int n1, 
			GridBool debug_split  )
{
  int newnode, newnode0, newnode1;
  double existing_min;
  int try;
  int iequ, nequ;
  int node, nodes;
  int worstnode;
  double jacdet, worstjacdet;
  

  newnode = gridSplitEdgeForce( grid, queue, n0, n1, FALSE );
  if (EMPTY != newnode) return newnode;

  /* punch out if this is not a critical edge */
  if (0 == gridParentGeometry(grid, n0, n1) ) return EMPTY;

  existing_min = gridMinInsertCost(grid);
  gridSetMinInsertCost(grid,1.0e-4);
  
  newnode = gridSplitEdgeForce( grid, queue, n0, n1, debug_split );
  gridSetMinInsertCost(grid,existing_min);

  return newnode;
}

int gridSplitEdgeIfNear(Grid *grid, int n0, int n1, double *xyz)
{
  int i;
  int newnode;
  double oldXYZ[3], xyz0[3], xyz1[3];
  double edgeXYZ[3], edgeLength, edgeDir[3];
  double newEdge[3], edgePosition, radius, radiusVector[3];

  gridNodeXYZ(grid, n0, xyz0);
  gridNodeXYZ(grid, n1, xyz1);
  for (i=0;i<3;i++) edgeXYZ[i] = xyz1[i]   - xyz0[i];
  edgeLength = sqrt ( edgeXYZ[0]*edgeXYZ[0] + 
		      edgeXYZ[1]*edgeXYZ[1] +
		      edgeXYZ[2]*edgeXYZ[2] );
  for (i=0;i<3;i++) edgeDir[i] = edgeXYZ[i]/edgeLength;
  for (i=0;i<3;i++) newEdge[i] = xyz[i] - xyz0[i];
  edgePosition =  ( newEdge[0]*edgeDir[0] + 
		    newEdge[1]*edgeDir[1] + 
		    newEdge[2]*edgeDir[2] ) / edgeLength;
  for (i=0;i<3;i++) radiusVector[i] = newEdge[i] -
		      edgePosition*edgeLength*edgeDir[i];
  radius = sqrt ( radiusVector[0]*radiusVector[0] + 
		  radiusVector[1]*radiusVector[1] + 
		  radiusVector[2]*radiusVector[2] );
      
  if ( edgePosition > -0.05 && edgePosition < 0.05 && 
       radius < 0.05*edgeLength) {
    newnode = n0;
    gridNodeXYZ(grid, newnode, oldXYZ );
    gridSetNodeXYZ(grid, newnode, xyz );
    if ( gridNegCellAroundNode(grid, newnode ) ){
      gridSetNodeXYZ(grid, newnode, oldXYZ );
    }else{
      return newnode;
    }
  }

  if ( edgePosition > 0.95 && edgePosition < 1.05 && 
       radius < 0.05*edgeLength) {
    newnode = n1;
    gridNodeXYZ(grid, newnode, oldXYZ );
    gridSetNodeXYZ(grid, newnode, xyz );
    if ( gridNegCellAroundNode(grid, newnode ) ){
      gridSetNodeXYZ(grid, newnode, oldXYZ );
    }else{
      return newnode;
    }
  }

  if ( edgePosition > 0.0 && edgePosition < 1.0 && 
       radius < 0.001*edgeLength) {
    newnode = gridSplitEdgeRatio(grid, NULL, n0, n1, edgePosition);
    return newnode;
  }

  return EMPTY;
}

int gridSplitFaceAt(Grid *grid, int *face_nodes, double *xyz)
{
  int newnode;
  int face, faceId, cell0, cell1;
  int nodes[4], nodes0[4], nodes1[4], newnodes[4], nodetemp;
  double xyz0[3], xyz1[3], xyz2[3], bary[3];
  double U[3], V[3], baryU, baryV, newU[3], newV[3];
  int n, i;

  cell0=gridFindOtherCellWith3Nodes(grid,
				    face_nodes[0],face_nodes[1],face_nodes[2],
				    EMPTY );
  if (EMPTY==cell0) return EMPTY;

  cell1=gridFindOtherCellWith3Nodes(grid,
				    face_nodes[0],face_nodes[1],face_nodes[2],
				    cell0 );

  face = gridFindFace(grid, face_nodes[0], face_nodes[1], face_nodes[2] );
  if ( EMPTY == cell1 && EMPTY == face ) return EMPTY;
  faceId = EMPTY;
  gridFace(grid, face, nodes, &faceId); // only need to set faceId

  //printf("gridSplitFaceAt faceId %d\n",faceId);

  /* set up nodes so that nodes[0]-nodes[2] is face 
     and nodes[0]-nodes[3] is the cell
     for the case of nodes0 or nodes0 and nodes1 */

  if (grid != gridCell(grid, cell0, nodes) ) return EMPTY;
  nodes0[0] = face_nodes[0]; nodes0[1] = face_nodes[1];
  if (grid != gridOrient(grid, nodes, nodes0 ) ) return EMPTY;
  if ( nodes0[2] != face_nodes[2] ) {
    nodetemp = nodes0[2];
    nodes0[2] = nodes0[3];
    nodes0[3] = nodetemp;
    nodetemp = nodes0[0];
    nodes0[0] = nodes0[1];
    nodes0[1] = nodetemp;
  }

  if (grid == gridCell(grid, cell1, nodes) ) {
    nodes1[0] = face_nodes[0]; nodes1[1] = face_nodes[1];
    if (grid != gridOrient(grid, nodes, nodes1 ) ) return EMPTY;
    if ( nodes1[2] != face_nodes[2] ) {
      nodetemp = nodes1[2];
      nodes1[2] = nodes1[3];
      nodes1[3] = nodetemp;
      nodetemp = nodes1[0];
      nodes1[0] = nodes1[1];
      nodes1[1] = nodetemp;
    }
  }

  // if (0.0>=gridVolume(grid, nodes? ) ) return EMPTY;
  
  newnode = gridAddNode(grid, xyz[0], xyz[1], xyz[2] );
  if ( newnode == EMPTY ) return EMPTY;

  gridSetMapMatrixToAverageOfNodes3(grid, newnode, 
				    face_nodes[0],face_nodes[1],face_nodes[2] );
  gridSetAuxToAverageOfNodes3(grid, newnode,
			      face_nodes[0], face_nodes[1], face_nodes[2] );

  if ( EMPTY != face ) {
    gridNodeXYZ(grid, face_nodes[0], xyz0);
    gridNodeXYZ(grid, face_nodes[1], xyz1);
    gridNodeXYZ(grid, face_nodes[2], xyz2);
    gridTriangularBarycentricCoordinate3D( xyz0, xyz1, xyz2, xyz, bary );
    for (i=0;i<3;i++) {
      U[i] = gridNodeU(grid, nodes0[i], faceId);
      V[i] = gridNodeV(grid, nodes0[i], faceId);
    }
    baryU = gridDotProduct(bary,U);
    baryV = gridDotProduct(bary,V);

    if ( grid != gridRemoveFace(grid, face ) ) return EMPTY;

    for (i=0;i<3;i++){
      for (n=0;n<4;n++) newnodes[n] = nodes0[n];
      newnodes[i] = newnode; 
      for (n=0;n<3;n++) newU[n] = U[n];
      for (n=0;n<3;n++) newV[n] = V[n];
      newU[i] = baryU;
      newV[i] = baryV;
      if (EMPTY == gridAddFaceUV(grid, 
				 newnodes[0], newU[0], newV[0], 
				 newnodes[1], newU[1], newV[1],  
				 newnodes[2], newU[2], newV[2],  
				 faceId ) ) return EMPTY;
    }

    if (grid!=gridProjectNodeToFace(grid, newnode, faceId )) return EMPTY;

  }

  if ( grid != gridRemoveCell(grid, cell0 ) ) return EMPTY;
  for (i=0;i<3;i++){
    for (n=0;n<4;n++) newnodes[n] = nodes0[n];
    newnodes[i] = newnode; 
    if (EMPTY == gridAddCell(grid, newnodes[0], newnodes[1], 
			     newnodes[2], newnodes[3] ) ) return EMPTY;
  }

  if ( EMPTY != cell1 ) {
    if ( grid != gridRemoveCell(grid, cell1 ) ) return EMPTY;
    for (i=0;i<3;i++){
      for (n=0;n<4;n++) newnodes[n] = nodes1[n];
      newnodes[i] = newnode; 
      if (EMPTY == gridAddCell(grid, newnodes[0], newnodes[1], 
			       newnodes[2], newnodes[3] ) ) return EMPTY;
    }
  }

  if ( gridNegCellAroundNode(grid, newnode) ) {
    gridCollapseEdge(grid, NULL, nodes0[0], newnode, 0.0 );
    return EMPTY;
  }else{
    return newnode;
  }
}

int gridSplitCellAt(Grid *grid, int cell, double *xyz)
{
  int newnode;
  int nodes[4], newnodes[4];
  int n, i;

  if (grid != gridCell(grid, cell, nodes) ) return EMPTY;

  //if (0.0>=gridVolume(grid, nodes ) ) return EMPTY;

  newnode = gridAddNode(grid, xyz[0], xyz[1], xyz[2] );
  if ( newnode == EMPTY ) return EMPTY;

  gridSetMapMatrixToAverageOfNodes4(grid, newnode, 
				    nodes[0], nodes[1], nodes[2], nodes[3] );
  gridSetAuxToAverageOfNodes4(grid, newnode,
			      nodes[0], nodes[1], nodes[2], nodes[3]);

  for (i=0;i<4;i++){
    for (n=0;n<4;n++) newnodes[n] = nodes[n];
    newnodes[i] = newnode; 
    if (1.0e-12 > gridVolume(grid, newnodes ) ) {
      gridRemoveNode(grid, newnode );
      return EMPTY;
    }
  }

  if ( grid != gridRemoveCell(grid, cell ) ) return EMPTY;

  for (i=0;i<4;i++){
    for (n=0;n<4;n++) newnodes[n] = nodes[n];
    newnodes[i] = newnode; 
    if (EMPTY == gridAddCell(grid, newnodes[0], newnodes[1], 
		     	           newnodes[2], newnodes[3] ) ) return EMPTY;
  }

  return newnode;
}

int gridInsertInToGeomEdge(Grid *grid, double *xyz )
{
  int edge, maxedge, edgeId, nodes[2];
  int newnode;

  newnode = EMPTY;
  edge = 0;
  maxedge = gridMaxEdge(grid);
  while ( EMPTY == newnode && edge < maxedge ) {
    if (grid == gridEdge(grid, edge, nodes, &edgeId) ){
      newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[1],xyz);
    }
    edge++;
  }

  if ( newnode != EMPTY ) gridProjectNode(grid, newnode );

  return newnode;
}

int gridInsertInToGeomFace(Grid *grid, double *xyz )
{
  int foundFace;
  int face, maxface, faceId, nodes[3];
  double xyz0[3], xyz1[3], xyz2[3];
  double edge0[3], edge1[3], edge2[3];
  double leg0[3], leg1[3], leg2[3];
  double norm[3], norm0[3], norm1[3], norm2[3];
  double normLength, unit[3], normDistance;
  int newnode;

  GridBool edgeSplit;

  newnode = EMPTY;
  foundFace = EMPTY;
  edgeSplit = FALSE;
  face = 0;
  maxface = gridMaxFace(grid);
  while ( EMPTY == foundFace && !edgeSplit && face < maxface ) {
    if (grid == gridFace(grid, face, nodes, &faceId) ){
      /* first try putting it on an edge */
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[1],xyz);
	edgeSplit = (EMPTY != newnode);
      }
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[1],nodes[2],xyz);
	edgeSplit = (EMPTY != newnode);
      }
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[2],nodes[0],xyz);
	edgeSplit = (EMPTY != newnode);
      }
      if (!edgeSplit) {
	/* then try putting it on the face */
	gridNodeXYZ(grid, nodes[0], xyz0);
	gridNodeXYZ(grid, nodes[1], xyz1);
	gridNodeXYZ(grid, nodes[2], xyz2);
	gridSubtractVector(xyz1, xyz0, edge0);
	gridSubtractVector(xyz2, xyz1, edge1);
	gridSubtractVector(xyz0, xyz2, edge2);
	gridSubtractVector(xyz, xyz0, leg0);
	gridSubtractVector(xyz, xyz1, leg1);
	gridSubtractVector(xyz, xyz2, leg2);
	gridCrossProduct(edge0,edge1,norm);
	gridCrossProduct(edge0,leg0,norm0);
	gridCrossProduct(edge1,leg1,norm1);
	gridCrossProduct(edge2,leg2,norm2);
	normLength = sqrt( gridDotProduct( norm, norm ) );
	unit[0] = norm[0] / normLength;
	unit[1] = norm[1] / normLength;
	unit[2] = norm[2] / normLength;
	normDistance = gridDotProduct( unit, leg0 );
	if (FALSE) {
	  printf("f%d X%8.5f Y%8.5f Z%8.5f %6.3f 0 %6.3f 1 %6.3f 2 %6.3f\n",
		 face, xyz[0], xyz[1], xyz[2], normDistance,
		 gridDotProduct( norm, norm0 ),
		 gridDotProduct( norm, norm1 ),
		 gridDotProduct( norm, norm2 ) );
	}
	if ( gridDotProduct( norm, norm0 ) > 0.00 &&
	     gridDotProduct( norm, norm1 ) > 0.00 &&
	     gridDotProduct( norm, norm2 ) > 0.00 &&
	     ABS(normDistance) < 0.01*normLength ){
	  foundFace = face;
	}
      }
    }
    face++;
  }

  if ( edgeSplit ) return newnode;
  if ( EMPTY == foundFace ) return EMPTY;

  newnode = gridSplitFaceAt(grid, nodes, xyz);
  if (EMPTY == newnode) return EMPTY;

  return newnode;
}

int gridInsertInToVolume(Grid *grid, double *xyz)
{
  int cell, maxcell, nodes[4], foundCell, newnode;
  GridBool edgeSplit;
  double xyz0[3], xyz1[3], xyz2[3], xyz3[3];
  double edge0[3], edge1[3];
  double leg0[3], leg1[3];
  double norm0[3], norm1[3], norm2[3], norm3[3];

  newnode = EMPTY;
  foundCell = EMPTY;
  edgeSplit = FALSE;
  cell = 0;
  maxcell = gridMaxCell(grid);
  while ( EMPTY == foundCell && !edgeSplit && cell < maxcell ) {
    if (grid == gridCell(grid, cell, nodes) ){
      /* first try putting it on an edge */
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[1],xyz);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[2],xyz);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[0],nodes[3],xyz);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[1],nodes[2],xyz);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[1],nodes[3],xyz);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit){
	newnode = gridSplitEdgeIfNear(grid,nodes[2],nodes[3],xyz);
	edgeSplit = (EMPTY != newnode);}
      if (!edgeSplit) {
	/* then try putting it in the cell -or splitting faces?-*/
	gridNodeXYZ(grid, nodes[0], xyz0);
	gridNodeXYZ(grid, nodes[1], xyz1);
	gridNodeXYZ(grid, nodes[2], xyz2);
	gridNodeXYZ(grid, nodes[3], xyz3);

	gridSubtractVector(xyz3, xyz1, edge0);
	gridSubtractVector(xyz2, xyz1, edge1);
	gridCrossProduct(edge0,edge1,norm0);

	gridSubtractVector(xyz2, xyz0, edge0);
	gridSubtractVector(xyz3, xyz0, edge1);
	gridCrossProduct(edge0,edge1,norm1);

	gridSubtractVector(xyz3, xyz0, edge0);
	gridSubtractVector(xyz1, xyz0, edge1);
	gridCrossProduct(edge0,edge1,norm2);

	gridSubtractVector(xyz1, xyz0, edge0);
	gridSubtractVector(xyz2, xyz0, edge1);
	gridCrossProduct(edge0,edge1,norm3);

	gridSubtractVector(xyz, xyz0, leg0);
	gridSubtractVector(xyz, xyz1, leg1);

	if (FALSE) {	
	  printf("c%d X%8.5f Y%8.5f Z%8.5f 0 %6.3f 1 %6.3f 2 %6.3f 3 %6.3f\n",
		 cell, xyz[0], xyz[1], xyz[2],
		 gridDotProduct( leg1, norm0 ),
		 gridDotProduct( leg0, norm1 ),
		 gridDotProduct( leg0, norm2 ),
		 gridDotProduct( leg0, norm3 ) );
	  printf("xyz0 X%8.5f Y%8.5f Z%8.5f\n",xyz0[0],xyz0[1],xyz0[2]);
	  printf("xyz1 X%8.5f Y%8.5f Z%8.5f\n",xyz1[0],xyz1[1],xyz1[2]);
	  printf("xyz2 X%8.5f Y%8.5f Z%8.5f\n",xyz2[0],xyz2[1],xyz2[2]);
	  printf("xyz3 X%8.5f Y%8.5f Z%8.5f\n",xyz3[0],xyz3[1],xyz3[2]);
	  printf("norm0 X%8.5f Y%8.5f Z%8.5f\n",norm0[0],norm0[1],norm0[2]);
	}

	if ( gridDotProduct( leg1, norm0 ) > 0.00 &&
	     gridDotProduct( leg0, norm1 ) > 0.00 &&
	     gridDotProduct( leg0, norm2 ) > 0.00 &&
	     gridDotProduct( leg0, norm3 ) > 0.00 ){
	  foundCell = cell;
	} 
      }
    }
    cell++;
  }

  if ( edgeSplit ) return newnode;

  if (EMPTY == foundCell) return EMPTY;

  return gridSplitCellAt(grid,foundCell,xyz);
}

Grid *gridCollapseEdge(Grid *grid, Queue *queue, int n0, int n1, 
		       double requestedratio )
{
  int i, face0, face1;
  int requiredRatio;
  double ratio;
  double xyz0[3], xyz1[3], xyz[3];
  double dummy_xyz[3];
  GridBool volumeEdge;
  int iequ, equ0, equ1;

  int gap0, gap1, faceId0, faceId1;
  double n0Id0uv[2], n1Id0uv[2], newId0uv[2];
  double n0Id1uv[2], n1Id1uv[2], newId1uv[2]; 
  int edge, edgeId;
  double t0, t1, t;

  if ( !gridValidNode(grid, n0) || !gridValidNode(grid, n1) ) return NULL; 

  /* logic to make sure collapse is valid w.r.t. boundary connectivity */
  if ( gridGeometryNode(grid, n1) ) return NULL;
  if ( gridGeometryEdge(grid, n1) && EMPTY == gridFindEdge(grid, n0, n1) ) 
    return NULL;

  if ( NULL == gridEquator( grid, n0, n1) ) return NULL;

  volumeEdge = gridContinuousEquator(grid);

  if ( volumeEdge && gridGeometryFace(grid, n0) && gridGeometryFace(grid, n1))
    return NULL;

  for(iequ=0;iequ<gridNEqu(grid)-1;iequ++) {
    equ0 = gridEqu(grid,iequ);
    equ1 = gridEqu(grid,iequ+1);
    if (EMPTY != gridFindFace(grid,n0,equ0,equ1) && 
	EMPTY != gridFindFace(grid,n1,equ0,equ1) ) return NULL;
  }

  /* override user specified ratio to keep nodes on boundary */
  ratio = requestedratio;
  
  requiredRatio = EMPTY;
  if ( volumeEdge && gridGeometryFace(grid, n1) ) requiredRatio = 1;    
  if ( volumeEdge && gridGeometryFace(grid, n0) ) requiredRatio = 0;
  if ( gridGeometryEdge(grid, n0) && 
       !gridGeometryEdge(grid, n1) ) requiredRatio = 0;
  if ( gridGeometryNode(grid, n0) ) requiredRatio = 0;
  
  if (0 == requiredRatio) { if (0.99 < ratio) return NULL; ratio = 0.0; }
  if (1 == requiredRatio) { if (0.01 > ratio) return NULL; ratio = 1.0; }

  /* determine new node location from ratio */
  if ( NULL == gridNodeXYZ( grid, n0, xyz0) ) return NULL;
  if ( NULL == gridNodeXYZ( grid, n1, xyz1) ) return NULL;
  for (i=0 ; i<3 ; i++) xyz[i] = (1.0-ratio) * xyz0[i] + ratio * xyz1[i];

  /* interpolate parameters */
  face0 = face1 = faceId0 = faceId1 = EMPTY;
  newId0uv[0] = newId0uv[1] = newId1uv[0] = newId1uv[1] = DBL_MAX;
  edge = edgeId = EMPTY;
  t0 = t1 = DBL_MAX;
  if ( !volumeEdge ) {

    edge = gridFindEdge(grid,n0,n1);
    if ( edge != EMPTY ) {
      edgeId = gridEdgeId(grid,n0,n1);
      gridNodeT(grid, n0, edgeId, &t0 );
      gridNodeT(grid, n1, edgeId, &t1 );
      t = (1.0-ratio) * t0 + ratio * t1;
      if ( gridSurfaceNodeConstrained(grid) ) {
	gridEvaluateOnEdge(grid, edgeId, t, xyz );
      }
      gridSetNodeT(grid, n0, edgeId, t);
      gridSetNodeT(grid, n1, edgeId, t);
    }

    gap0 = gridEqu(grid,0);
    gap1 = gridEqu(grid,gridNGem(grid));
    face0 = gridFindFace(grid, n0, n1, gap0 );
    face1 = gridFindFace(grid, n0, n1, gap1 );
    faceId0 = gridFaceId(grid, n0, n1, gap0 );
    faceId1 = gridFaceId(grid, n0, n1, gap1 );
    if ( faceId0 == EMPTY || faceId1 == EMPTY ) {
      printf("ERROR %s: %d: collapse faceId empty: %d %d\n",__FILE__,__LINE__,
	     faceId0, faceId1);
      return NULL;
    }
    gridNodeUV(grid,n0,faceId0,n0Id0uv);
    gridNodeUV(grid,n1,faceId0,n1Id0uv);
    gridNodeUV(grid,n0,faceId1,n0Id1uv);
    gridNodeUV(grid,n1,faceId1,n1Id1uv);
    newId0uv[0] = (1.0-ratio) * n0Id0uv[0] + ratio * n1Id0uv[0];
    newId0uv[1] = (1.0-ratio) * n0Id0uv[1] + ratio * n1Id0uv[1];
    newId1uv[0] = (1.0-ratio) * n0Id1uv[0] + ratio * n1Id1uv[0];
    newId1uv[1] = (1.0-ratio) * n0Id1uv[1] + ratio * n1Id1uv[1];

    if ( gridSurfaceNodeConstrained(grid) ) {
      if ( EMPTY != edge ) {
	/* update uv parameters only for faces next to geom edge */
	gridResolveOnFace(grid, faceId0, newId0uv, xyz, dummy_xyz);
	gridResolveOnFace(grid, faceId1, newId1uv, xyz, dummy_xyz);
      } else {
	/* assume id0==id1 or fake geometry */
	gridEvaluateOnFace(grid, faceId0, newId0uv, xyz);
      }
    }

    gridSetNodeUV(grid, n0, faceId0, newId0uv[0], newId0uv[1]);
    gridSetNodeUV(grid, n0, faceId1, newId1uv[0], newId1uv[1]);
    gridSetNodeUV(grid, n1, faceId0, newId0uv[0], newId0uv[1]);
    gridSetNodeUV(grid, n1, faceId1, newId1uv[0], newId1uv[1]);
  }

  /* set node location (this will include any surface geometry constraints) */
  gridSetNodeXYZ( grid, n0, xyz);
  gridSetNodeXYZ( grid, n1, xyz); 
 
  /* if this is not a valid configuration set everything back */
  if (( gridMinARAroundNodeExceptGem( grid, n0 ) < gridMinInsertCost(grid) ) || 
      ( gridMinARAroundNodeExceptGem( grid, n1 ) < gridMinInsertCost(grid) ) ||
      ( gridMinARAroundNodeExceptGemRecon( grid, n0, n1 ) < 
	gridMinInsertCost(grid) ) ||
      ( gridMinARAroundNodeExceptGemRecon( grid, n1, n0 ) < 
	gridMinInsertCost(grid) )  ) {
    if ( edgeId != EMPTY ) {
      gridSetNodeT(grid, n0, edgeId, t0 );
      gridSetNodeT(grid, n1, edgeId, t1 );
    }
    if ( EMPTY != faceId0) {
      gridSetNodeUV( grid, n0, faceId0, n0Id0uv[0], n0Id0uv[1]);
      gridSetNodeUV( grid, n1, faceId0, n1Id0uv[0], n1Id0uv[1]);
    }
    if ( EMPTY != faceId1) {
      gridSetNodeUV( grid, n0, faceId1, n0Id1uv[0], n0Id1uv[1]);
      gridSetNodeUV( grid, n1, faceId1, n1Id1uv[0], n1Id1uv[1]);
    }
    gridSetNodeXYZ( grid, n0, xyz0);
    gridSetNodeXYZ( grid, n1, xyz1);
    return NULL;
  }

  gridRemoveGem( grid );
  
  gridReconnectAllCell(grid, n1, n0 );
  if ( !volumeEdge ) {
    gridRemoveFace(grid, face0 );
    gridRemoveFace(grid, face1 );
    gridReconnectAllFace(grid, n1, n0);
  }
  if ( edge != EMPTY ) {
    gridRemoveEdge(grid,edge);
    gridReconnectAllEdge(grid, n1, n0 );
  }

  if ( volumeEdge && gridGeometryFace(grid, n1) ) {
    gridReconnectAllFace(grid, n1, n0 );
    gridReconnectAllEdge(grid, n1, n0 );
    /* CHECK ME - do edges need to be reconnected too??? */
  }

  gridRemoveNode(grid, n1);

  return grid;
}

Grid *gridFreezeGoodNodes(Grid *grid, double goodAR, 
			  double minLength, double maxLength )
{
  int n0, n1;
  double ratio;
  double ar;

  for ( n0=0; n0<gridMaxNode(grid); n0++ ) { 
    if ( gridValidNode( grid, n0) && !gridNodeFrozen( grid, n0 ) ) {
      if ( NULL == gridLargestRatioEdge( grid, n0, &n1, &ratio) ) 
	return NULL;
      if ( ratio < maxLength ) {
	if ( NULL == gridSmallestRatioEdge( grid, n0, &n1, &ratio) ) 
	  return NULL;
	if ( ratio > minLength ) { 
	  gridNodeAR(grid,n0,&ar);
	  if ( grid == gridProjectNode(grid, n0 ) &&
	       ar > goodAR ) gridFreezeNode(grid,n0);
	}
      }
    }
  }
  return grid;
}

Grid *gridVerifyEdgeExists(Grid *grid, int n0, int n1 )
{
  int i0, i1, nodes0[4], nodes1[4];
  GridBool gotIt;
  AdjIterator it0, it1;

  gotIt=gridCellEdge( grid, n0, n1 );

  if( !gotIt && n0 != EMPTY && n1 != EMPTY ) {
    for ( it0 = adjFirst(gridCellAdj(grid),n0); 
	  adjValid(it0) && !gotIt; 
	  it0 = adjNext(it0) ) {
      gridCell( grid, adjItem(it0), nodes0 );
      for(i0=0;i0<4&&!gotIt;i0++){
	for ( it1 = adjFirst(gridCellAdj(grid),nodes0[i0]); 
	      adjValid(it1) && !gotIt; 
	      it1 = adjNext(it1) ) {
	  gridCell( grid, adjItem(it1), nodes1 );
	  for(i1=0;i1<4&&!gotIt;i1++){
	    if (nodes1[i1] == n1 && !gridNodeFrozen( grid, nodes0[i0] )) {
	      gridCollapseEdge(grid, NULL, n0, nodes0[i0], 0.0);
	      gotIt = gridCellEdge( grid, n0, n1 );
	      if (!gotIt) gridCollapseEdge(grid, NULL, n1, nodes0[i0], 0.0);
	      gotIt = gridCellEdge( grid, n0, n1 );
	    }	    
	  }     
	}
      }
    }
  }

  if (gotIt) return grid;
  return NULL;
}

Grid *gridVerifyFaceExists(Grid *grid, int n0, int n1, int n2 )
{
  GridBool gotIt;

  gotIt=gridCellFace( grid, n0, n1, n2 );

  if (gotIt) return grid;
  return NULL;
}

