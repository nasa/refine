
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef __APPLE__       /* Not needed on Mac OS X */
#include <float.h>
#else
#include <values.h>
#endif
#include "intersect.h"
#include "plan.h"
#include "gridshape.h"
#include "gridmath.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "gridinsert.h"
#include "gridcad.h"

static GridBool gridThisEdgeCanBeModifiedInThisPhase( Grid *grid, 
					       int node0, int node1 )
{
  switch( gridPhase(grid) ) {
  case (gridALL_PHASE) :
    return TRUE;
    break;
  case (gridEDGE_PHASE) :
    return ( 0 > gridParentGeometry(grid, node0, node1) );
    break;
  case (gridFACE_PHASE) :
    return ( 0 < gridParentGeometry(grid, node0, node1) );
    break;
  case (gridVOL_PHASE) :
    return ( 0 == gridParentGeometry(grid, node0, node1) );
    break;
  default : 
    printf("%s: %d: gridThisEdgeCanBeModifiedInThisPhase: phase %d unknown?\n",
	   __FILE__,__LINE__,gridPhase(grid));
    return FALSE;
  }
}

Grid *gridAdapt(Grid *grid, double minLength, double maxLength )
{
  int ranking, conn, nodes[2];
  int report, nnodeAdd, nnodeRemove;
  double length, ratio;
  int newnode;
  Plan *plan;

  gridCreateConn(grid);
  plan = planCreate( gridNConn(grid)/2, MAX(gridNConn(grid)/10,1000) );
  for(conn=0;conn<gridNConn(grid);conn++) {
    gridConn2Node(grid,conn,nodes);
    if ( gridThisEdgeCanBeModifiedInThisPhase(grid,nodes[0],nodes[1]) ) {
      length = gridEdgeRatio(grid,nodes[0],nodes[1]);
      if ( length <= minLength ) planAddItemWithPriority( plan, conn, length );
    }
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
    if ( gridThisEdgeCanBeModifiedInThisPhase(grid,nodes[0],nodes[1]) ) {
      length = gridEdgeRatio(grid,nodes[0],nodes[1]);
      if ( length >= maxLength ) planAddItemWithPriority( plan, conn, length );
    }
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
  GridBool undo;

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
	gridEvaluateOnFace(grid, faceId1, newId1uv, xyz);
	if ( faceId0 == faceId1 ) {
	  newId0uv[0] = newId1uv[0];
	  newId0uv[1] = newId1uv[1];
	}
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

  undo = ( (grid != gridInterpolateMap2(grid, n0, n1, ratio, newnode ) ) ||
	   (grid != gridInterpolateAux2(grid, n0, n1, ratio, newnode ) ) );

  if (!undo) {
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
    undo = ( minAR < gridMinInsertCost(grid) );
  }

  /* if the worst cell is not good enough then undo the split and return */
  if ( undo ) {
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
  gridFace(grid, face, nodes, &faceId); /* only need to set faceId */

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
    gridBarycentricCoordinateTri( xyz0, xyz1, xyz2, xyz, bary );
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

/* node1 is removed */
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

  /* no-op for NULL queue, to suppress compiler warnings */
  queueResetCurrentTransaction(queue);

  if ( !gridValidNode(grid, n0) || !gridValidNode(grid, n1) ) return NULL; 

  /* logic to make sure collapse is valid w.r.t. boundary connectivity */
  if ( gridGeometryNode(grid, n1) ) return NULL;

  if ( gridGeometryBetweenFace(grid, n1) && 
       !gridGeometryBetweenFace(grid, n0)) return NULL;

  if ( gridGeometryBetweenFace(grid, n1) && 
       EMPTY == gridFindEdge(grid, n0, n1) ) 
    return NULL;

  if ( gridGeometryEdge(grid, n1) && gridGeometryEdge(grid, n0) &&
       EMPTY == gridEdgeId(grid,n0,n1) ) return NULL;

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
  if ( gridGeometryBetweenFace(grid, n0) && 
       !gridGeometryBetweenFace(grid, n1) ) requiredRatio = 0;
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
	gridEvaluateOnFace(grid, faceId1, newId1uv, xyz);
	if ( faceId0 == faceId1 ) {
	  newId0uv[0] = newId1uv[0];
	  newId0uv[1] = newId1uv[1];
	}
	/* one node on face one node on edge */
	if ( gridGeometryEdge( grid, n0 ) ) {
	  int edgenodes[2];
	  edge = adjItem(adjFirst(gridEdgeAdj(grid),n0));
	  gridEdge( grid, edge, edgenodes, &edgeId );
	  gridNodeT(grid, n0, edgeId, &t );
	  gridEvaluateOnEdge(grid, edgeId, t, xyz );
	  edge = edgeId = EMPTY;
	}
      }
    }

    gridSetNodeUV(grid, n0, faceId0, newId0uv[0], newId0uv[1]);
    gridSetNodeUV(grid, n0, faceId1, newId1uv[0], newId1uv[1]);
    gridSetNodeUV(grid, n1, faceId0, newId0uv[0], newId0uv[1]);
    gridSetNodeUV(grid, n1, faceId1, newId1uv[0], newId1uv[1]);
  }

  /* do not move gometry nodes */
  if ( gridGeometryNode( grid, n0 ) ) gridNodeXYZ( grid, n0, xyz);

  /* set node location (this will include any surface geometry constraints) */
  gridSetNodeXYZ( grid, n0, xyz);
  gridSetNodeXYZ( grid, n1, xyz); 
 
  /* if this is not a valid configuration set everything back */
  if (( gridMinARAroundNodeExceptGem( grid, n0 ) < gridMinInsertCost(grid) ) || 
      ( gridMinARAroundNodeExceptGem( grid, n1 ) < gridMinInsertCost(grid) ) ||
      ( gridMinARAroundNodeExceptGemRecon( grid, n0, n1 ) < 
	gridMinInsertCost(grid) ) ||
      ( gridMinARAroundNodeExceptGemRecon( grid, n1, n0 ) < 
	gridMinInsertCost(grid) ) ||
      !gridReconnectionOfAllFacesOK(grid, n1, n0) ) {
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

  gridInterpolateMap2(grid, n0, n1, ratio, n0);
  gridInterpolateAux2(grid, n0, n1, ratio, n0);

  gridRemoveNode(grid, n1);

  return grid;
}

