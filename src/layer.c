
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <math.h>
#include <values.h>
#include "master_header.h"
#include "layer.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "grid.h"
#include "adj.h"

typedef struct Normal Normal;
struct Normal {
  int constrained;
  int root, tip;
  double direction[3];
  bool terminated;
};

typedef struct Front Front;
struct Front {
  int globalNode[3];
  int normal[3];
  int constrainedSide[3];
  int parentEdge[3];
};

struct Layer {
  Grid *grid;
  int nfront;
  Front *front;
  int nFrontParent, *frontParent;
  int nnormal;
  Normal *normal;
  int *globalNode2Normal;
  int nConstrainingGeometry, *constrainingGeometry;
  Adj *adj;
};

Layer *layerCreate( Grid *grid )
{
  Layer *layer;
  layer = malloc(sizeof(Layer));
  layer->grid = grid;
  gridAttachNodeSorter( grid, layerSortGlobalNodes, layer );
  layer->nfront=0;
  layer->front=NULL;
  layer->nFrontParent=0;
  layer->frontParent=NULL;
  layer->nnormal=0;
  layer->normal=NULL;
  layer->globalNode2Normal=NULL;
  layer->nConstrainingGeometry=0;
  layer->constrainingGeometry=NULL;
  layer->adj=NULL;
  return layer;
}

Layer *formAdvancingFront( Grid *grid, char *project )
{
  Layer *layer;
  int i, nbc, bc[3];
  bool box, plate, om6;
  
  box = (NULL != strstr( project, "box"));
  plate = (NULL != strstr( project, "plate"));
  om6 = (NULL != strstr( project, "om6"));

  if (box) printf("string %s has box.\n",project);
  if (plate) printf("string %s has plate.\n",project);
  if (om6) printf("string %s has om6.\n",project);

  bc[0]=1;
  bc[1]=2;
  bc[2]=3;
  if(box) nbc = 1;
  if(om6) nbc = 2;
  if(plate) nbc = 3;
  /* skip freezing
  printf("freezing distant volume nodes.\n");
  gridFreezeAll(grid);
  for (i=0;i<nbc;i++){  
    printf("thaw bc %d.\n",bc[i]);
    gridThawNearBC(grid,0.5,bc[i]);
    gridFreezeBCFace(grid,bc[i]);
  }
  */
  printf("make advancing layer object.\n");
  layer = layerCreate(grid);
  printf("make advancing layer front.\n");
  layerMakeFront(layer,nbc,bc);
  printf("make advancing layer front normals.\n");
  layerMakeNormal(layer);
  if (box) {
    layerConstrainNormal(layer,-9);
    layerConstrainNormal(layer,-10);
    layerConstrainNormal(layer,-11);
    layerConstrainNormal(layer,-12);
    layerConstrainNormal(layer,3);
    layerConstrainNormal(layer,4);
    layerConstrainNormal(layer,5);
    layerConstrainNormal(layer,6);
  }
  if (plate) {
    layerTerminateNormalWithX(layer,-0.0001);
    layerTerminateNormalWithX(layer,1.0001);

    /*
      face and four edges
      04 13 12 07 11
      05 15 11 02 14
      06 17 14 10 16

      14 20 26 16 08

      07 20 09 19 18
      08 19 04 22 21
      09 24 23 22 05

      13 12 28 24 06
    */
    layerConstrainNormal(layer,-11);
    layerConstrainNormal(layer,-14);
    layerConstrainNormal(layer,-16);

    layerConstrainNormal(layer,-20);

    layerConstrainNormal(layer,-19);
    layerConstrainNormal(layer,-22);
    layerConstrainNormal(layer,-24);

    layerConstrainNormal(layer,-12);
    /* y1 */
    layerConstrainNormal(layer,4);
    layerConstrainNormal(layer,5);
    layerConstrainNormal(layer,6);
    /* y0 */
    layerConstrainNormal(layer,7);
    layerConstrainNormal(layer,8);
    layerConstrainNormal(layer,9);
    /* inout */
    layerConstrainNormal(layer,13);
    layerConstrainNormal(layer,14);
  }
  if (om6) {
    layerConstrainNormal(layer,5);
  }
  printf("make advancing layer front normals visible to front.\n");
  layerVisibleNormals(layer);
  return layer;
}

Grid *layerGrid(Layer *layer)
{
  return layer->grid;
}

void layerFree(Layer *layer)
{
  gridDetachNodeSorter( layer->grid );
  if (layer->adj != NULL) adjFree(layer->adj);
  if (layer->constrainingGeometry != NULL) free(layer->constrainingGeometry);
  if (layer->globalNode2Normal != NULL) free(layer->globalNode2Normal);
  if (layer->normal != NULL) free(layer->normal);
  if (layer->frontParent != NULL) free(layer->frontParent);
  if (layer->front != NULL) free(layer->front);
  free(layer);
}

void layerSortGlobalNodes(void *voidLayer, int *o2n)
{
  Layer *layer = (Layer *)voidLayer;
  int i, front, normal;

  for (front = 0 ; front < layerNFront(layer) ; front++ )
    for (i=0;i<3;i++) if (EMPTY != layer->front[front].globalNode[i])
      layer->front[front].globalNode[i] = 
	o2n[layer->front[front].globalNode[i]];
  
  for (normal = 0 ; normal < layerNNormal(layer) ; normal++ ) {
    if (EMPTY != layer->normal[normal].root)
      layer->normal[normal].root = o2n[layer->normal[normal].root];
    if (EMPTY != layer->normal[normal].tip)
      layer->normal[normal].tip = o2n[layer->normal[normal].tip];
  }

}

int layerNFront(Layer *layer)
{
  return layer->nfront;
}

int layerNNormal(Layer *layer)
{
  return layer->nnormal;
}

int layerMaxNode(Layer *layer)
{
  return gridMaxNode(layer->grid);
}

Layer *layerMakeFront(Layer *layer, int nbc, int *bc)
{
  int i, ibc, face, nface, id, nodes[3], ifront;

  layer->nfront=0;

  layer->nFrontParent = nbc;
  layer->frontParent = malloc( layer->nFrontParent * sizeof(int) );
  for(ibc=0;ibc<layer->nFrontParent;ibc++) layer->frontParent[ibc] = bc[ibc];

  for (ibc=0;ibc<nbc;ibc++){
    nface =0;
    for(face=0;face<gridMaxFace(layer->grid);face++){
      if (layer->grid == gridFace(layer->grid,face,nodes,&id) &&
	  id==bc[ibc] ) nface++;
    }
    layer->nfront += nface;
    //printf("boundary %d with %d faces added to the %d face front.\n",
    //   bc[ibc],nface,layer->nfront);
  }
  
  layer->front = malloc( 3 * layer->nfront * sizeof(Front) );

  ifront =0;
  for (ibc=0;ibc<nbc;ibc++){
    for(face=0;face<gridMaxFace(layer->grid);face++){
      if (layer->grid == gridFace(layer->grid,face,nodes,&id) &&
	  id==bc[ibc] ) {
	layer->front[ifront].globalNode[0] = nodes[0];
	layer->front[ifront].globalNode[1] = nodes[1];
	layer->front[ifront].globalNode[2] = nodes[2];
	for (i=0;i<3;i++){
	  layer->front[ifront].normal[i] = EMPTY;
	  layer->front[ifront].constrainedSide[i] = 0;
	  layer->front[ifront].parentEdge[i] = 0;
	}
	ifront++;
      }
    }
  }
  
  return layer;
}

bool layerParentFace(Layer *layer, int faceId )
{
  int i;

  for (i=0;i<layer->nFrontParent;i++) 
    if (faceId == layer->frontParent[i]) return TRUE;

  return FALSE;
}

Layer *layerFront(Layer *layer, int front, int *nodes )
{

  if (front < 0 || front >= layerNFront(layer)) return NULL;
  nodes[0] = layer->front[front].globalNode[0];
  nodes[1] = layer->front[front].globalNode[1];
  nodes[2] = layer->front[front].globalNode[2];
  
  return layer;
}

Layer *layerFrontDirection(Layer *layer, int front, double *direction )
{
  int i, *nodes;
  double node0[3], node1[3], node2[3];
  double edge1[3], edge2[3], norm[3], length; 
  
  if (front < 0 || front >= layerNFront(layer) ) return NULL;
  nodes = layer->front[front].globalNode;

  if (layer->grid != gridNodeXYZ( layer->grid, nodes[0], node0 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[1], node1 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[2], node2 )) return NULL;

  for (i = 0 ; i < 3 ; i++ ){
    edge1[i] = node1[i] - node0[i];
    edge2[i] = node2[i] - node0[i];
  }

  norm[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1]; 
  norm[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2]; 
  norm[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0]; 
  length = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);

  if (length > 0.0) {
    for ( i=0;i<3;i++) direction[i] = norm[i]/length;
  }else{
    for ( i=0;i<3;i++) direction[i] = 0.0;
  }

  return layer;
}

Layer *layerMakeNormal(Layer *layer)
{
  int i, front, normal, globalNode;
  double direction[3], *norm, length;

  if (layerNFront(layer)==0) return NULL;
  layer->globalNode2Normal = malloc(layerMaxNode(layer)*sizeof(int));
  for (i=0;i<layerMaxNode(layer);i++) layer->globalNode2Normal[i]=EMPTY;
  normal = 0;
  for (front=0;front<layerNFront(layer);front++){
    for(i=0;i<3;i++){
      globalNode = layer->front[front].globalNode[i];
      if (EMPTY == layer->globalNode2Normal[globalNode] ){
	layer->globalNode2Normal[globalNode]=normal;
	layer->front[front].normal[i]=normal;
	normal++;
      }else{
	layer->front[front].normal[i]=layer->globalNode2Normal[globalNode];
      }
    }
  }
  layer->nnormal=normal;

  layer->normal = malloc( layer->nnormal * sizeof(Normal));
  for(normal=0;normal<layer->nnormal;normal++){ 
    layer->normal[normal].constrained = 0;
  }

  for(i=0;i<layerMaxNode(layer);i++){
    normal = layer->globalNode2Normal[i];
    if (normal!=EMPTY) {
      layer->normal[normal].root = i;
      layer->normal[normal].tip = EMPTY;
      layer->normal[normal].direction[0] = 0.0;
      layer->normal[normal].direction[1] = 0.0;
      layer->normal[normal].direction[2] = 0.0;
      layer->normal[normal].terminated = FALSE;
    }
  }

  layer->adj = adjCreate( layer->nnormal,layerNFront(layer)*3  );
  for (front=0;front<layerNFront(layer);front++){
    for(i=0;i<3;i++){
      normal = layer->front[front].normal[i];
      adjRegister( layer->adj, normal, front );
      layerFrontDirection(layer,front,direction);
      layer->normal[normal].direction[0] += direction[0];
      layer->normal[normal].direction[1] += direction[1];
      layer->normal[normal].direction[2] += direction[2];
    }
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    norm = layer->normal[normal].direction;
    length = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    if (length > 0.0) {
      for ( i=0;i<3;i++) norm[i] = norm[i]/length;
    }else{
      for ( i=0;i<3;i++) norm[i] = 0.0;
    }
  }

  return layer;
}

Layer *layerFrontNormals(Layer *layer, int front, int *normals )
{
  if (layerNNormal(layer) == 0 ) return NULL;
  if (front < 0 || front >= layerNFront(layer)) return NULL;
  normals[0] = layer->front[front].normal[0];
  normals[1] = layer->front[front].normal[1];
  normals[2] = layer->front[front].normal[2];
  
  return layer;
}

int layerNormalRoot(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return 0;
  return layer->normal[normal].root;
}

int layerNormalDeg(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return 0;
  return adjDegree(layer->adj, normal);
}

Layer *layerNormalFronts(Layer *layer, int normal, int nfront, int *fronts )
{
  int i;
  AdjIterator it;
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;
  i=0;
  for ( it = adjFirst(layer->adj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
    if (i>=nfront) return NULL;
    fronts[i] = adjItem(it);
    i++;
  }
  return layer;
}

Layer *layerNormalDirection(Layer *layer, int normal, double *direction )
{
  int i;
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  for ( i=0;i<3;i++) direction[i] = layer->normal[normal].direction[i];

  return layer;
}

Layer *layerVisibleNormals(Layer *layer)
{
  int normal, iter, front, i;
  double *dir, norm[3], mindir[3], dot, mindot, radian, length; 
  AdjIterator it;
  int minFront, lastFront;

  if (layerNNormal(layer) == 0 ) return NULL;

  for (normal=0;normal<layerNNormal(layer);normal++){
    lastFront = EMPTY;
    radian = 0.01;
    for (iter=0;iter<1000 && radian > 1.0e-15;iter++){
      dir = layer->normal[normal].direction;
      mindot = 2.0;
      mindir[0]=dir[0];
      mindir[1]=dir[1];
      mindir[2]=dir[2];
      minFront = EMPTY;
      for ( it = adjFirst(layer->adj,normal); 
	    adjValid(it); 
	    it = adjNext(it) ){
	front = adjItem(it);
	layerFrontDirection(layer,front,norm);
	dot = norm[0]*dir[0] + norm[1]*dir[1] + norm[2]*dir[2];
	if (dot<mindot) {
	  mindot = dot;
	  mindir[0]=norm[0];
	  mindir[1]=norm[1];
	  mindir[2]=norm[2];
	  minFront = front;
	}
      }
      if (minFront != lastFront) {
	radian = radian * 0.5;
	lastFront = minFront;
	//printf("normal %d, dot %f rad %e front %d\n",normal,mindot,radian,minFront);
      }
      dir[0] += radian*mindir[0];
      dir[1] += radian*mindir[1];
      dir[2] += radian*mindir[2];
      length = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
      if (length > 0.0) {
	for ( i=0;i<3;i++) dir[i] = dir[i]/length;
      }else{
	for ( i=0;i<3;i++) dir[i] = 0.0;
      }
    }
    if (mindot <= 0.0 ) 
      printf("ERROR: %s, %d, Invisible normal %d, min dot product %f\n",
	     __FILE__, __LINE__, normal, mindot);
  }

  return layer;
}

Layer *layerConstrainNormal(Layer *layer, int edgeface )
{
  int edgeId, faceId, face, nodes[3], id, i, normal;
  int n0, n1, normal0, normal1;
  int nCurveNode, *curve;
  int exist;

  if (layerNNormal(layer) == 0 ) return NULL;
  if (edgeface == 0 ) return NULL;
  
  if (edgeface > 0) {
    faceId = edgeface;
    for(face=0;face<gridMaxFace(layer->grid);face++){
      if (layer->grid == gridFace(layer->grid,face,nodes,&id) && id==faceId ){
	for(i=0;i<3;i++){
	  normal = layer->globalNode2Normal[nodes[i]];
	  if (normal != EMPTY) {
	    if ( layer->normal[normal].constrained >= 0) {
	      layer->normal[normal].constrained=faceId;
	    }
	  }
	}
	for(n0=0;n0<3;n0++){
	  n1 = n0 + 1; if ( n1 > 2 ) n1 = 0;
	  normal0 = layer->globalNode2Normal[nodes[n0]];
	  normal1 = layer->globalNode2Normal[nodes[n1]];
	  layerConstrainFrontSide( layer, normal0, normal1, faceId );
	}
      }
    }
  }else{
    edgeId = -edgeface;
    nCurveNode = gridGeomEdgeSize( layer->grid, edgeId );
    curve = malloc( nCurveNode * sizeof(int) );
    gridGeomEdge( layer->grid, edgeId, curve );
    for ( i=0; i<nCurveNode; i++){ 
      normal = layer->globalNode2Normal[curve[i]];
      if (normal != EMPTY) layer->normal[normal].constrained=-edgeId;
    }
    free(curve);
  }
  
  if ( !layerConstrainingGeometry(layer, edgeface) ) {
    if ( layer->constrainingGeometry == NULL ) {
      layer->nConstrainingGeometry = 1;
      layer->constrainingGeometry = 
	malloc( layer->nConstrainingGeometry * sizeof(int) );
      layer->constrainingGeometry[0] = edgeface;
    }else{
      layer->nConstrainingGeometry++;
      layer->constrainingGeometry = 
	realloc( layer->constrainingGeometry, 
		 layer->nConstrainingGeometry * sizeof(int) );
      layer->constrainingGeometry[layer->nConstrainingGeometry-1] = edgeface;
    }
  }

  return layer;
}

bool layerConstrainingGeometry(Layer *layer, int edgeface )
{
  int i;

  if ( edgeface == 0 ) return FALSE;

  for (i=0;i<layer->nConstrainingGeometry;i++)
    if (layer->constrainingGeometry[i]==edgeface) return TRUE;
  
  return FALSE;
}

int layerConstrained(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return 0;
  return layer->normal[normal].constrained;
}

Layer *layerConstrainFrontSide(Layer *layer, int normal0, int normal1, int bc )
{
  AdjIterator it;
  int i0, i1, front, side;
 
  if (normal0 < 0 || normal0 >= layerNNormal(layer) ) return NULL;
  if (normal1 < 0 || normal1 >= layerNNormal(layer) ) return NULL;

  for ( it = adjFirst(layer->adj,normal0); 
	adjValid(it); 
	it = adjNext(it) ){
    front = adjItem(it);
    for (i1=0;i1<3;i1++) {
      if (normal1 == layer->front[front].normal[i1]) {
	for (i0=0;i0<3;i0++) {
	  if (normal0 == layer->front[front].normal[i0]) {
	    side = MIN(i0,i1);
	    if ( side == 0 && 2 == MAX(i0,i1) ) side = 2;
	    layer->front[front].constrainedSide[side]=bc;
	  }
	}
      }   
    }
  }
  return layer;
}

int layerConstrainedSide(Layer *layer, int front, int side )
{
  if (front < 0 || front >= layerNFront(layer) ) return 0;

  if (side < 0 || side > 2 ) return 0;

  return layer->front[front].constrainedSide[side];
}

int layerNConstrainedSides(Layer *layer, int faceId )
{
  int front, i, nside;

  if (faceId==0) return 0;
  nside = 0;
  for (front=0;front<layerNFront(layer);front++){
    for(i=0;i<3;i++){
      if (layerConstrainedSide(layer,front,i)==faceId) nside++;
    }
  }
  return nside;
}

Layer *layerFindParentEdges(Layer *layer)
{
  int front, side, n0, n1, edgeId;

  if (layerNFront(layer) == 0 ) return NULL;
  if (layerNNormal(layer) == 0 ) return NULL;

  for (front=0;front<layerNFront(layer);front++){
    for(side=0;side<3;side++){
      n0 = side;
      n1 = side+1; if (n1>2) n1 = 0;
      n0 = layer->front[front].globalNode[n0];
      n1 = layer->front[front].globalNode[n1];
      edgeId = gridEdgeId(layer->grid,n0,n1);
      if (EMPTY != edgeId) layer->front[front].parentEdge[side]=edgeId;
    }
  }
  return layer;
}

Layer *layerSetParentEdge(Layer *layer, int normal0, int normal1, int edgeId )
{
  AdjIterator it;
  int i0, i1, front, side;
 
  if (normal0 < 0 || normal0 >= layerNNormal(layer) ) return NULL;
  if (normal1 < 0 || normal1 >= layerNNormal(layer) ) return NULL;

  for ( it = adjFirst(layer->adj,normal0); 
	adjValid(it); 
	it = adjNext(it) ){
    front = adjItem(it);
    for (i1=0;i1<3;i1++) {
      if (normal1 == layer->front[front].normal[i1]) {
	for (i0=0;i0<3;i0++) {
	  if (normal0 == layer->front[front].normal[i0]) {
	    side = MIN(i0,i1);
	    if ( side == 0 && 2 == MAX(i0,i1) ) side = 2;
	    layer->front[front].parentEdge[side]=edgeId;
	  }
	}
      }   
    }
  }
  return layer;
}

int layerParentEdge(Layer *layer, int front, int side )
{
  if (front < 0 || front >= layerNFront(layer) ) return 0;

  if (side < 0 || side > 2 ) return 0;

  return layer->front[front].parentEdge[side];
}

int layerNParentEdgeSegments(Layer *layer, int edgeId )
{
  int front, i, nSegments;

  if (edgeId==0) return 0;
  nSegments = 0;
  for (front=0;front<layerNFront(layer);front++){
    for(i=0;i<3;i++){
      if (layerParentEdge(layer,front,i)==edgeId) nSegments++;
    }
  }
  return nSegments;
}

Layer *layerTerminateNormal(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  layer->normal[normal].terminated = TRUE;

  return layer;
}

bool layerNormalTerminated(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return FALSE;

  return layer->normal[normal].terminated;
}

int layerNActiveNormal(Layer *layer )
{
  int normal, nActive;

  nActive=0;

  for ( normal=0 ; normal < layerNNormal(layer); normal++ ) 
    if ( !layer->normal[normal].terminated ) nActive++;

  return nActive;
}

Layer *layerAdvance(Layer *layer, double height )
{
  Grid *grid = layer->grid;
  int normal, normal0, normal1, root, tip, faceId, edgeId, i;
  int cell, node;
  int front, normals[3], n[6], side[2];
  double xyz[3];

  if (layerNNormal(layer) == 0 ) return NULL;

  for (normal=0;normal<layerNNormal(layer);normal++){
    root = layer->normal[normal].root;
    gridFreezeNode( grid, root );
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    if (layer->normal[normal].terminated){
      layer->normal[normal].tip = layer->normal[normal].root;
    }else{
      root = layer->normal[normal].root;
      gridNodeXYZ(grid,root,xyz);
      for(i=0;i<3;i++)xyz[i]=xyz[i]+height*layer->normal[normal].direction[i];
      tip = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);
      if ( EMPTY == tip) return NULL;
      layer->normal[normal].tip = tip;
      gridReconnectCellUnlessFrozen(grid, root, tip);
      faceId = layerConstrained(layer,normal);
      if (0 > faceId) {
	edgeId = -faceId;
	gridReconnectEdgeUnlessFrozen(grid, edgeId, root, tip);
      }
      gridCopySpacing(grid, root, tip );
      gridFreezeNode( grid, tip );
    }
  }

  /* reconnect faces for constrained frontside */
  for (front=0;front<layerNFront(layer);front++){
    for (i=0;i<3;i++){
      faceId = layerConstrainedSide(layer, front, i);
      if (faceId > 0) {
	normal0 = i;
	normal1 = i+1; if (normal1>2) normal1 = 0;
	normal0 = layer->front[front].normal[normal0];
	normal1 = layer->front[front].normal[normal1];
	gridReconnectFaceUnlessFrozen(grid, faceId, 
				      layer->normal[normal0].root, 
				      layer->normal[normal0].tip);
	gridReconnectFaceUnlessFrozen(grid, faceId, 
				      layer->normal[normal1].root, 
				      layer->normal[normal1].tip);
      }
    }    
  }

  // advance edges
  for (normal=0;normal<layerNNormal(layer);normal++) {
    edgeId = -layerConstrained(layer,normal);
    if (edgeId > 0) {
      root = layer->normal[normal].root;
      tip  = layer->normal[normal].tip;
      /* note that tip has been set to root on terminated normals */
      if (root != tip) gridAddEdge(grid,root,tip,edgeId,DBL_MAX,DBL_MAX);
    }

  }

  for (front=0;front<layerNFront(layer);front++){
    layerFrontNormals(layer, front, normals);
    for (i=0;i<3;i++) 
      layer->front[front].globalNode[i] = layer->normal[normals[i]].tip;

    /* note that tip has been set to root on terminated normals */
    /* the if (n[0]!=n[3]) checks are for layer termiantion */

    // advance faces
    for (i=0;i<3;i++){
      faceId = layerConstrainedSide(layer, front, i);
      if (faceId > 0) {
	side[1] = i;
	side[0] = i+1; if (side[0]>2) side[0] = 0;
	side[0] = normals[side[0]];
	side[1] = normals[side[1]];
	n[0] = layer->normal[side[0]].root;
	n[1] = layer->normal[side[1]].root;
	n[2] = layer->normal[side[0]].tip;
	n[3] = layer->normal[side[1]].tip;
	if (side[0]<side[1]){
	  if (n[1]!=n[3]) gridAddFace(grid,n[0],n[1],n[3],faceId);
	  if (n[0]!=n[2]) gridAddFace(grid,n[0],n[3],n[2],faceId);
	}else{
	  if (n[0]!=n[2]) gridAddFace(grid,n[0],n[1],n[2],faceId);
	  if (n[1]!=n[3]) gridAddFace(grid,n[2],n[1],n[3],faceId);
	}
      }
    }

    // advance cells
    /* pg. 82-85 of Garimella Thesis*/
    /* sort so that normals[0] is the smallest normal id*/
    if (normals[1]<normals[0] && normals[1]<normals[2]){
      normal = normals[1];
      normals[1] = normals[2];
      normals[2] = normals[0];
      normals[0] = normal;
    }
    if (normals[2]<normals[0] && normals[2]<normals[1]){
      normal = normals[2];
      normals[2] = normals[1];
      normals[1] = normals[0];
      normals[0] = normal;
    }

    for (i=0;i<3;i++){
      n[i]   = layer->normal[normals[i]].root;
      n[i+3] = layer->normal[normals[i]].tip;
    }

    if (normals[2]<normals[1]){
      if (n[0]!=n[3]) gridAddCell(grid, n[0], n[4], n[5], n[3]);
      if (n[2]!=n[5]) gridAddCell(grid, n[2], n[0], n[4], n[5]);
      if (n[1]!=n[4]) gridAddCell(grid, n[2], n[0], n[1], n[4]);
    }else{
      if (n[0]!=n[3]) gridAddCell(grid, n[0], n[4], n[5], n[3]);
      if (n[1]!=n[4]) gridAddCell(grid, n[0], n[1], n[5], n[4]);
      if (n[2]!=n[5]) gridAddCell(grid, n[2], n[0], n[1], n[5]);
    }
    
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    faceId = layerConstrained(layer,normal);
    if (0 > faceId) {
      edgeId = -faceId;
      gridProjectNodeToEdge(grid, layer->normal[normal].root, edgeId );
      // need to do faces to get the UV vals but have xyz.
    }
    layer->normal[normal].root = layer->normal[normal].tip;
    layer->normal[normal].tip = EMPTY;
    if (0 < faceId) {
      gridProjectNodeToFace(grid, layer->normal[normal].root, faceId );
    }
    if (0 > faceId) {
      edgeId = -faceId;
      gridProjectNodeToEdge(grid, layer->normal[normal].root, edgeId );
      // need to do faces to get the UV vals but have xyz.
    }
    gridFreezeNode(grid,layer->normal[normal].root);
  }

  return layer;
}

Layer *layerWiggle(Layer *layer, double height )
{
  Grid *grid = layer->grid;
  int normal, root, faceId, edgeId, i;
  double xyz[3];

  if (layerNNormal(layer) == 0 ) return NULL;

  for (normal=0;normal<layerNNormal(layer);normal++){
    if ( !layerNormalTerminated(layer,normal) ) {
      root = layer->normal[normal].root;
      gridNodeXYZ(grid,root,xyz);
      for(i=0;i<3;i++)xyz[i]=xyz[i]+height*layer->normal[normal].direction[i];
      gridSetNodeXYZ(grid, root, xyz);
      faceId = layerConstrained(layer,normal);
      if (0 < faceId) {
	gridProjectNodeToFace(grid, layer->normal[normal].root, faceId );
      }
      if (0 > faceId) {
	edgeId = -faceId;
	gridProjectNodeToEdge(grid, layer->normal[normal].root, edgeId );
      }
    }
  }

  return layer;
}

Layer *layerSmoothLayerNeighbors(Layer *layer )
{
  int normal;

  for (normal=0;normal<layerNNormal(layer);normal++)
    gridSmoothNearNode(layer->grid,layer->normal[normal].root);

  return layer;
}


Layer *layerTerminateNormalWithSpacing(Layer *layer, double spacing)
{
  int normal, nterm;;
  if (layerNNormal(layer) == 0 ) return NULL;

  nterm = 0;
  for (normal=0;normal<layerNNormal(layer);normal++){
    if (gridSpacing(layer->grid, layer->normal[normal].root ) < spacing ) {
      nterm++;
      layerTerminateNormal(layer, normal);
    }
  }
  printf("normals %d of %d terminated\n",nterm,layerNNormal(layer) );
  return layer;
}

Layer *layerTerminateNormalWithX(Layer *layer, double x)
{
  int normal, nterm;
  double xyz[3];
  
  if (layerNNormal(layer) == 0 ) return NULL;

  nterm = 0;
  for (normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(layer->grid, layer->normal[normal].root, xyz);
    if (x > 0.0 ) {
      if (xyz[0]>x) { layerTerminateNormal(layer, normal); nterm++; }
    }else{
      if (xyz[0]<x) { layerTerminateNormal(layer, normal); nterm++; }     
    }
  }
  printf("normals %d of %d terminated\n",nterm,layerNNormal(layer) );
  return layer;
}

Layer *layerInsertPhantomFront(Layer *layer, double dz )
{
  Grid *grid = layer->grid;
  int normal, faceId, edgeId, newnode;
  double xyz[3], spacing, m;

  layerVisibleNormals(layer);

  for (normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(grid,layer->normal[normal].root,xyz);
    spacing = gridSpacing(grid, layer->normal[normal].root );
    layer->normal[normal].root = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);
    m = 1.0/spacing/spacing;
    gridSetMap(grid,layer->normal[normal].root,m,0.0,0.0,m,0.0,m); 
    layer->normal[normal].tip = EMPTY;
  }

  layerWiggle(layer, dz );

  for (normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(grid,layer->normal[normal].root,xyz);
    spacing = gridSpacing(grid, layer->normal[normal].root );
    if (0 == layerConstrained(layer,normal)){
      newnode = gridInsertInToVolume(grid, xyz[0], xyz[1], xyz[2]);
    }else{
      faceId = layerConstrained(layer,normal);
      if (0 < faceId) {
	gridForceNodeToFace(grid, layer->normal[normal].root, faceId );
      }
      if (0 > faceId) {
	edgeId = -faceId;
	gridForceNodeToEdge(grid, layer->normal[normal].root, edgeId );
      }
      newnode = gridInsertInToGeomFace(grid, xyz[0], xyz[1], xyz[2]);
    }
    if (EMPTY == newnode) 
      printf("Could not insert node norm %8d x %10.5f y %10.5f z %10.5f\n",
	     normal, xyz[0], xyz[1], xyz[2]);
    layer->normal[normal].tip = newnode;
    gridFreezeNode( grid, newnode );
    m = 1.0/spacing/spacing;
    gridSetMap(grid,newnode,m,0.0,0.0,m,0.0,m); 
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    gridRemoveNode(grid,layer->normal[normal].root);
    layer->normal[normal].root = EMPTY;
  }

  return layer;
}

Layer *layerVerifyPhantomEdges(Layer *layer)
{
  Grid *grid = layer->grid;
  int front, normals[3], nline, ngot, n0, n1;

  nline = 0;
  ngot  = 0;
  for(front=0;front<layerNFront(layer);front++){
    layerFrontNormals(layer, front, normals );
    n0 = layer->normal[normals[0]].tip;
    n1 = layer->normal[normals[1]].tip;
    if ( n0 != EMPTY && n1 != EMPTY ) {
      nline++;
      if (grid == gridVerifyEdgeExists(grid, n0, n1) ){
	ngot++;
      }else{
	printf("line %d<->%d is missing.\n",n0,n1);
      }
    }
    n0 = layer->normal[normals[1]].tip;
    n1 = layer->normal[normals[2]].tip;
    if ( n0 != EMPTY && n1 != EMPTY ) {
      nline++;
      if (grid == gridVerifyEdgeExists(grid, n0, n1) ){
	ngot++;
      }else{
	printf("line %d<->%d is missing.\n",n0,n1);
      }
    }
    n0 = layer->normal[normals[2]].tip;
    n1 = layer->normal[normals[0]].tip;
    if ( n0 != EMPTY && n1 != EMPTY ) {
      nline++;
      if (grid == gridVerifyEdgeExists(grid, n0, n1) ){
	ngot++;
      }else{
	printf("line %d<->%d is missing.\n",n0,n1);
      }
    }
  }
  printf("lines %d of %d verified.\n",ngot,nline);
 
  return layer;
}

Layer *layerVerifyPhantomFaces(Layer *layer)
{
  Grid *grid = layer->grid;
  int front, normals[3], nface, ngot, n0, n1, n2;

  nface = 0;
  ngot  = 0;
  for(front=0;front<layerNFront(layer);front++){
    layerFrontNormals(layer, front, normals );
    n0 = layer->normal[normals[0]].tip;
    n1 = layer->normal[normals[1]].tip;
    n2 = layer->normal[normals[2]].tip;
    if ( n0 != EMPTY && n1 != EMPTY && n2 != EMPTY ) {
      nface++;
      if (grid == gridVerifyFaceExists(grid, n0, n1, n2) ){
	ngot++;
      }else{
	printf("face %d %d %d is missing.\n",n0,n1,n2);
      }
    }
  }
  printf("faces %d of %d verified.\n",ngot,nface);
 
  return layer;
}
