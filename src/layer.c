
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
#include "grid.h"
#include "adj.h"

typedef struct Normal Normal;
struct Normal {
  int constrained;
  int faceId1, faceId2;
  int root, tip;
  double direction[3];
  bool terminated;
};

typedef struct Front Front;
struct Front {
  int globalNode[3];
  int normal[3];
};

struct Layer {
  Grid *grid;
  int nfront;
  Front *front;
  int nnormal;
  Normal *normal;
  int *globalNode2Normal;
  Adj *adj;
};

Layer *layerCreate( Grid *grid )
{
  Layer *layer;
  layer = malloc(sizeof(Layer));
  layer->grid = grid;
  layer->nfront=0;
  layer->front=NULL;
  layer->nnormal=0;
  layer->normal=NULL;
  layer->globalNode2Normal=NULL;
  layer->adj=NULL;
  return layer;
}

void layerFree(Layer *layer)
{
  if (layer->adj != NULL) adjFree(layer->adj);
  if (layer->globalNode2Normal != NULL) free(layer->globalNode2Normal);
  if (layer->normal != NULL) free(layer->normal);
  if (layer->front != NULL) free(layer->front);
  free(layer);
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
  int ibc, face, nface, id, nodes[3], ifront;

  layer->nfront=0;

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
	ifront++;
      }
    }
  }
  
  return layer;
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
    layer->normal[normal].faceId1 = EMPTY;
    layer->normal[normal].faceId2 = EMPTY;
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

Layer *layerConstrainNormal(Layer *layer, int bc )
{
  int face, nodes[3], id, i, normal;
  int edge, nCurveNode, *curve;
  int exist;
  if (layerNNormal(layer) == 0 ) return NULL;
  
  if (bc > 0) {
    for(face=0;face<gridMaxFace(layer->grid);face++){
      if (layer->grid == gridFace(layer->grid,face,nodes,&id) &&
	  id==bc ) {
	for(i=0;i<3;i++){
	  normal = layer->globalNode2Normal[nodes[i]];
	  if (normal != EMPTY) {
	    if ( layer->normal[normal].constrained >= 0) {
	      layer->normal[normal].constrained=bc;
	    }else{
	      if (layer->normal[normal].faceId1 == EMPTY ||
		  layer->normal[normal].faceId1 == bc ) {
		layer->normal[normal].faceId1 = bc;
	      }else{
		layer->normal[normal].faceId2 = bc;
	      }
	    }
	  }
	}
      }
    }
  }else{
    edge = -bc;
    nCurveNode = gridGeomEdgeSize( layer->grid, edge );
    curve = malloc( nCurveNode * sizeof(int) );
    gridGeomEdge( layer->grid, edge, curve );
    for ( i=0; i<nCurveNode; i++){ 
      normal = layer->globalNode2Normal[curve[i]];
      if (normal != EMPTY) layer->normal[normal].constrained=bc;
    }
    free(curve);
  }
  
  return layer;
}

int layerConstrained(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return 0;
  return layer->normal[normal].constrained;
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
  int normal, root, tip, faceId, edgeId, i;
  int cell, node;
  int front, normals[3], n[6], side[2];
  double xyz[3];
  AdjIterator it;  

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
      if (0 < faceId) {
	gridReconnectFaceUnlessFrozen(grid, faceId, root, tip);
      }
      if (0 > faceId) {
	edgeId = -faceId;
	gridReconnectEdgeUnlessFrozen(grid, edgeId, root, tip);
	gridReconnectFaceUnlessFrozen(grid, layer->normal[normal].faceId1, root, tip);
	gridReconnectFaceUnlessFrozen(grid, layer->normal[normal].faceId2, root, tip);
      }
      gridCopySpacing(grid, root, tip );
      gridFreezeNode( grid, tip );
    }
  }

  for (front=0;front<layerNFront(layer);front++){
    for (i=0;i<3;i++) layer->front[front].globalNode[i] = 
			layer->normal[layer->front[front].normal[i]].root;
    layerFrontNormals(layer, front, normals);
    /* pg. 82-85 of Garimella Thesis*/
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

    /* note that tip has been set to root on terminated normals */
    /* the if (n[0]!=n[3]) checks are for layer termiantion */

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
    
    edgeId = -layerConstrained(layer,normals[0]);
    if (edgeId > 0 && n[0]!=n[3]) 
      gridAddEdge(grid,n[0],n[3],edgeId,DBL_MAX,DBL_MAX);

    edgeId = -layerConstrained(layer,normals[1]);
    if (edgeId > 0 && n[1]!=n[4]) 
      gridAddEdge(grid,n[1],n[5],edgeId,DBL_MAX,DBL_MAX);

    edgeId = -layerConstrained(layer,normals[2]);
    if (edgeId > 0 && n[2]!=n[4]) 
      gridAddEdge(grid,n[2],n[5],edgeId,DBL_MAX,DBL_MAX);

    if (0 != layerConstrained(layer,normals[0]) && 
	0 != layerConstrained(layer,normals[1]) ){
      side[0] = normals[1];
      side[1] = normals[0];

      faceId = MAX(layerConstrained(layer,side[0]),
		   layerConstrained(layer,side[1]));
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
    if (0 != layerConstrained(layer,normals[1]) && 
	0 != layerConstrained(layer,normals[2]) ){
      side[0] = normals[2];
      side[1] = normals[1];

      faceId = MAX(layerConstrained(layer,side[0]),
		   layerConstrained(layer,side[1]));
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
    if (0 != layerConstrained(layer,normals[2]) && 
	0 != layerConstrained(layer,normals[0]) ){
      side[0] = normals[0];
      side[1] = normals[2];

      faceId = MAX(layerConstrained(layer,side[0]),
		   layerConstrained(layer,side[1]));
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

  for (normal=0;normal<layerNNormal(layer);normal++){
    layer->normal[normal].root = layer->normal[normal].tip;
    layer->normal[normal].tip = EMPTY;
    faceId = layerConstrained(layer,normal);
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

Layer *layerInsertPhantomFront(Layer *layer)
{
  Grid *grid = layer->grid;
  int normal, faceId, edgeId, newnode;
  double xyz[3];

  layerVisibleNormals(layer);

  for (normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(grid,layer->normal[normal].root,xyz);
    layer->normal[normal].root = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);
  }

  layerWiggle(layer, 0.22 );

  for (normal=0;normal<layerNNormal(layer);normal++){
    if (0 != layerConstrained(layer,normal)){
      faceId = layerConstrained(layer,normal);
      if (0 < faceId) {
	gridForceNodeToFace(grid, layer->normal[normal].root, faceId );
      }
      if (0 > faceId) {
	edgeId = -faceId;
	gridForceNodeToEdge(grid, layer->normal[normal].root, edgeId );
      }
      gridNodeXYZ(grid,layer->normal[normal].root,xyz);
      newnode = gridInsertInToGeomFace(grid, xyz[0], xyz[1], xyz[2]);
      printf("insert node %10d x %10.5f y %10.5f z %10.5f\n",
	     newnode,xyz[0], xyz[1], xyz[2]);
      if (EMPTY == newnode) printf("Could not insert node %d\n",normal);
      gridFreezeNode( grid, newnode );
    }
  }
 
  for (normal=0;normal<layerNNormal(layer);normal++){
    gridRemoveNode(grid,layer->normal[normal].root);
  }

  return layer;
}
