
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "layer.h"

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
  int *globalNode2Normal;
};

Layer *layerCreate( Grid *grid )
{
  Layer *layer;
  layer = malloc(sizeof(Layer));
  layer->grid = grid;
  layer->nfront=0;
  layer->front=NULL;
  layer->nnormal=0;
  layer->globalNode2Normal=NULL;
  return layer;
}

void layerFree(Layer *layer)
{
  if (layer->globalNode2Normal != NULL) free(layer->globalNode2Normal);
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
    printf("boundary %d with %d faces added to the %d face front.\n",
	   bc[ibc],nface,layer->nfront);
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

Layer *layerMakeNormal(Layer *layer)
{
  int i, front, normal, globalNode;
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
  return layer;
}

Layer *layerFrontNormal(Layer *layer, int front, int *normals )
{
  if (layerNNormal(layer) == 0 ) return NULL;
  if (front < 0 || front >= layerNFront(layer)) return NULL;
  normals[0] = layer->front[front].normal[0];
  normals[1] = layer->front[front].normal[1];
  normals[2] = layer->front[front].normal[2];
  
  return layer;
}

