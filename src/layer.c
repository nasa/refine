
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "layer.h"

struct Layer {
  int nfront;
  int *front;
  Grid *grid;
};

Layer *layerCreate( Grid *grid )
{
  Layer *layer;
  layer = malloc(sizeof(Layer));
  layer->nfront=0;
  layer->front=NULL;
  layer->grid = grid;
  return layer;
}

void layerFree(Layer *layer)
{
  if (layer->front != NULL) free(layer->front);
  free(layer);
}

int layerNFront(Layer *layer)
{
  return layer->nfront;
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
  
  layer->front = malloc( 3 * layer->nfront * sizeof(int) );

  ifront =0;
  for (ibc=0;ibc<nbc;ibc++){
    for(face=0;face<gridMaxFace(layer->grid);face++){
      if (layer->grid == gridFace(layer->grid,face,nodes,&id) &&
	  id==bc[ibc] ) {
	layer->front[0+3*ifront] = nodes[0];
	layer->front[1+3*ifront] = nodes[1];
	layer->front[2+3*ifront] = nodes[2];
	ifront++;
      }
    }
  }
  
  return layer;
}

Layer *layerFront(Layer *layer, int front, int *nodes )
{

  if (front < 0 || front >= layer->nfront) return NULL;
  nodes[0] = layer->front[0+3*front];
  nodes[1] = layer->front[1+3*front];
  nodes[2] = layer->front[2+3*front];
  
  return layer;
}
