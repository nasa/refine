
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
  Grid *grid;
};

Layer *layerCreate( Grid *grid )
{
  Layer *layer;
  layer = malloc(sizeof(Layer));
  layer->grid = grid;
  return layer;
}

void layerFree(Layer *layer)
{
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
  return layer;
}
