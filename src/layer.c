
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
};

Layer *layerCreate(void)
{
  return malloc(sizeof(Layer));
}

void layerFree(Layer *layer)
{
  free(layer);
}

int layerNFront(Layer *layer)
{
  return layer->nfront;
}
