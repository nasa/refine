
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef LAYER_H
#define LAYER_H

#include "master_header.h"
#include "grid.h"

BEGIN_C_DECLORATION

typedef struct Layer Layer;

Layer *layerCreate(Grid *);
void layerFree(Layer *);
int layerNFront(Layer *);
int layerMaxNode(Layer *);
Layer *layerMakeFront(Layer *, int nbc, int *bc);
Layer *layerFront(Layer *, int front, int *nodes);

END_C_DECLORATION

#endif /* LAYER_H */
