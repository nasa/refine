
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef MESHERXREBUILD_H
#define MESHERXREBUILD_H

#include "master_header.h"

#include "layer.h"

BEGIN_C_DECLORATION

Layer *layerRebuildInterior(Layer *, int vol);

Layer *layerRebuildEdges(Layer *, int vol);
Layer *layerRebuildFaces(Layer *, int vol);
Layer *layerRebuildVolume(Layer *, int vol);

END_C_DECLORATION

#endif /* MESHERXREBUILD_H */
