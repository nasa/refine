
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef MESHERX_H
#define MESHERX_H

#include "master_header.h"

#include "grid.h"
#include "layer.h"

BEGIN_C_DECLORATION

int MesherX_DiscretizeVolume( int maxNodes, double scale, char *project, 
			      bool mixedElement );

Layer *layerFormAdvancingLayerWithCADGeomBCS( int volumeId, Grid * );

int layerTerminateNormalWithBGSpacing(Layer *, double ratio);

Layer *layerRebuildEdges(Layer *, int vol);
Layer *layerRebuildFaces(Layer *, int vol);
Layer *layerRebuildVolume(Layer *, int vol);

END_C_DECLORATION

#endif /* MESHERX_H */
