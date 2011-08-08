

/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */



#ifndef MESHERXINIT_H
#define MESHERXINIT_H

#include "refine_defs.h"
#include "grid.h"
#include "layer.h"

BEGIN_C_DECLORATION

Layer *mesherxInit(int vol, int maxNodes);
Layer *layerFormAdvancingLayerWithCADGeomBCS( int vol, Grid *grid );

END_C_DECLORATION

#endif /* MESHERXINIT_H */
