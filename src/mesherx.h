
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

#include "layer.h"

BEGIN_C_DECLORATION

int MesherX_DiscretizeVolume( int maxNodes, double scale, char *project, 
			      bool mixedElement,
			      bool blendElement,
			      bool qualityImprovement,
			      bool copyGridY,
			      bool bil );

Layer *layerComputeNormalRateWithBGSpacing(Layer *, double finalRatio);
Layer *layerComputeInitialCellHeightWithBGSpacing(Layer *, double finalRatio);

Layer *layerCreateWakeWithBGSpacing(Layer *, double *origin, double *direction, 
				    double length );

int layerTerminateNormalWithBGSpacing(Layer *, 
				      double normalRatio, double edgeRatio);

END_C_DECLORATION

#endif /* MESHERX_H */
