
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef MESHERX_HLA_H
#define MESHERX_HLA_H

#include "master_header.h"

#include "layer.h"

BEGIN_C_DECLORATION

int MesherX_DiscretizeVolumeHLA( int maxNodes,
				 double scale,
				 double scalev,
				 char *project, 
                                char *outputName,
			      bool mixedElement,
			      bool blendElement,
			      bool qualityImprovement,
			      bool bil );

Layer *layerComputeNormalRateWithBGSpacing2(Layer *, double finalRatio);


END_C_DECLORATION

#endif /* MESHERX_HLA_H */
