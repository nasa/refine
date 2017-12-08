
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#ifndef MESHERX_H
#define MESHERX_H

#include "refine_defs.h"

#include "layer.h"

BEGIN_C_DECLORATION

int MesherX_DiscretizeVolume( int maxNodes, double scale, char *project, 
			      GridBool mixedElement,
			      GridBool blendElement,
			      GridBool qualityImprovement,
			      GridBool copyGridY,
			      GridBool bil );

Layer *layerComputeNormalRateWithBGSpacing(Layer *, double finalRatio);
Layer *layerComputeInitialCellHeightWithBGSpacing(Layer *, double finalRatio);

Layer *layerCreateWakeWithBGSpacing(Layer *, double *origin, double *direction, 
				    double length );

int layerTerminateNormalWithBGSpacing(Layer *, 
				      double normalRatio, double edgeRatio);

END_C_DECLORATION

#endif /* MESHERX_H */
