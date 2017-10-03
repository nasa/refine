
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

#ifndef GRIDMPI_H

#include "refine_defs.h"
#include "queue.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridIdentityNodeGlobal(Grid *g, int offset );
Grid *gridSetAllLocal(Grid *g );
Grid *gridSetGhost(Grid *g, int node );

Grid *gridParallelAdapt(Grid *g, Queue *q, 
		       double minLength, double maxLength );
Grid *gridParallelPreProject(Grid *g, Queue *q );
int gridParallelEdgeSplit(Grid *g, Queue *q, int node0, int node1 );
Grid *gridParallelEdgeCollapse(Grid *g, Queue *q, int node0, int node1 );

Grid *gridParallelSmooth(Grid *grid, GridBool localOnly,
			 double optimizationLimit, double laplacianLimit,
                         GridBool smoothOnSurface );
Grid *gridParallelRelaxNegativeCells(Grid *grid, 
				     GridBool localOnly, 
				     GridBool smoothOnSurface );

Grid *gridParallelRelaxNegativeFaceAreaUV(Grid *grid, 
					  GridBool localOnly );

Grid *gridParallelSwap(Grid *grid, Queue *queue, double ARlimit );
Grid *gridParallelEdgeSwap(Grid *g, Queue *q, int node0, int node1 );
Grid *gridApplyQueue(Grid *g, Queue *q );

Grid *gridGhostDataCountByPartition(Grid *g, int total_number_of_partitions,
				    int *partition_data_count);

END_C_DECLORATION

#endif /* GRIDMPI_H */
