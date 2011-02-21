
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */

#ifndef GRIDINSERT_H
#define GRIDINSERT_H

#include "refine_defs.h"
#include "queue.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridAdapt(Grid *g, double minLength, double maxLength );

int gridSplitEdge(Grid *g, int n0, int n1 );
int gridSplitEdgeRatio(Grid *g, Queue *q, int n0, int n1, double ratio);

int gridSplitFaceAt(Grid *g, int *face_nodes, double *xyz);
int gridSplitCellAt(Grid *g, int cell, double *xyz);

Grid *gridCollapseEdgeToAnything(Grid *g, Queue *q, int n0, int n1);
/* node1 is removed */
Grid *gridCollapseEdge(Grid *g, Queue *q, int n0, int n1,
		       double requestedRatio );

END_C_DECLORATION

#endif /* GRIDINSERT_H */
