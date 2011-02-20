
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
Grid *gridAdapt2(Grid *g );
Grid *gridAdaptBasedOnConnRankings(Grid *g );
Grid *gridAdaptVolumeEdges(Grid *g);

int gridSplitEdge(Grid *g, int n0, int n1 );
int gridSplitEdgeRatio(Grid *g, Queue *q, int n0, int n1, double ratio);
int gridSplitEdgeRepeat(Grid *g, Queue *q, int n0, int n1, GridBool debug_split );
int gridSplitEdgeForce(Grid *g, Queue *q, int n0, int n1, GridBool debug_split );
int gridSplitEdgeIfNear(Grid *g, int n0, int n1, double *xyz);
int gridSplitFaceAt(Grid *g, int *face_nodes, double *xyz);
int gridSplitCellAt(Grid *g, int cell, double *xyz);
int gridInsertInToGeomEdge(Grid *g, double *xyz);
int gridInsertInToGeomFace(Grid *g, double *xyz);
int gridInsertInToVolume(Grid *g, double *xyz);

Grid *gridCollapseEdgeToAnything(Grid *g, Queue *q, int n0, int n1);
/* node1 is removed */
Grid *gridCollapseEdge(Grid *g, Queue *q, int n0, int n1,
		       double requestedRatio );

Grid *gridFreezeGoodNodes(Grid *g, double goodAR, 
			  double minLength, double maxLength );

Grid *gridVerifyEdgeExists(Grid *g, int n0, int n1);
Grid *gridVerifyFaceExists(Grid *g, int n0, int n1, int n2);

Grid *gridCollapseWedgeCells(Grid *g);

Grid *gridSplitSliverCell(Grid *g, Queue *q, int sliver_cell);

END_C_DECLORATION

#endif /* GRIDINSERT_H */
