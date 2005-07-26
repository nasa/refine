
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDINSERT_H
#define GRIDINSERT_H

#include "refine_defs.h"
#include "queue.h"
#include "grid.h"

BEGIN_C_DECLORATION
Grid *gridThrash(Grid *g );
Grid *gridRemoveAllNodes(Grid *g );

Grid *gridAdapt(Grid *g, double minLength, double maxLength );
Grid *gridAdaptBasedOnConnRankings(Grid *g );
Grid *gridAdaptLongShort(Grid *g, double minLength, double maxLength,
			 GridBool debug_split );

int gridSplitEdge(Grid *g, int n0, int n1 );
int gridSplitEdgeRatio(Grid *g, Queue *q, int n0, int n1, double ratio);
int gridSplitEdgeRepeat(Grid *g, Queue *q, int n0, int n1 );
int gridSplitEdgeForce(Grid *g, Queue *q, int n0, int n1 );
int gridSplitEdgeIfNear(Grid *g, int n0, int n1, 
			double newX, double newY, double newZ);
int gridSplitFaceAt(Grid *g, int face,  
		    double newX, double newY, double newZ);
int gridSplitCellAt(Grid *g, int cell,  
		    double newX, double newY, double newZ);
int gridInsertInToGeomEdge(Grid *g, double newX, double newY, double newZ);
int gridInsertInToGeomFace(Grid *g, double newX, double newY, double newZ);
int gridInsertInToVolume(Grid *g, double newX, double newY, double newZ);

Grid *gridCollapseEdgeToAnything(Grid *g, Queue *q, int n0, int n1);
Grid *gridCollapseEdge(Grid *g, Queue *q, int n0, int n1,
		       double requestedRatio );

Grid *gridFreezeGoodNodes(Grid *g, double goodAR, 
			  double minLength, double maxLength );

Grid *gridVerifyEdgeExists(Grid *g, int n0, int n1);
Grid *gridVerifyFaceExists(Grid *g, int n0, int n1, int n2);

END_C_DECLORATION

#endif /* GRIDINSERT_H */
