
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDINSERT_H
#define GRIDINSERT_H

#include "master_header.h"
#include "grid.h"

BEGIN_C_DECLORATION
Grid *gridThrash(Grid *g );
Grid *gridAdapt(Grid *g, double minLength, double maxLength );
int gridSplitEdge(Grid *g, int n0, int n1 );
Grid *gridCollapseEdge(Grid *g, int n0, int n1, double ratio );

Grid *gridFreezeGoodNodes(Grid *g, double goodAR, 
			  double minLength, double maxLength );

END_C_DECLORATION

#endif /* GRIDINSERT_H */
