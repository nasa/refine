
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDFILLER_H
#define GRIDFILLER_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

Grid *gridLoadPart( char *modeler, char *project, int maxnode );
Grid *gridFillFromPart( int vol, int maxnode );
int gridSavePart( Grid *grid, char *project );
int gridSavePartExplicitly( Grid *grid, char *project );

END_C_DECLORATION

#endif /* GRIDFILLER_H */
