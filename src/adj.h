
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef ADJ_H
#define ADJ_H

#include "master_header.h"

BEGIN_C_DECLORATION

typedef struct Adj Adj;

Adj *adjCreate(int nnode, int perNode);
void adjFree(Adj *adj);

int adjNNode(Adj *adj);

END_C_DECLORATION

#endif /* ADJ_H */
