
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
typedef struct NodeItem NodeItem;
typedef NodeItem * AdjIterator;

Adj *adjCreate( int nnode, int perNode );
void adjFree( Adj *adj );

int adjNNode( Adj *adj );
Adj *adjRegister( Adj *adj, int node, int item );
Adj *adjRemove( Adj *adj, int node, int item );

bool adjValid( AdjIterator iterator );
bool adjMore( AdjIterator iterator );
AdjIterator adjFirst( Adj *adj, int node );
int adjItem( AdjIterator iterator );
AdjIterator adjNext( AdjIterator iterator );

AdjIterator adjGetCurrent( Adj *adj );
Adj *adjSetCurrent( Adj *adj, AdjIterator iterator );

bool adjExists( Adj *adj, int node, int item );
int adjDegree( Adj *adj, int node );

END_C_DECLORATION

#endif /* ADJ_H */
