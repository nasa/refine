
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  


#ifndef ADJ_H
#define ADJ_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

typedef struct Adj Adj;
typedef struct NodeItem NodeItem;
typedef NodeItem * AdjIterator;

struct NodeItem {
  int item;
  NodeItem *next;
};

struct Adj {
  int nnode, nadj, chunkSize;
  NodeItem *node2item;
  NodeItem **first;
  NodeItem *current;
  NodeItem *blank;
};

Adj *adjCreate( int nnode, int nadj, int chunkSize );
void adjFree( Adj *adj );

int adjNNode( Adj *adj );
int adjNAdj( Adj *adj );
int adjChunkSize( Adj *adj );

Adj *adjRealloc( Adj *adj, int nnode );

Adj *adjRegister( Adj *adj, int node, int item );
Adj *adjRemove( Adj *adj, int node, int item );

#define adjValid(iterator) (iterator!=NULL)
#define adjMore(iterator) ((iterator!=NULL)&&(iterator->next != NULL))
#define adjFirst(adj,node) \
((node < 0 || node >= adj->nnode)?NULL:adj->first[node])

#define adjItem(iterator) (iterator==NULL?EMPTY:iterator->item)
#define adjNext(iterator) (iterator==NULL?NULL:iterator->next)

AdjIterator adjGetCurrent( Adj *adj );
Adj *adjSetCurrent( Adj *adj, AdjIterator iterator );

GridBool adjExists( Adj *adj, int node, int item );
int adjDegree( Adj *adj, int node );

END_C_DECLORATION

#endif /* ADJ_H */
