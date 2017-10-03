
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
