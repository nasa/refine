
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "adj.h"

static void adjAllocateAndInitNode2Item(Adj *adj)
{
  int i;
  adj->node2item = (NodeItem *)malloc( adj->nadj * sizeof(NodeItem));

  for ( i=0 ; i < (adj->nadj-1) ; i++ ) {
    adj->node2item[i].item = EMPTY;
    adj->node2item[i].next = &adj->node2item[i+1];
  }
  adj->node2item[adj->nadj-1].item = EMPTY;
  adj->node2item[adj->nadj-1].next = NULL;

  adj->blank = adj->node2item;
}

Adj* adjCreate( int nnode, int nadj, int chunkSize )
{
  int node;
  Adj *adj;

  adj = malloc( sizeof(Adj) );

  adj->nnode     = MAX(nnode,1);
  adj->nadj      = MAX(nadj,1);
  adj->chunkSize = MAX(chunkSize,1);

  adjAllocateAndInitNode2Item(adj);

  adj->current = NULL;

  adj->first = malloc( adj->nnode * sizeof(NodeItem*) );

  for ( node=0 ; node<adj->nnode; node++ ) adj->first[node] = NULL; 

  return adj;
}

void adjFree( Adj *adj )
{
  free( adj->first );
  free( adj );
}

int adjNNode( Adj *adj )
{
  return adj->nnode;
}

int adjNAdj( Adj *adj )
{
  return adj->nadj;
}

int adjChunkSize( Adj *adj )
{
  return adj->chunkSize;
}

Adj *adjRealloc( Adj *adj, int nnode )
{
  NodeItem *remove;
  int node, oldSize, newSize;
  oldSize = adj->nnode;
  newSize = MAX(nnode,1);
  if ( oldSize > newSize) {
    for ( node=newSize ; node<oldSize; node++ ) {
      while ( NULL != adj->first[node] ) {
	remove = adj->first[node];
	adj->first[node] = remove->next;
	remove->item = EMPTY;
	remove->next = adj->blank;
	adj->blank = remove;
      }
    }
  }
  adj->nnode = newSize;
  adj->first = realloc( adj->first, adj->nnode * sizeof(NodeItem*) );
  for ( node=oldSize ; node<adj->nnode; node++ ) adj->first[node] = NULL; 
  return adj;
}

Adj *adjRegister( Adj *adj, int node, int item )
{
  NodeItem *oldnode2item;
  NodeItem *oldfirst;
  NodeItem *new;
  int copynode;

  if (node>=adj->nnode || node<0) return NULL;

  if (NULL == adj->blank) {
    adj->nadj += adj->chunkSize;
    oldnode2item = adj->node2item;
    adjAllocateAndInitNode2Item(adj);
    for ( copynode = 0 ; copynode < adj->nnode ; copynode++ ) {
      oldfirst = adj->first[copynode];
      adj->first[copynode] = NULL;
      while ( NULL != oldfirst ) {
	new = adj->blank;
	adj->blank = adj->blank->next;
	new->next = adj->first[copynode];
	adj->first[copynode] = new;
	new->item = oldfirst->item;
	oldfirst = oldfirst->next;
      }
    }
    free(oldnode2item);
  }

  adj->current = NULL;

  new = adj->blank;
  adj->blank = adj->blank->next;
  new->next = adj->first[node];
  adj->first[node] = new;
  new->item = item;
  return adj;
}

Adj* adjRemove(Adj *adj, int node, int item)
{
  AdjIterator it;
  NodeItem  *remove, *previous;

  remove = NULL;
  previous = NULL;

  for ( it = adjFirst(adj,node); adjValid(it); it = adjNext(it) ) {
    if (adjItem(it)==item) {
      remove = it;
      break;
    }else{
      previous = it;
    }
  }

  if (remove == NULL) return NULL;
 
  if ( previous == NULL ) {
    adj->first[node] = remove->next;
  }else{
    previous->next = remove->next;
  }

  remove->item = EMPTY;
  remove->next = adj->blank;
  adj->blank = remove;

  return adj;
}

AdjIterator adjGetCurrent( Adj *adj ){
  return adj->current;
}

Adj *adjSetCurrent( Adj *adj, AdjIterator iterator ){
  adj->current = iterator;
  return adj;
}

GridBool adjExists( Adj *adj, int node, int item )
{
  AdjIterator it;
  GridBool exist;
  exist = FALSE;
  for ( it = adjFirst(adj,node); 
	!exist && adjValid(it); 
	it = adjNext(it)) 
    exist = (item == adjItem(it));
  return exist;
}

int adjDegree(Adj *adj, int node )
{
  AdjIterator it;
  int degree;
  degree =0;
  for ( it = adjFirst(adj,node) ; adjValid(it); it = adjNext(it)) degree++;
  return degree;
}
