
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "adj.h"

struct NodeItem {
  int item;
  NodeItem *next;
};

struct Adj {
  int nnode;
  NodeItem *node2item;
  NodeItem **first;
  NodeItem *current;
  NodeItem *blank;
};

Adj* adjCreate( int nnode, int perNode )
{
  int i, nlist;
  Adj *adj;

  adj = malloc( sizeof(Adj) );

  adj->nnode   = nnode;
  nlist = adj->nnode * perNode;

  adj->node2item = (NodeItem *)malloc( nlist * sizeof(NodeItem));

  for ( i=0 ; i<nlist-1 ; i++ ) { // pointer majic?
    adj->node2item[i].item = EMPTY;
    adj->node2item[i].next = &adj->node2item[i+1];
  }
  adj->node2item[nlist-1].item = EMPTY;
  adj->node2item[nlist-1].next = NULL;

  adj->blank = adj->node2item;

  adj->current = NULL;

  adj->first = malloc( perNode * adj->nnode * sizeof(NodeItem) );

  for ( i=0 ; i<adj->nnode; i++ ) adj->first[i] = NULL; 

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

Adj *adjRegister( Adj *adj, int node, int item )
{
  NodeItem *new;
  if (node>=adj->nnode) return NULL;
  if (adj->blank == NULL) return NULL;
  new = adj->blank;
  adj->blank = adj->blank->next;
  new->next = adj->first[node];
  adj->first[node] = new;
  new->item = item;
  return adj;
}

Adj* adjRemove(Adj *adj, int node, int item)
{
  NodeItem *it, *remove, *previous;
  remove = NULL;

  for ( it = adjFirst(adj,node); adjValid(it); it = adjNext(it) ) 
    if (adjItem(it)==item) remove = it;

  if (remove == NULL) return NULL;
 
  previous = NULL;

  for ( it = adjFirst(adj,node); adjValid(it); it = adjNext(it) ) 
    if (it != NULL && it->next == remove) previous = it;
  
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

bool adjValid( NodeItem *iterator )
{
  return (iterator != NULL);
}

bool adjMore( NodeItem *iterator )
{
  return ( (iterator != NULL) && (iterator->next != NULL) );
}

NodeItem *adjFirst( Adj *adj, int node )
{
  if ( node < adj->nnode ) return adj->first[node];
  return NULL;
}

int adjItem( NodeItem *iterator )
{
  if ( iterator == NULL ) return EMPTY;
  return iterator->item;
}

NodeItem *adjNext( NodeItem *iterator )
{
  if ( iterator != NULL ) return iterator->next;
  return NULL;
}

NodeItem *adjGetCurrent( Adj *adj ){
  return adj->current;
}

Adj *adjSetCurrent( Adj *adj, NodeItem *iterator ){
  adj->current = iterator;
  return adj;
}

bool adjExists( Adj *adj, int node, int item )
{
  NodeItem *it;
  bool exist;
  exist = FALSE;
  for ( it = adjFirst(adj,node); 
	!exist && adjValid(it); 
	it = adjNext(it)) 
    exist = (item == adjItem(it));
  return exist;
}

int adjDegree(Adj *adj, int node )
{
  NodeItem *it;
  int degree;
  degree =0;
  for ( it = adjFirst(adj,node) ; adjValid(it); it = adjNext(it)) degree++;
  return degree;
}
