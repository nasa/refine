
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

Adj* adjCreate( int nnode, int nadj )
{
  int i;
  Adj *adj;

  adj = malloc( sizeof(Adj) );

  adj->nnode   = MAX(nnode,1);
  nadj         = MAX(nadj,1);

  adj->node2item = (NodeItem *)malloc( nadj * sizeof(NodeItem));

  for ( i=0 ; i<nadj-1 ; i++ ) { // pointer majic?
    adj->node2item[i].item = EMPTY;
    adj->node2item[i].next = &adj->node2item[i+1];
  }
  adj->node2item[nadj-1].item = EMPTY;
  adj->node2item[nadj-1].next = NULL;

  adj->blank = adj->node2item;

  adj->current = NULL;

  adj->first = malloc( adj->nnode * sizeof(NodeItem*) );

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
  AdjIterator it;
  NodeItem  *remove, *previous;
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

bool adjValid( AdjIterator iterator )
{
  return (iterator != NULL);
}

bool adjMore( AdjIterator iterator )
{
  return ( (iterator != NULL) && (iterator->next != NULL) );
}

AdjIterator adjFirst( Adj *adj, int node )
{
  if ( node < adj->nnode ) return adj->first[node];
  return NULL;
}

int adjItem( AdjIterator iterator )
{
  if ( iterator == NULL ) return EMPTY;
  return iterator->item;
}

AdjIterator adjNext( AdjIterator iterator )
{
  if ( iterator != NULL ) return iterator->next;
  return NULL;
}

AdjIterator adjGetCurrent( Adj *adj ){
  return adj->current;
}

Adj *adjSetCurrent( Adj *adj, AdjIterator iterator ){
  adj->current = iterator;
  return adj;
}

bool adjExists( Adj *adj, int node, int item )
{
  AdjIterator it;
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
  AdjIterator it;
  int degree;
  degree =0;
  for ( it = adjFirst(adj,node) ; adjValid(it); it = adjNext(it)) degree++;
  return degree;
}
