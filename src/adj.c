
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "adj.h"

typedef struct NODE2ITEM NODE2ITEM;

struct NODE2ITEM {
  int item;
  NODE2ITEM *next;
};

struct Adj {
  int nnode;
  NODE2ITEM *node2item;
  NODE2ITEM **first;
  NODE2ITEM *current;
  NODE2ITEM *blank;
};

Adj* adjCreate( int nnode, int perNode )
{
  int i, nlist;
  Adj *adj;

  adj = malloc( sizeof(Adj) );

  adj->nnode   = nnode;
  nlist = adj->nnode * perNode;

  adj->node2item = (NODE2ITEM *)malloc( nlist * sizeof(NODE2ITEM));

  for ( i=0 ; i<nlist-1 ; i++ ) { // pointer majic?
    adj->node2item[i].item = EMPTY;
    adj->node2item[i].next = &adj->node2item[i+1];
  }
  adj->node2item[nlist-1].item = EMPTY;
  adj->node2item[nlist-1].next = NULL;

  adj->blank = adj->node2item;

  adj->current = NULL;

  adj->first = malloc( perNode * adj->nnode * sizeof(NODE2ITEM) );

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
  NODE2ITEM *new;
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
  NODE2ITEM *remove, *previous;
  remove = NULL;

  for ( adjFirst(adj,node); adjValid(adj); adjNext(adj) ) 
    if (adjCurrent(adj)==item) 
      remove = adj->current;

  if (remove == NULL) return NULL;
 
  previous = NULL;

  for ( adjFirst(adj,node); adjValid(adj); adjNext(adj) ) 
    if (adj->current != NULL && adj->current->next == remove) 
      previous = adj->current;
  
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

bool adjValid( Adj *adj )
{
  return (adj->current != NULL);
}

bool adjMore( Adj *adj )
{
  return ( (adj->current != NULL) && (adj->current->next != NULL) );
}

Adj *adjFirst( Adj *adj, int node )
{
  if ( node < adj->nnode ) {
    adj->current = adj->first[node];
  }else{
    adj->current = NULL;
    return NULL;
  }
  return adj;
}

int adjCurrent( Adj *adj )
{
  if (adj->current == NULL ) return EMPTY;
  return adj->current->item;
}

void adjNext( Adj *adj )
{
  if ( adj->current != NULL ) adj->current = adj->current->next;
}

bool adjExists( Adj *adj, int node, int item )
{
  bool exist;
  exist = FALSE;
  for ( adjFirst(adj,node); 
	!exist && adjValid(adj); 
	adjNext(adj)) 
    exist = (item == adjCurrent(adj));
  return exist;
}

int adjDegree(Adj *adj, int node )
{
  int degree;
  degree =0;
  for ( adjFirst(adj,node) ; adjValid(adj); adjNext(adj)) degree++;
  return degree;
}
