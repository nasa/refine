
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include "adj.h"

struct Adj {
  int nnode;
};

Adj* adjCreate( int nnode, int perNode )
{
  Adj *adj;
  adj = malloc(sizeof(Adj));
  adj->nnode   = nnode;
  return adj;
}

void adjFree(Adj *adj)
{
  free(adj);
}

int adjNNode(Adj *adj)
{
  return adj->nnode;
}
