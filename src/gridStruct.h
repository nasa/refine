
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDSTRUCT_H
#define GRIDSTRUCT_H

#include "master_header.h"

BEGIN_C_DECLORATION

#define MAXDEG 200

struct Grid {
  int maxnode, nnode;
  double *xyz;

  int maxcell, ncell;
  int blankc2n;
  int *c2n;
  Adj *cellAdj;

  int maxface, nface;
  int blankf2n;
  int *f2n;
  int *faceId;
  double *faceU, *faceV;
  Adj *faceAdj;

  int maxedge, nedge;
  int blanke2n;
  int *e2n;
  int *edgeId;
  double *edgeT;
  Adj *edgeAdj;

  int ngem;
  int gem[MAXDEG];

  int nequ;
  int equ[MAXDEG];
};

END_C_DECLORATION

#endif /* GRIDSTRUCT_H */
