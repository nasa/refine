
/* Michael A. Park
 * Computational AeroSciences Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:Mike.Park@NASA.Gov 
 */
  
/* $Id$ */

#ifndef GRIDEDGER_H
#define GRIDEDGER_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

typedef struct GridEdger GridEdger;
struct GridEdger {
  Grid *grid;
  void *gridRubyVALUEusedForGC;

  int edgeId;

  int ideal_nodes;
  double *t;

  int total_unused;
  int *unused;
};

GridEdger *gridedgerCreate(Grid *, int edgeId);
Grid *gridedgerGrid(GridEdger *);
void gridedgerFree(GridEdger *);
void gridedgerPack(void *voidGridEdger, 
		  int nnode, int maxnode, int *nodeo2n,
		  int ncell, int maxcell, int *cello2n,
		  int nface, int maxface, int *faceo2n,
		  int nedge, int maxedge, int *edgeo2n);
void gridedgerSortNode(void *voidGridEdger, int maxnode, int *o2n);
void gridedgerReallocator(void *voidGridEdger, int reallocType, 
			 int lastSize, int newSize);
void gridedgerGridHasBeenFreed(void *voidGridEdger );

int gridedgerEdgeId(GridEdger *);
int gridedgerIdealNodes(GridEdger *);
GridEdger *gridedgerIdealNodeT(GridEdger *, int node, double *t );

GridEdger *gridedgerSupportingSegment(GridEdger *, double t, double *segment );

GridEdger *gridedgerDiscreteSegmentAndRatio(GridEdger *, double segment, 
					    int *discrete_segment, 
					    double *segment_ratio );
GridEdger *gridedgerSegmentT(GridEdger *, double segment, double *t );
GridEdger *gridedgerSegmentMap(GridEdger *, double segment, double *map );
GridEdger *gridedgerLengthToS(GridEdger *, double segment, double length, 
			      double *next_s );

GridEdger *gridedgerDiscretize(GridEdger *, double length );
GridEdger *gridedgerDiscretizeEvenly(GridEdger * );

GridEdger *gridedgerInsert(GridEdger *);
GridEdger *gridedgerRemoveUnused(GridEdger *);

END_C_DECLORATION

#endif /* GRIDEDGER_H */
