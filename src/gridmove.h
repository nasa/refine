
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRIDMOVE_H
#define GRIDMOVE_H

#include "refine_defs.h"
#include "grid.h"

BEGIN_C_DECLORATION

typedef struct GridMove GridMove;
struct GridMove {
  Grid *grid;
  void *gridRubyVALUEusedForGC;
  double *displacement;
  GridBool *specified;

  int *c2e;
  int nsprings, *springs;
  double *xyz;
  double *k;
  double *source;

  double *ksum;
  double *kxyz;

  int *rowStart;
  int *compRow;
  double *a, *lu;
  double *dxyz;
};

GridMove *gridmoveCreate(Grid *);
Grid *gridmoveGrid(GridMove *);
void gridmoveFree(GridMove *);
void gridmovePack(void *voidGridMove, 
		  int nnode, int maxnode, int *nodeo2n,
		  int ncell, int maxcell, int *cello2n,
		  int nface, int maxface, int *faceo2n,
		  int nedge, int maxedge, int *edgeo2n);
void gridmoveSortNode(void *voidGridMove, int maxnode, int *o2n);
void gridmoveReallocator(void *voidGridMove, int reallocType, 
			 int lastSize, int newSize);
void gridmoveGridHasBeenFreed(void *voidGridMove );

GridMove *gridmoveDisplace(GridMove *, int node, double *displace);
GridMove *gridmoveDisplacement(GridMove *, int node, double *displacement);
GridBool gridmoveSpecified(GridMove *, int node);

GridMove *gridmoveCellFaceNormals(GridMove *, double *xyz, int *nodes, 
				  double normals[4][3]);
GridMove *gridmoveSpringConstant(GridMove *, double *xyz, int nsprings, 
				 double *k, int *springs, int *c2e);

GridMove *gridmoveSpringRelaxationStartUp(GridMove *);
GridMove *gridmoveSpringRelaxationStartStep(GridMove *, double position);
GridMove *gridmoveSpringRelaxationSubIteration(GridMove *, double *residual2);
GridMove *gridmoveSpringRelaxationShutDown(GridMove *);

GridMove *gridmoveSpringRelaxation(GridMove *, int nsteps, int subIterations);

GridMove *gridmoveComputeC2E(GridMove *, int *nedge, int *c2e);
GridMove *gridmoveComputeSpringsWithC2E(GridMove *, int *c2e, int *springs);
GridMove *gridmoveSprings(GridMove *, int *nsprings, int **springs);

GridMove *gridmoveApplyDisplacements(GridMove *);
GridMove *gridmoveProjectionDisplacements(GridMove *);

GridMove *gridmoveDataLeadingDimension(GridMove *, int *ndim );
GridMove *gridmoveLoadFortranNodeData(GridMove *, int nnode, 
				     int *nodes, double *data);
GridMove *gridmoveSetFortranNodeData(GridMove *, int nnode, 
				    int *nodes, double *data);

GridMove *gridmoveInitializeCompRow(GridMove *);
int gridmoveRowStart(GridMove *, int row);
int gridmoveNNZ(GridMove *);
int gridmoveRowNode(GridMove *, int entry);
int gridmoveRowEntry(GridMove *, int row, int node);

GridMove *gridmoveElasticRelaxationStartUp(GridMove *);
GridMove *gridmoveElasticRelaxationStartStep(GridMove *, double position);
GridMove *gridmoveElasticRelaxationSubIteration(GridMove *, double *residual2);
GridMove *gridmoveElasticRelaxationShutDown(GridMove *);

GridMove *gridmoveElasticRelaxationDumpA(GridMove *);

GridMove *gridmoveElasticRelaxation(GridMove *, int nsteps, int subIterations);

END_C_DECLORATION

#endif /* GRIDMOVE_H */
