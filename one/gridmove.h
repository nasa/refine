
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

  int relaxationScheme;

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
  double *lastxyz;
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
GridMove *gridmoveInitializeMPITest(GridMove * );
GridMove *gridmoveCompleteMPITest(GridMove * );
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
GridMove *gridmoveElasticRelax_Old_StartStep(GridMove *, double position);
GridMove *gridmoveElasticRelaxationStartStep(GridMove *, double position);
GridMove *gridmoveElasticRelaxationSubIteration(GridMove *, double *residual2);
GridMove *gridmoveElasticRelaxationEvenOdd(GridMove *, double *residual2);
GridMove *gridmoveElasticRelaxationShutDown(GridMove *);

GridMove *gridmoveElasticRelaxationDumpA(GridMove *);

GridMove *gridmoveElasticRelaxation(GridMove *, int nsteps, int subIterations);

#define gridmoveEMPTY_SCHEME (EMPTY)
#define gridmoveELASTIC_SCHEME (1)
#define gridmoveSPRING_SCHEME (2)
GridMove *gridmoveRelaxationStartUp(GridMove *, int relaxationScheme );
GridMove *gridmoveRelaxationStartStep(GridMove *, double position);
GridMove *gridmoveRelaxationSubIteration(GridMove *, double *residual2);
GridMove *gridmoveRelaxationShutDown(GridMove *);

GridMove *gridmoveRelaxation(GridMove *, int relaxationScheme, 
			     int nsteps, int subIterations);

END_C_DECLORATION

#endif /* GRIDMOVE_H */
