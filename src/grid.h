
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRID_H
#define GRID_H

#include "master_header.h"

BEGIN_C_DECLORATION

typedef struct Grid Grid;

Grid *gridCreate(int maxnode, int maxcell, int maxface, int maxedge );
Grid *gridImport(int maxnode, int nnode, 
		 int maxface, int nface, 
		 int maxcell, int ncell,
		 int maxedge,
		 double *xyz, int *f2n, int *faceId, int *c2n );
Grid *gridImportFAST( char *filename );
Grid *gridExport(Grid *g, int *nnode, int *nface, int *ncell,
		 double **xyz, int **f2n, int **faceId, int **c2n );
void gridFree(Grid *g);

int gridMaxNode(Grid *g);
int gridNNode(Grid *g);
int gridMaxCell(Grid *g);
int gridNCell(Grid *g);
int gridMaxFace(Grid *g);
int gridNFace(Grid *g);
int gridMaxEdge(Grid *g);
int gridNEdge(Grid *g);
int gridCellDegree(Grid *g, int nodeIndex);

Grid *gridAddCell(Grid *g, int n0, int n1, int n2, int n3 );
Grid *gridRemoveCell(Grid *g, int cellId );
Grid *gridCell(Grid *g, int cellId, int *nodes );

Grid *gridAddFace(Grid *g, int n0, int n1, int n2, int faceId );
Grid *gridRemoveFace(Grid *g, int face );
int gridFindFace(Grid *g, int n0, int n1, int n2 );
int gridFaceId(Grid *g, int n0, int n1, int n2 );

Grid *gridAddEdge(Grid *g, int n0, int n1, 
		  int edgeId, double t0, double t1 );
Grid *gridRemoveEdge(Grid *g, int edge );
int gridFindEdge(Grid *g, int n0, int n1 );
int gridEdgeId(Grid *g, int n0, int n1 );

int gridGeomCurveSize( Grid *g, int edgeId, int startNode );
Grid *gridGeomCurve( Grid *g, int edgeId, int startNode, int *curve );
Grid *gridGeomCurveT( Grid *g, int edgeId, int startNode, double *curve );

Grid *gridMakeGem(Grid *g, int n0, int n1 );
int gridNGem(Grid *g );
int gridGem(Grid *g, int index );

Grid *gridOrient(Grid *g, int *cell, int *nodes );
Grid *gridEquator(Grid *g, int n0, int n1 );
int gridNEqu(Grid *g );
int gridEqu(Grid *g, int index );
Grid *gridThrash(Grid *g );

Grid *gridSplitEdge(Grid *g, int n0, int n1 );

int gridAddNode(Grid *g, double x, double y, double z );

double gridVolume(Grid *g, int *nodes );
double gridAR(Grid *g, int *nodes );
double gridMinVolume(Grid *g);
double gridMinAR(Grid *g);

int gridFindCellWithFace(Grid *g, int face );
bool gridRightHandedFace(Grid *g, int face );
bool gridRightHandedBoundary(Grid *g );

END_C_DECLORATION

#endif /* GRID_H */
