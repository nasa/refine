
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
Grid *gridExportFAST(Grid *g, char *filename );
Grid *gridExport(Grid *g, int *nnode, int *nface, int *ncell,
		 double **xyz, int **f2n, int **faceId, int **c2n );
Grid *gridImportAdapt(Grid *g, char *filename );
Grid *gridPack(Grid *g);
Grid *gridSortNodeGridEx(Grid *g);
void gridFree(Grid *g);

Grid *gridWriteTecplotSurfaceZone(Grid *g);

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
Grid *gridAddFaceUV(Grid *g, 
		    int n0, double u0, double v0,
		    int n1, double u1, double v1,
		    int n2, double u2, double v2, int faceId );
Grid *gridRemoveFace(Grid *g, int face );
int gridFindFace(Grid *g, int n0, int n1, int n2 );
int gridFaceId(Grid *g, int n0, int n1, int n2 );

Grid *gridNodeUV(Grid *g, int node, int faceId, double *uv );
Grid *gridSetNodeUV(Grid *g, int node, int faceId, double u, double v );
double gridNodeU(Grid *grid, int node, int faceId);
double gridNodeV(Grid *grid, int node, int faceId);
Grid *gridNodeT(Grid *g, int node, int edgeId, double *t );
Grid *gridSetNodeT(Grid *g, int node, int edgeId, double t );

Grid *gridAddEdge(Grid *g, int n0, int n1, 
		  int edgeId, double t0, double t1 );
Grid *gridRemoveEdge(Grid *g, int edge );
int gridFindEdge(Grid *g, int n0, int n1 );
int gridEdgeId(Grid *g, int n0, int n1 );

int gridGeomCurveSize( Grid *g, int edgeId, int startNode );
Grid *gridGeomCurve( Grid *g, int edgeId, int startNode, int *curve );
Grid *gridGeomCurveT( Grid *g, int edgeId, int startNode, double *curve );

bool gridNodeFrozen( Grid *g, int node );
Grid *gridFreezeNode( Grid *g, int node );
Grid *gridThawNode( Grid *g, int node );
Grid *gridFreezeAll( Grid *g );
Grid *gridThawAll( Grid *g );

Grid *gridMakeGem(Grid *g, int n0, int n1 );
int gridNGem(Grid *g );
int gridGem(Grid *g, int index );

Grid *gridOrient(Grid *g, int *cell, int *nodes );
Grid *gridEquator(Grid *g, int n0, int n1 );
int gridNEqu(Grid *g );
int gridEqu(Grid *g, int index );

int gridAddNode(Grid *g, double x, double y, double z );
Grid *gridRemoveNode(Grid *g, int node );
bool gridValidNode(Grid *g, int node );
Grid *gridNodeXYZ(Grid *g, int node, double *xyz );

int gridFindCellWithFace(Grid *g, int face );

int gridNGeomNode(Grid *g);
Grid *gridSetNGeomNode(Grid *g, int nGeomNode);

int gridNGeomEdge(Grid *g);
Grid *gridSetNGeomEdge(Grid *g, int nGeomEdge);
Grid *gridAddGeomEdge(Grid *g, int edge, int n0, int n1);
int gridGeomEdgeSize( Grid *g, int edge );
Grid *gridGeomEdge( Grid *g, int edge, int *curve );

bool gridGeometryNode(Grid *g, int node);
bool gridGeometryEdge(Grid *g, int node);
bool gridGeometryFace(Grid *g, int node);

END_C_DECLORATION

#endif /* GRID_H */
