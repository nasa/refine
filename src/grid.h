
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef GRID_H
#define GRID_H

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <values.h>
#include "master_header.h"
#include "adj.h"
#include "line.h"

BEGIN_C_DECLORATION

#define MAXDEG 200

typedef struct Prism Prism;
struct Prism {
  int nodes[6];
};

typedef struct Pyramid Pyramid;
struct Pyramid {
  int nodes[5];
};

typedef struct Quad Quad;
struct Quad {
  int nodes[4];
  int faceId;
};

typedef struct Grid Grid;
struct Grid {
  int maxnode, nnode;
  int blanknode;
  double *xyz;
  double *map;
  bool *frozen;

  int *nodeGlobal;
  int *part;
  int nsorted;
  int *sortedGlobal, *sortedLocal;

  int maxcell, ncell;
  int blankc2n;
  int *c2n;
  Adj *cellAdj;
  int *cellGlobal;

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

  int nprism, maxprism;
  Prism *prism;
  int *prismDeg;

  int npyramid, maxpyramid;
  Pyramid *pyramid;

  int nquad, maxquad;
  Quad *quad;

  int partId;
  int globalNNode;
  int globalNCell;

  int nGeomNode;
  int nGeomEdge;
  int nGeomFace;
  int *geomEdge;

  int ngem;
  int gem[MAXDEG];

  int nequ;
  int equ[MAXDEG];

  int degAR;
  double AR[MAXDEG];
  double dARdX[3*MAXDEG];

  FILE *tecplotFile;

  void (*renumberFunc)(void *renumberData, int *o2n);
  void *renumberData;

  Lines *lines;
};

Grid *gridCreate(int maxnode, int maxcell, int maxface, int maxedge );
Grid *gridImport(int maxnode, int nnode, 
		 int maxface, int nface, 
		 int maxcell, int ncell,
		 int maxedge,
		 double *xyz, int *f2n, int *faceId, int *c2n );
Grid *gridImportFAST( char *filename );
Grid *gridExportFAST(Grid *g, char *filename );
Grid *gridExportAFLR3(Grid *g, char *filename );
Grid *gridExport(Grid *g, int *nnode, int *nface, int *ncell,
		 double **xyz, int **f2n, int **faceId, int **c2n );
Grid *gridImportAdapt(Grid *g, char *filename );
Grid *gridAttachNodeSorter(Grid *g, 
			   void (*renumberFunc)(void *renumberData, int *o2n),
			   void *renumberData );
Grid *gridDetachNodeSorter(Grid *g );
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
int gridNPrism(Grid *g);
int gridNPyramid(Grid *g);
int gridNQuad(Grid *g);
int gridPartId(Grid *g);
Grid *gridSetPartId(Grid *g, int partId );
int gridGlobalNNode(Grid *g);
Grid *gridSetGlobalNNode(Grid *g, int nglobal );
int gridGlobalNCell(Grid *g);
Grid *gridSetGlobalNCell(Grid *g, int nglobal );
int gridCellDegree(Grid *g, int nodeIndex);
int gridCellGlobal(Grid *g, int cellIndex);
Grid *gridSetCellGlobal(Grid *g, int cellIndex, int globalIndex );

int gridAddCell(Grid *g, int n0, int n1, int n2, int n3 );
Grid *gridRemoveCell(Grid *g, int cellId );
#define gridCellAdj(grid) (NULL==grid?NULL:grid->cellAdj)
Grid *gridReconnectAllCell(Grid *g, int oldNode, int newNode );
Grid *gridCell(Grid *g, int cellId, int *nodes );
bool gridCellValid(Grid *g, int cellId );
bool gridCellEdge(Grid *g, int node0, int node1 );
bool gridCellFace(Grid *g, int node0, int node1, int node2 );
int gridFindOtherCellWith3Nodes(Grid *g, int node0, int node1, int node2,
				int currentCell );
int gridFindCellWithFace(Grid *g, int face );
int gridFindCell(Grid *g, int *nodes );
Grid *gridCheckCellConnections(Grid *g);
Grid *gridDeleteThawedCells(Grid *g);

int gridAddFace(Grid *g, int n0, int n1, int n2, int faceId );
int gridAddFaceUV(Grid *g, 
		  int n0, double u0, double v0,
		  int n1, double u1, double v1,
		  int n2, double u2, double v2, int faceId );
Grid *gridRemoveFace(Grid *g, int face );
#define gridFaceAdj(grid) (NULL==grid?NULL:grid->faceAdj)
int gridFindFace(Grid *g, int n0, int n1, int n2 );
int gridFaceId(Grid *g, int n0, int n1, int n2 );
Grid *gridReconnectAllFace(Grid *g, int oldNode, int newNode );
Grid *gridFace(Grid *g, int face, int *nodes, int *id );
Grid *gridDeleteThawedFaces(Grid *g, int faceId );
int gridNThawedFaces(Grid *g, int faceId );

Grid *gridNodeUV(Grid *g, int node, int faceId, double *uv );
Grid *gridSetNodeUV(Grid *g, int node, int faceId, double u, double v );
double gridNodeU(Grid *grid, int node, int faceId);
double gridNodeV(Grid *grid, int node, int faceId);
Grid *gridNodeT(Grid *g, int node, int edgeId, double *t );
Grid *gridSetNodeT(Grid *g, int node, int edgeId, double t );

int gridAddEdge(Grid *g, int n0, int n1, 
		int edgeId, double t0, double t1 );
Grid *gridRemoveEdge(Grid *g, int edge );
#define gridEdgeAdj(grid) (NULL==grid?NULL:grid->edgeAdj)
int gridFindEdge(Grid *g, int n0, int n1 );
int gridEdgeId(Grid *g, int n0, int n1 );
Grid *gridReconnectAllEdge(Grid *g, int oldNode, int newNode );
Grid *gridEdge(Grid *g, int edge, int *nodes, int *edgeId );
Grid *gridDeleteThawedEdgeSegments(Grid *g, int edgeId );
int gridNThawedEdgeSegments(Grid *g, int edgeId );

int gridGeomCurveSize( Grid *g, int edgeId, int startNode );
Grid *gridGeomCurve( Grid *g, int edgeId, int startNode, int *curve );
Grid *gridGeomCurveT( Grid *g, int edgeId, int startNode, double *curve );

int gridNFrozen( Grid *g );
bool gridNodeFrozen( Grid *g, int node );
Grid *gridFreezeNode( Grid *g, int node );
Grid *gridThawNode( Grid *g, int node );
Grid *gridFreezeAll( Grid *g );
Grid *gridThawAll( Grid *g );
Grid *gridThawSphere( Grid *g, double x, double y, double z, double r );
Grid *gridThawNearBC( Grid *g, double r, int faceId );
Grid *gridFreezeBCFace( Grid *g, int faceId );

Grid *gridMakeGem(Grid *g, int n0, int n1 );
int gridNGem(Grid *g );
int gridGem(Grid *g, int index );
Grid *gridRemoveGem(Grid *g);
bool gridGemIsLocal(Grid *g);

Grid *gridOrient(Grid *g, int *cell, int *nodes );
Grid *gridEquator(Grid *g, int n0, int n1 );
int gridNEqu(Grid *g );
int gridEqu(Grid *g, int index );
bool gridContinuousEquator(Grid *g);
Grid *gridCycleEquator( Grid *g );

int gridAddNode(Grid *g, double x, double y, double z );
Grid *gridRemoveNode(Grid *g, int node );
#define gridValidNode(grid,node) \
(node>-1 && node<grid->maxnode && DBL_MAX!=grid->xyz[0+3*node])


Grid *gridNodeXYZ(Grid *g, int node, double *xyz );
#define gridNodeXYZPointer(grid, node) (&grid->xyz[3*node])
#define gridNodeXYZEntry(grid, node, ixyz) (grid->xyz[ixyz+3*node])
Grid *gridSetNodeXYZ(Grid *g, int node, double *xyz );
Grid *gridDeleteNodesNotUsed(Grid *);

int gridNodeGlobal(Grid *g, int node );
int gridGlobal2Local(Grid *g, int global );
Grid *gridSetNodeGlobal(Grid *g, int node, int global );
int gridNodePart(Grid *g, int node );
Grid *gridSetNodePart(Grid *g, int node, int part );
bool gridNodeLocal(Grid *g, int node );
bool gridNodeGhost(Grid *g, int node );

int gridNGeomNode(Grid *g);
Grid *gridSetNGeomNode(Grid *g, int nGeomNode);
int gridNGeomEdge(Grid *g);
Grid *gridSetNGeomEdge(Grid *g, int nGeomEdge);
int gridNGeomFace(Grid *g);
Grid *gridSetNGeomFace(Grid *g, int nGeomFace);

Grid *gridAddGeomEdge(Grid *g, int edgeId, int n0, int n1);
int gridGeomEdgeStart( Grid *g, int edgeId );
int gridGeomEdgeEnd( Grid *g, int edgeId );
int gridGeomEdgeSize( Grid *g, int edgeId );
Grid *gridGeomEdge( Grid *g, int edgeId, int *curve );
int gridFrozenEdgeEndPoint( Grid *g, int edgeId, int startNode );

bool gridGeometryNode(Grid *g, int node);
bool gridGeometryEdge(Grid *g, int node);
bool gridGeometryFace(Grid *g, int node);

Grid *gridAddPrism(Grid *g, int n0, int n1, int n2, int n3, int n4, int n5);
Grid *gridPrism(Grid *g, int prismIndex, int *nodes);

Grid *gridAddPyramid(Grid *g, int n0, int n1, int n2, int n3, int n4);
Grid *gridPyramid(Grid *g, int pyramidIndex, int *nodes);

Grid *gridAddQuad(Grid *g, int n0, int n1, int n2, int n3, int faceId );
Grid *gridQuad(Grid *g, int quadIndex, int *nodes, int *faceId );

Grid *gridMap(Grid *g, int node, double *map);
#define gridMapPointer(grid, node) (&grid->map[6*node])

Grid *gridSetMap(Grid *g, int node,
		 double m11, double m12, double m13,
		             double m22, double m23,
		                         double m33);

int gridStoredARDegree( Grid *g );
Grid *gridClearStoredAR( Grid *g );
Grid *gridAddStoredAR( Grid *g, double AR, double *dARdX );
double gridStoredAR( Grid *g, int index );
Grid *gridStoredARDerivative( Grid *g, int index, double *dARdX );

#define gridLines(grid) (NULL==grid?NULL:grid->lines)
Grid *gridFreezeLinesNodes(Grid *g);
Grid *gridReportLinesLocation(Grid *g);

Grid *gridCopyAboutY0(Grid *g);

END_C_DECLORATION

#endif /* GRID_H */
