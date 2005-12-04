
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
#include "refine_defs.h"
#include "adj.h"
#include "line.h"
#include "queue.h"

BEGIN_C_DECLORATION

#define MAXDEG 200
#define MAXFACEIDDEG 100
#define MAXEDGEIDDEG 100

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
  GridBool *frozen;
  int naux;
  double *aux;

  Grid *child;
  Grid *parent;

  int *child_reference;
  int *nodeGlobal;
  int *part;
  int nsorted;
  int *sortedGlobal, *sortedLocal;
  int maxUnusedNodeGlobal, nUnusedNodeGlobal, *unusedNodeGlobal;

  int maxcell, ncell;
  int blankc2n;
  int *c2n;
  Adj *cellAdj;
  int *cellGlobal;
  int maxUnusedCellGlobal, nUnusedCellGlobal, *unusedCellGlobal;

  GridBool constrain_surface_node;

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
  int *geomNode;
  int *geomEdge;

  int ngem;
  int gem[MAXDEG];

  int nequ;
  int equ[2*MAXDEG];

  int nconn;
  int *cell2conn;
  int *conn2node;

  int costDegree;
  double storedCost[MAXDEG];
  double storedCostDerivative[3*MAXDEG];

  FILE *tecplotGeomFile;
  FILE *tecplotScalarFile;

  int costFunction;
  int costConstraint;
  double min_allowed_insert_cost;
  double min_allowed_surface_smooth_cost;

  void (*packFunc)(void *packData, 
		   int nnode, int maxnode, int *nodeo2n,
		   int ncell, int maxcell, int *cello2n,
		   int nface, int maxface, int *faceo2n,
		   int nedge, int maxedge, int *edgeo2n);
  void *packData;

  void (*renumberFunc)(void *renumberData, int maxnode, int *o2n);
  void *renumberData;

  void (*reallocFunc)(void *reallocData, int reallocType, 
		      int lastSize, int newSize);
  void *reallocData;

  void (*freeNotificationFunc)(void *freeNotificationData);
  void *freeNotificationData;

  Lines *lines;

  int phase;

  int model; /* CAPRI 2 Model Id */
};

Grid *gridCreate(int maxnode, int maxcell, int maxface, int maxedge );
Grid *gridImport(int maxnode, int nnode, 
		 int maxface, int nface, 
		 int maxcell, int ncell,
		 int maxedge,
		 double *xyz, int *f2n, int *faceId, int *c2n );
Grid *gridDup(Grid *g);
Grid *gridImportFAST( char *filename );
Grid *gridExportFAST(Grid *g, char *filename );
Grid *gridExportAFLR3(Grid *g, char *filename );
Grid *gridExport(Grid *g, int *nnode, int *nface, int *ncell,
		 double **xyz, int **f2n, int **faceId, int **c2n );
Grid *gridImportAdapt(Grid *g, char *filename );

Grid *gridAttachPacker(Grid *g, 
		       void (*packFunc)(void *packData, 
					int nnode, int maxnode, int *nodeo2n,
					int ncell, int maxcell, int *cello2n,
					int nface, int maxface, int *faceo2n,
					int nedge, int maxedge, int *edgeo2n),
		       void *packData );
Grid *gridDetachPacker(Grid *g );

Grid *gridAttachNodeSorter(Grid *g, 
			   void (*renumberFunc)(void *renumberData, 
						int maxnode, int *o2n),
			   void *renumberData );
Grid *gridDetachNodeSorter(Grid *g );

Grid *gridAttachReallocator(Grid *g, 
			    void (*reallocFunc)(void *reallocData, 
						int reallocType, 
						int lastSize, int newSize), 
			    void *reallocData );
Grid *gridDetachReallocator(Grid *g);
#define gridREALLOC_NODE (1)
#define gridREALLOC_EDGE (2)
#define gridREALLOC_FACE (3)
#define gridREALLOC_CELL (4)
#define gridREALLOC_PRISM (5)
#define gridREALLOC_PYRAMID (6)
#define gridREALLOC_QUAD (7)

Grid *gridAttachFreeNotifier(Grid *g, void (*freeNotificationFunc)
			     (void *freeNotificationData),
			     void *freeNotificationData);
Grid *gridDetachFreeNotifier(Grid *g);


Grid *gridPack(Grid *g);
Grid *gridSortNodeGridEx(Grid *g);
Grid *gridSortNodeFUN3D(Grid *g, int *nnodes0);
Grid *gridRenumber(Grid *g, int *o2n);
void gridFree(Grid *g);

Grid *gridWriteTecplotSurfaceGeom(Grid *g, char *filename );
Grid *gridWriteTecplotGeomFaceUV(Grid *g, char *filename, int id );
Grid *gridWriteTecplotTriangleZone(Grid *g, char *filename,
				   int nnode, double *xyz,
				   int nface, int *f2n);
Grid *gridWriteTecplotComment(Grid *g, char *comment );
Grid *gridWriteTecplotCellGeom(Grid *g, int *nodes, double *scalar,
			       char *filename );
Grid *gridWriteTecplotEquator(Grid *g, int n0, int n1, char *filename );
Grid *gridWriteTecplotEquatorFaces(Grid *g, int n0, int n1, char *filename );
Grid *gridCloseTecplotGeomFile(Grid *g);
/* Warning, call gridSortNodeGridEx before calculating scalar 
 * in gridWriteTecplotSurfaceScalar to avoid a renumbering bug. */
Grid *gridWriteTecplotSurfaceScalar(Grid *g, char *filename, double *scalar );
Grid *gridCloseTecplotScalarFile(Grid *g);
Grid *gridWriteVTK(Grid *g, char *filename );

#define gridMaxNode(grid) (grid->maxnode)
#define gridNNode(grid)   (grid->nnode)
#define gridMaxCell(grid) (grid->maxcell)
#define gridNCell(grid)   (grid->ncell)
#define gridMaxFace(grid) (grid->maxface)
#define gridNFace(grid)   (grid->nface)
#define gridMaxEdge(grid) (grid->maxedge)
#define gridNEdge(grid)   (grid->nedge)
int gridNPrism(Grid *g);
int gridNPyramid(Grid *g);
int gridNQuad(Grid *g);
int gridNAux(Grid *g);
Grid *gridSetNAux(Grid *g, int naux);
double gridAux(Grid *g, int node, int aux);
Grid *gridSetAux(Grid *g, int node, int aux, double value);
Grid *gridInterpolateAux2(Grid *g, int node0, int node1, double ratio,
			  int target);
Grid *gridSetAuxToAverageOfNodes2(Grid *g, int avgNode,
				  int n0, int n1 );
Grid *gridSetAuxToAverageOfNodes3(Grid *g, int avgNode,
				  int n0, int n1, int n2 );
Grid *gridSetAuxToAverageOfNodes4(Grid *g, int avgNode,
				  int n0, int n1, int n2, int n3 );
int gridPartId(Grid *g);
Grid *gridSetPartId(Grid *g, int partId );
int gridGlobalNNode(Grid *g);
Grid *gridSetGlobalNNode(Grid *g, int nglobal );
int gridGlobalNCell(Grid *g);
Grid *gridSetGlobalNCell(Grid *g, int nglobal );
int gridNUnusedNodeGlobal(Grid *g );
int gridNUnusedCellGlobal(Grid *g );
Grid *gridGetUnusedNodeGlobal(Grid *g, int *unused );
Grid *gridGetUnusedCellGlobal(Grid *g, int *unused );
Grid *gridJoinUnusedNodeGlobal(Grid *g, int global );
Grid *gridJoinUnusedCellGlobal(Grid *g, int global );
Grid *gridEliminateUnusedNodeGlobal(Grid *g );
Grid *gridEliminateUnusedCellGlobal(Grid *g );
int gridCellDegree(Grid *g, int nodeIndex);
int gridCellGlobal(Grid *g, int cellIndex);
Grid *gridSetCellGlobal(Grid *g, int cellIndex, int globalIndex );
Grid *gridGlobalShiftCell(Grid *g, int oldncellg, int newncellg, 
			  int celloffset );
#define gridCellHasLocalNode(grid,nodes) ( \
gridNodeLocal(grid,(nodes)[0]) || \
gridNodeLocal(grid,(nodes)[1]) || \
gridNodeLocal(grid,(nodes)[2]) || \
gridNodeLocal(grid,(nodes)[3]) )
#define gridCellHasGhostNode(grid,nodes) ( \
gridNodeGhost(grid,(nodes)[0]) || \
gridNodeGhost(grid,(nodes)[1]) || \
gridNodeGhost(grid,(nodes)[2]) || \
gridNodeGhost(grid,(nodes)[3]) )
#define gridFaceHasLocalNode(grid,n0,n1,n2) ( \
gridNodeLocal(grid,n0) || \
gridNodeLocal(grid,n1) || \
gridNodeLocal(grid,n2) )
#define gridFaceHasGhostNode(grid,n0,n1,n2) ( \
gridNodeGhost(grid,n0) || \
gridNodeGhost(grid,n1) || \
gridNodeGhost(grid,n2) )
#define gridEdgeHasLocalNode(grid,n0,n1) ( \
gridNodeLocal(grid,n0) || \
gridNodeLocal(grid,n1) )
#define gridEdgeHasGhostNode(grid,n0,n1) ( \
gridNodeGhost(grid,n0) || \
gridNodeGhost(grid,n1) )

int gridAddCell(Grid *g, int n0, int n1, int n2, int n3 );
int gridAddCellAndQueue(Grid *g, Queue *, int n0, int n1, int n2, int n3 );
int gridAddCellWithGlobal(Grid *g, int n0, int n1, int n2, int n3, int global );
Grid *gridRemoveCell(Grid *g, int cellId );
Grid *gridRemoveCellAndQueue(Grid *g, Queue *, int cellId );
Grid *gridRemoveCellWithOutGlobal(Grid *g, int cellId );
#define gridCellAdj(grid) (NULL==grid?NULL:grid->cellAdj)
Grid *gridReconnectAllCell(Grid *g, int oldNode, int newNode );
Grid *gridCell(Grid *g, int cellId, int *nodes );
GridBool gridCellValid(Grid *g, int cellId );
GridBool gridCellEdge(Grid *g, int node0, int node1 );
GridBool gridCellFace(Grid *g, int node0, int node1, int node2 );
int gridFindOtherCellWith3Nodes(Grid *g, int node0, int node1, int node2,
				int currentCell );
int gridFindCellWithFace(Grid *g, int face );
int gridFindCell(Grid *g, int *nodes );
int gridFindEnclosingCell(Grid *g, int starting_guess, 
			  double *target, double *bary );
Grid *gridBarycentricCoordinate(Grid *g, double *xyz0, double *xyz1, 
				double *xyz2, double *xyz3, 
				double *target, double *bary );
Grid *gridDeleteThawedCells(Grid *g);

#define gridNConn(grid) (grid->nconn)
int gridCell2Conn(Grid *g, int cell, int index );
Grid *gridConn2Node(Grid *g, int conn, int *nodes );
Grid *gridCreateConn(Grid *g );
Grid *gridEraseConn(Grid *g );

Grid *gridConstrainSurfaceNode(Grid *g);
Grid *gridUnconstrainSurfaceNode(Grid *g);
GridBool gridSurfaceNodeConstrained(Grid *g);

int gridAddFace(Grid *g, int n0, int n1, int n2, int faceId );
int gridAddFaceUV(Grid *g, 
		  int n0, double u0, double v0,
		  int n1, double u1, double v1,
		  int n2, double u2, double v2, int faceId );
int gridAddFaceUVAndQueue(Grid *g, Queue *,
		  int n0, double u0, double v0,
		  int n1, double u1, double v1,
		  int n2, double u2, double v2, int faceId );
Grid *gridRemoveFace(Grid *g, int face );
Grid *gridRemoveFaceAndQueue(Grid *g, Queue *, int face );
#define gridFaceAdj(grid) (NULL==grid?NULL:grid->faceAdj)
int gridFindFace(Grid *g, int n0, int n1, int n2 );
int gridFaceId(Grid *g, int n0, int n1, int n2 );
int gridFindFaceWithNodesUnless(Grid *g, int n0, int n1, int unless_face );
Grid *gridReconnectAllFace(Grid *g, int oldNode, int newNode );
GridBool gridReconnectionOfAllFacesOK(Grid *g, int oldNode, int newNode );
Grid *gridFace(Grid *g, int face, int *nodes, int *id );
Grid *gridDeleteThawedFaces(Grid *g, int faceId );
int gridNThawedFaces(Grid *g, int faceId );

Grid *gridNodeUV(Grid *g, int node, int faceId, double *uv );
Grid *gridSetNodeUV(Grid *g, int node, int faceId, double u, double v );
int gridNodeFaceIdDegree(Grid *g, int node);
Grid *gridNodeFaceId(Grid *g, int node, int maxId, int *ids, int *id );

double gridNodeU(Grid *grid, int node, int faceId);
double gridNodeV(Grid *grid, int node, int faceId);
Grid *gridNodeT(Grid *g, int node, int edgeId, double *t );
Grid *gridSetNodeT(Grid *g, int node, int edgeId, double t );
int gridNodeEdgeIdDegree(Grid *g, int node);
Grid *gridNodeEdgeId(Grid *g, int node, int maxId, int *ids, int *id );

int gridAddEdge(Grid *g, int n0, int n1, 
		int edgeId, double t0, double t1 );
int gridAddEdgeAndQueue(Grid *g, Queue *, int n0, int n1, 
		int edgeId, double t0, double t1 );
int gridAddEdgeInGlobal(Grid *g, int g0, int g1, 
			int edgeId, double t0, double t1 );
Grid *gridRemoveEdge(Grid *g, int edge );
Grid *gridRemoveEdgeAndQueue(Grid *g, Queue *, int edge );
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
GridBool gridNodeFrozen( Grid *g, int node );
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
Grid *gridRemoveGemAndQueue(Grid *g, Queue*);
GridBool gridGemIsAllLocal(Grid *g);
GridBool gridNodeNearGhost(Grid *g, int node );

Grid *gridFaceOppositeCellNode(Grid *g, int *nodes, int node, int *face );
Grid *gridOrient(Grid *g, int *cell, int *nodes );
Grid *gridEquator(Grid *g, int n0, int n1 );
int gridNEqu(Grid *g );
int gridEqu(Grid *g, int index );
GridBool gridContinuousEquator(Grid *g);
Grid *gridCycleEquator( Grid *g );

int gridAddNode(Grid *g, double x, double y, double z );
int gridAddNodeWithGlobal(Grid *g, double x, double y, double z, int global );
Grid *gridRemoveNode(Grid *g, int node );
Grid *gridRemoveNodeWithOutGlobal(Grid *g, int node );
#define gridValidNode(grid,node) \
(node>-1 && node<grid->maxnode && DBL_MAX!=grid->xyz[0+3*node])


Grid *gridNodeXYZ(Grid *g, int node, double *xyz );
#define gridNodeXYZPointer(grid, node) (&grid->xyz[3*node])
#define gridNodeXYZEntry(grid, node, ixyz) (grid->xyz[ixyz+3*node])
Grid *gridSetNodeXYZ(Grid *g, int node, double *xyz );
Grid *gridDeleteNodesNotUsed(Grid *);

int gridNodeGlobal(Grid *g, int node );
Grid *gridCreateSortedGlobal(Grid *g );
int gridGlobal2Local(Grid *g, int global );
Grid *gridSetNodeGlobal(Grid *g, int node, int global );
Grid *gridGlobalShiftNode(Grid *g, int oldnnodeg, int newnnodeg, 
			  int nodeoffset );
Grid *gridRenumberGlobalNodes(Grid *g, int nnode, int *n2o);
int gridNodePart(Grid *g, int node );
Grid *gridSetNodePart(Grid *g, int node, int part );
GridBool gridNodeLocal(Grid *g, int node );
GridBool gridNodeGhost(Grid *g, int node );

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

GridBool gridGeometryNode(Grid *g, int node);
GridBool gridGeometryEdge(Grid *g, int node);
GridBool gridGeometryFace(Grid *g, int node);
GridBool gridGeometryBetweenFace(Grid *g, int node);
/* returns -edgeId, faceId, or 0 for no geometry parent */
int gridParentGeometry(Grid *g, int node0, int node1);

Grid *gridAddPrism(Grid *g, int n0, int n1, int n2, int n3, int n4, int n5);
Grid *gridPrism(Grid *g, int prismIndex, int *nodes);

Grid *gridAddPyramid(Grid *g, int n0, int n1, int n2, int n3, int n4);
Grid *gridPyramid(Grid *g, int pyramidIndex, int *nodes);

Grid *gridAddQuad(Grid *g, int n0, int n1, int n2, int n3, int faceId );
Grid *gridQuad(Grid *g, int quadIndex, int *nodes, int *faceId );

Grid *gridMap(Grid *g, int node, double *map);
#define gridMapPointerAllocated(grid) ((NULL != grid->map))
#define gridMapPointer(grid, node) (&grid->map[6*node])

Grid *gridSetMap(Grid *g, int node,
		 double m11, double m12, double m13,
		             double m22, double m23,
		                         double m33);

Grid *gridInterpolateMap2(Grid *g, int node0, int node1, double ratio,
			  int target);

#define gridCostFunction(grid) (grid->costFunction)
Grid *gridSetCostFunction(Grid *g, int costFunction);
#define gridCOST_FCN_MEAN_RATIO            (0)
#define gridCOST_FCN_ASPECT_RATIO          (1)
#define gridCOST_FCN_EDGE_LENGTH           (2)
#define gridCOST_FCN_JAC_SCALED_MEAN_RATIO (3)

#define gridCostConstraint(grid) (grid->costConstraint)
Grid *gridSetCostConstraint(Grid *g, int costConstraint);
#define gridCOST_CNST_VOLUME (0x01)
#define gridCOST_CNST_VALID  (0x02)
#define gridCOST_CNST_AREAUV (0x04)

Grid *gridSetMinInsertCost(Grid *g, double min_cost );
double gridMinInsertCost(Grid *g );

Grid *gridSetMinSurfaceSmoothCost(Grid *g, double min_cost );
double gridMinSurfaceSmoothCost(Grid *g );

int gridStoredCostDegree( Grid *g );
Grid *gridClearStoredCost( Grid *g );
Grid *gridStoreCost( Grid *g, double cost, double *costDerivative );
double gridStoredCost( Grid *g, int index );
Grid *gridStoredCostDerivative( Grid *g, int index, double *costDerivative );

#define gridLines(grid) (NULL==grid?NULL:grid->lines)
Grid *gridFreezeLinesNodes(Grid *g);
Grid *gridReportLinesLocation(Grid *g);

Grid *gridCopyAboutY0(Grid *g, int symmetryFaceId, int mirrorAux );

Grid *gridReportZeroDegreeNodes(Grid *g);

int gridPhase(Grid *g);
Grid *gridSetPhase(Grid *g, int phase);

#define gridALL_PHASE  (0)
#define gridEDGE_PHASE (1)
#define gridFACE_PHASE (2)
#define gridVOL_PHASE  (3)

Grid *gridCacheCurrentGridAndMap(Grid *g);


END_C_DECLORATION

#endif /* GRID_H */
