
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <values.h>
#include "master_header.h"
#include "layer.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "grid.h"

typedef struct Normal Normal;
struct Normal {
  int constrained;
  int root, tip;
  double direction[3];
  double height;
  bool terminated;
};

typedef struct Triangle Triangle;
struct Triangle {
  int normal[3];
  int constrainedSide[3];
  int parentGeomEdge[3];
};

typedef struct Blend Blend;
struct Blend {
  int nodes[2];
  int normal[4];
  int edgeId[2];
  int oldnormal[4];
};

struct Layer {
  Grid *grid;
  int maxtriangle, ntriangle;
  Triangle *triangle;
  int maxParentGeomFace, nParentGeomFace, *ParentGeomFace;
  int maxblend, nblend;
  Blend *blend;
  int maxnormal, nnormal;
  Normal *normal;
  int *globalNode2Normal;
  int nConstrainingGeometry, *constrainingGeometry;
  Adj *adj;
  bool mixedElementMode;

  bool *cellInLayer;
  bool *faceInLayer;
};

Layer *layerCreate( Grid *grid )
{
  int i;
  Layer *layer;
  layer = malloc(sizeof(Layer));
  layer->grid = grid;
  gridAttachNodeSorter( grid, layerSortGlobalNodes, layer );
  layer->maxtriangle=0;
  layer->ntriangle=0;
  layer->triangle=NULL;
  layer->maxParentGeomFace=0;
  layer->nParentGeomFace=0;
  layer->ParentGeomFace=NULL;
  layer->maxblend=0;
  layer->nblend=0;
  layer->blend=NULL;
  layer->maxnormal=0;
  layer->nnormal=0;
  layer->normal=NULL;
  layer->globalNode2Normal=NULL;
  layer->nConstrainingGeometry=0;
  layer->constrainingGeometry=NULL;
  layer->adj=NULL;
  layer->mixedElementMode=FALSE;

  layer->cellInLayer = malloc(gridMaxCell(grid)*sizeof(bool));
  for (i=0;i<gridMaxCell(grid);i++) layer->cellInLayer[i] = FALSE;
  layer->faceInLayer = malloc(gridMaxFace(grid)*sizeof(bool));
  for (i=0;i<gridMaxFace(grid);i++) layer->faceInLayer[i] = FALSE;
  return layer;
}

Grid *layerGrid(Layer *layer)
{
  return layer->grid;
}

void layerFree(Layer *layer)
{
  free(layer->faceInLayer);
  free(layer->cellInLayer);
  gridDetachNodeSorter( layer->grid );
  if (layer->adj != NULL) adjFree(layer->adj);
  if (layer->constrainingGeometry != NULL) free(layer->constrainingGeometry);
  if (layer->globalNode2Normal != NULL) free(layer->globalNode2Normal);
  if (layer->normal != NULL) free(layer->normal);
  if (layer->ParentGeomFace != NULL) free(layer->ParentGeomFace);
  if (layer->blend != NULL) free(layer->blend);
  if (layer->triangle != NULL) free(layer->triangle);
  free(layer);
}

Layer *formAdvancingFront( Grid *grid, char *project )
{
  Layer *layer;
  int i, nbc, bc[2];
  bool box, om6, n12;
  
  box = (NULL != strstr( project, "box"));
  om6 = (NULL != strstr( project, "om6"));
  n12 = (NULL != strstr( project, "n12"));

  if (box) printf("string %s has box.\n",project);
  if (om6) printf("string %s has om6.\n",project);
  if (n12) printf("string %s has n12.\n",project);

  bc[0]=1;
  bc[1]=2;
  if(box) nbc = 1;
  if(om6) nbc = 2;
  if(n12){
    nbc=2;
    bc[0]=5;
    bc[1]=6;
  }
  printf("make advancing layer object.\n");
  layer = layerCreate(grid);
  printf("make advancing layer triangle.\n");
  layerPopulateAdvancingFront(layer,nbc,bc);
  if (box) {
    layerConstrainNormal(layer,-9);
    layerConstrainNormal(layer,-10);
    layerConstrainNormal(layer,-11);
    layerConstrainNormal(layer,-12);
    layerConstrainNormal(layer,3);
    layerConstrainNormal(layer,4);
    layerConstrainNormal(layer,5);
    layerConstrainNormal(layer,6);
  }
  if (om6) {
    layerConstrainNormal(layer,5);
  }
  if(n12){
    layerConstrainNormal(layer,1);
    layerConstrainNormal(layer,2);
  }
  printf("make advancing layer triangle normals visible to triangle.\n");
  layerVisibleNormals(layer,-1.0,-1.0);
  return layer;
}

void layerSortGlobalNodes(void *voidLayer, int *o2n)
{
  Layer *layer = (Layer *)voidLayer;
  int i, triangle, normal;

  for (normal = 0 ; normal < layerNNormal(layer) ; normal++ ) {
    if (EMPTY != layer->normal[normal].root)
      layer->normal[normal].root = o2n[layer->normal[normal].root];
    if (EMPTY != layer->normal[normal].tip)
      layer->normal[normal].tip = o2n[layer->normal[normal].tip];
  }

}

int layerMaxTriangle(Layer *layer)
{
  return layer->ntriangle;
}

int layerNTriangle(Layer *layer)
{
  return layer->ntriangle;
}

int layerNBlend(Layer *layer)
{
  return layer->nblend;
}

int layerMaxNormal(Layer *layer)
{
  return layer->maxnormal;
}

int layerNNormal(Layer *layer)
{
  return layer->nnormal;
}

int layerMaxNode(Layer *layer)
{
  return gridMaxNode(layer->grid);
}

Layer *layerPopulateAdvancingFront(Layer *layer, int nbc, int *bc)
{
  int i, ibc, face, id, nodes[3];
  int triangle;
  int globalNodeId;
  int normal;
  double direction[3], *norm;
  double length;
  Grid *grid;

  grid = layerGrid(layer);

  for(ibc=0;ibc<nbc;ibc++) layerAddParentGeomFace(layer,bc[ibc]);

  for (ibc=0;ibc<nbc;ibc++){
    for(face=0;face<gridMaxFace(layer->grid);face++){
      if (grid == gridFace(grid,face,nodes,&id) &&
	  id==bc[ibc] ) {
	layerAddTriangle(layer,nodes[0],nodes[1],nodes[2]);
	layer->faceInLayer[face]=TRUE;
      }
    }
  }

  layerBuildNormalTriangleAdjacency(layer);

  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    for(i=0;i<3;i++){
      normal = layer->triangle[triangle].normal[i];
      if (layer != layerTriangleDirection(layer,triangle,direction))
	printf("Error: layerPopulateAdvancingFront: %s: %d: %s\n",
	       __FILE__,__LINE__,"triangle direction");
      layer->normal[normal].direction[0] += direction[0];
      layer->normal[normal].direction[1] += direction[1];
      layer->normal[normal].direction[2] += direction[2];
    }
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    norm = layer->normal[normal].direction;
    length = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    if (length > 0.0) {
      for ( i=0;i<3;i++) norm[i] = norm[i]/length;
    }else{
      for ( i=0;i<3;i++) norm[i] = 0.0;
    }
  }
  
  return layer;
}

Layer *layerBuildNormalTriangleAdjacency(Layer *layer)
{
  int triangle, i, normals[3];

  if (NULL != layer->adj) adjFree(layer->adj);
  layer->adj = adjCreate( layerNNormal(layer), layerNTriangle(layer)*3  );

  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangleNormals(layer,triangle,normals);
    for (i=0;i<3;i++) adjRegister( layer->adj, normals[i], triangle );
  }

  return layer;
}

Layer *layerAddParentGeomFace(Layer *layer, int faceId )
{
  if (layerParentGeomFace(layer, faceId )) return NULL;

  if (layer->nParentGeomFace >= layer->maxParentGeomFace) {
    layer->maxParentGeomFace += 100;
    if (layer->ParentGeomFace == NULL) {
      layer->ParentGeomFace = malloc(layer->maxParentGeomFace*sizeof(int));
    }else{
      layer->ParentGeomFace = realloc(layer->ParentGeomFace,
				      layer->maxParentGeomFace*sizeof(int));
    }
  }

  layer->ParentGeomFace[layer->nParentGeomFace] = faceId;
  layer->nParentGeomFace++;

  return layer;
}

bool layerParentGeomFace(Layer *layer, int faceId )
{
  int i;

  for (i=0;i<layer->nParentGeomFace;i++) 
    if (faceId == layer->ParentGeomFace[i]) return TRUE;

  return FALSE;
}

Layer *layerAddTriangle(Layer *layer, int n0, int n1, int n2 )
{
  int i;
  Grid *grid;
  grid = layerGrid(layer);

  if (layer->ntriangle >= layer->maxtriangle) {
    layer->maxtriangle += 5000;
    if (layer->triangle == NULL) {
      layer->triangle = malloc(layer->maxtriangle*sizeof(Triangle));
    }else{
      layer->triangle = realloc(layer->triangle,layer->maxtriangle*sizeof(Triangle));
    }
  }

  layer->triangle[layer->ntriangle].normal[0] = layerUniqueNormalId(layer,n0);
  layer->triangle[layer->ntriangle].normal[1] = layerUniqueNormalId(layer,n1);
  layer->triangle[layer->ntriangle].normal[2] = layerUniqueNormalId(layer,n2);
  for (i=0;i<3;i++){
    layer->triangle[layer->ntriangle].constrainedSide[i] = 0;
    layer->triangle[layer->ntriangle].parentGeomEdge[i] = 0;
  }
  gridFreezeNode(grid,n0);
  gridFreezeNode(grid,n1);
  gridFreezeNode(grid,n2);

  layer->ntriangle++;

  return layer;
}

int layerForceTriangle(Layer *layer, int normal0, int normal1, int normal2 )
{
  int i;
  Grid *grid;
  grid = layerGrid(layer);

  if (layer->ntriangle >= layer->maxtriangle) {
    layer->maxtriangle += 5000;
    if (layer->triangle == NULL) {
      layer->triangle = malloc(layer->maxtriangle*sizeof(Triangle));
    }else{
      layer->triangle = realloc(layer->triangle,layer->maxtriangle*sizeof(Triangle));
    }
  }

  layer->triangle[layer->ntriangle].normal[0] = normal0;
  layer->triangle[layer->ntriangle].normal[1] = normal1;
  layer->triangle[layer->ntriangle].normal[2] = normal2;
  for (i=0;i<3;i++){
    layer->triangle[layer->ntriangle].constrainedSide[i] = 0;
    layer->triangle[layer->ntriangle].parentGeomEdge[i] = 0;
  }

  layer->ntriangle++;

  return (layer->ntriangle-1);
}

Layer *layerTriangle(Layer *layer, int triangle, int *nodes )
{

  if (triangle < 0 || triangle >= layerNTriangle(layer)) return NULL;
  nodes[0] = layerNormalRoot(layer,layer->triangle[triangle].normal[0]);
  nodes[1] = layerNormalRoot(layer,layer->triangle[triangle].normal[1]);
  nodes[2] = layerNormalRoot(layer,layer->triangle[triangle].normal[2]);
  
  return layer;
}

Layer *layerTriangleDirection(Layer *layer, int triangle, double *direction )
{
  int i, nodes[3];
  double node0[3], node1[3], node2[3];
  double edge1[3], edge2[3], norm[3], length; 
  
  if ( layer != layerTriangle(layer,triangle,nodes) ) return NULL;

  if (layer->grid != gridNodeXYZ( layer->grid, nodes[0], node0 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[1], node1 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[2], node2 )) return NULL;

  for (i = 0 ; i < 3 ; i++ ){
    edge1[i] = node1[i] - node0[i];
    edge2[i] = node2[i] - node0[i];
  }

  norm[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1]; 
  norm[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2]; 
  norm[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0]; 
  length = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);

  if (length > 0.0) {
    for ( i=0;i<3;i++) direction[i] = norm[i]/length;
  }else{
    for ( i=0;i<3;i++) direction[i] = 0.0;
  }

  return layer;
}

Layer *layerTriangleArea(Layer *layer, int triangle, double *area )
{
  int i, nodes[3];
  double node0[3], node1[3], node2[3];
  double edge1[3], edge2[3], norm[3], length; 
  
  if ( layer != layerTriangle(layer,triangle,nodes) ) return NULL;

  if (layer->grid != gridNodeXYZ( layer->grid, nodes[0], node0 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[1], node1 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[2], node2 )) return NULL;

  for (i = 0 ; i < 3 ; i++ ){
    edge1[i] = node1[i] - node0[i];
    edge2[i] = node2[i] - node0[i];
  }

  norm[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1]; 
  norm[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2]; 
  norm[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0]; 
  length = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);

  *area = length * 0.5;

  return layer;
}

Layer *layerTriangleCenter(Layer *layer, int triangle, double *center )
{
  int i, nodes[3];
  double node0[3], node1[3], node2[3];
  
  if ( layer != layerTriangle(layer,triangle,nodes) ) return NULL;

  if (layer->grid != gridNodeXYZ( layer->grid, nodes[0], node0 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[1], node1 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[2], node2 )) return NULL;

  for (i = 0 ; i < 3 ; i++ ) center[i] = (node0[i]+node1[i]+node2[i])/3.0;

  return layer;
}

Layer *layerTriangleMaxEdgeLength(Layer *layer, int triangle, double *length )
{
  int i, nodes[3];
  double node0[3], node1[3], node2[3];
  double edge0[3], edge1[3], edge2[3], maxLength; 
  
  if ( layer != layerTriangle(layer,triangle,nodes) ) return NULL;

  if (layer->grid != gridNodeXYZ( layer->grid, nodes[0], node0 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[1], node1 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[2], node2 )) return NULL;

  gridSubtractVector(node1,node0,edge0);
  gridSubtractVector(node2,node1,edge1);
  gridSubtractVector(node0,node2,edge2);

  maxLength = -1.0;
  maxLength = MAX(maxLength,sqrt(gridDotProduct(edge0,edge0)));
  maxLength = MAX(maxLength,sqrt(gridDotProduct(edge1,edge1)));
  maxLength = MAX(maxLength,sqrt(gridDotProduct(edge2,edge2)));
  
  *length = maxLength;

  return layer;
}

Layer *layerInitializeNormal(Layer *layer, int normal)
{
  if (normal < 0 || normal >= layerMaxNormal(layer) ) return NULL;
 
  layer->normal[normal].constrained = 0;
  layer->normal[normal].root = EMPTY;
  layer->normal[normal].tip = EMPTY;
  layer->normal[normal].direction[0] = 0.0;
  layer->normal[normal].direction[1] = 0.0;
  layer->normal[normal].direction[2] = 0.0;
  layer->normal[normal].height = 1.0;
  layer->normal[normal].terminated = FALSE;

  return layer;
}

int layerDuplicateNormal(Layer *layer, int normal)
{
  int newone, root;

  root = layerNormalRoot(layer,normal);
  if (EMPTY == root) return EMPTY;
 
  newone = layerAddNormal(layer, root);
  if (EMPTY == newone) return EMPTY;

  layer->normal[newone].constrained =  layer->normal[normal].constrained;
  layer->normal[newone].root =         layer->normal[normal].root;
  layer->normal[newone].tip =          layer->normal[normal].tip;
  layer->normal[newone].direction[0] = layer->normal[normal].direction[0];
  layer->normal[newone].direction[1] = layer->normal[normal].direction[1];
  layer->normal[newone].direction[2] = layer->normal[normal].direction[2];
  layer->normal[newone].height =       layer->normal[normal].height;
  layer->normal[newone].terminated =   layer->normal[normal].terminated;

  return newone;
}

int layerAddNormal(Layer *layer, int globalNodeId )
{
  int i; 

  if (globalNodeId < 0 || globalNodeId >= layerMaxNode(layer) ) return EMPTY;

  if (layer->nnormal >= layer->maxnormal) {
    layer->maxnormal += 5000;
    if (layer->normal == NULL) {
      layer->normal = malloc(layer->maxnormal*sizeof(Normal));
      layer->globalNode2Normal = malloc(layerMaxNode(layer)*sizeof(int));
      for (i=0;i<layerMaxNode(layer);i++) layer->globalNode2Normal[i]=EMPTY;
    }else{
      layer->normal = realloc(layer->normal,layer->maxnormal*sizeof(Normal));
    }
  }

  if (layer != layerInitializeNormal(layer, layer->nnormal)) return EMPTY;
  layer->normal[layer->nnormal].root = globalNodeId;
  layer->globalNode2Normal[globalNodeId] = layer->nnormal;

  layer->nnormal++;

  return (layer->nnormal - 1);
}

int layerUniqueNormalId(Layer *layer, int globalNodeId)
{
  if ( layer->globalNode2Normal == NULL || 
       layer->globalNode2Normal[globalNodeId] == EMPTY ) 
    return layerAddNormal(layer,globalNodeId);

  return layer->globalNode2Normal[globalNodeId];
}

Layer *layerTriangleNormals(Layer *layer, int triangle, int *normals )
{
  if (triangle < 0 || triangle >= layerNTriangle(layer)) return NULL;
  normals[0] = layer->triangle[triangle].normal[0];
  normals[1] = layer->triangle[triangle].normal[1];
  normals[2] = layer->triangle[triangle].normal[2];
  
  return layer;
}

int layerNormalRoot(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return EMPTY;
  return layer->normal[normal].root;
}

int layerNormalDeg(Layer *layer, int normal )
{
  if (NULL == layer->adj) return 0;
  return adjDegree(layer->adj, normal);
}

Layer *layerNormalTriangles(Layer *layer, int normal, int ntriangle, int *triangles )
{
  int i;
  AdjIterator it;
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;
  i=0;
  for ( it = adjFirst(layer->adj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
    if (i>=ntriangle) return NULL;
    triangles[i] = adjItem(it);
    i++;
  }
  return layer;
}

int layerPreviousTriangle(Layer *layer, int normal, int triangle )
{
  int nodes[3], previousRoot, previousTriangle;
  int root;
  AdjIterator it;

  if (layer != layerTriangle(layer, triangle, nodes ) ) return EMPTY;
  root = layerNormalRoot(layer,normal);

  previousRoot = EMPTY;
  if (root == nodes[0]) previousRoot = nodes[1];
  if (root == nodes[1]) previousRoot = nodes[2];
  if (root == nodes[2]) previousRoot = nodes[0];
  if (EMPTY == previousRoot) return EMPTY;

  for ( it = adjFirst(layer->adj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
    previousTriangle = adjItem(it);
    if ( triangle != previousTriangle ) {
      layerTriangle(layer, previousTriangle, nodes );
      if ( previousRoot == nodes[0] ||
	   previousRoot == nodes[1] ||
	   previousRoot == nodes[2] ) return previousTriangle;
    }
  }

  return EMPTY;
}

int layerNextTriangle(Layer *layer, int normal, int triangle )
{
  int nodes[3], nextRoot, nextTriangle;
  int root;
  AdjIterator it;

  if (layer != layerTriangle(layer, triangle, nodes ) ) return EMPTY;
  root = layerNormalRoot(layer,normal);

  nextRoot = EMPTY;
  if (root == nodes[0]) nextRoot = nodes[2];
  if (root == nodes[1]) nextRoot = nodes[0];
  if (root == nodes[2]) nextRoot = nodes[1];
  if (EMPTY == nextRoot) return EMPTY;

  for ( it = adjFirst(layer->adj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
    nextTriangle = adjItem(it);
    if ( triangle != nextTriangle ) {
      layerTriangle(layer, nextTriangle, nodes );
      if ( nextRoot == nodes[0] ||
	   nextRoot == nodes[1] ||
	   nextRoot == nodes[2] ) return nextTriangle;
    }
  }

  return EMPTY;
}

Layer *layerCommonEdge(Layer *layer, int triangle0, int triangle1, int *nodes)
{
  int nodes0[3], nodes1[3];
  int n0, n1, g0, g1, start, end;

  if (layer != layerTriangle(layer,triangle0,nodes0) ) return NULL;
  if (layer != layerTriangle(layer,triangle1,nodes1) ) return NULL;

  start = end = EMPTY;
  for ( n0 = 0 ; n0 < 3 ; n0++) {
    n1 = n0+1; if ( n1 >= 3 ) n1 = 0;
    g0 = nodes0[n0];
    g1 = nodes0[n1];
    if ( ( g0 == nodes1[0] || g0 == nodes1[1] || g0 == nodes1[2] ) &&
	 ( g1 == nodes1[0] || g1 == nodes1[1] || g1 == nodes1[2] ) ) {
      start = g0;
      end   = g1;
    } 
  }

  if ( start == EMPTY || end == EMPTY ) return NULL;

  nodes[0]=start;
  nodes[1]=end;

  return layer;
}

#define PI (3.14159265358979)
#define ConvertRadianToDegree(radian) ((radian)*57.2957795130823)
#define ConvertDegreeToRadian(degree) ((degree)*0.0174532925199433)

double layerEdgeAngle(Layer *layer, int triangle0, int triangle1 )
{
  double direction0[3], direction1[3];
  double dot, radian;
  int commonEdge[2];
  double xyzstart[3], xyzend[3], edge[3], cross[3];

  if ( triangle0 == triangle1 ) return -3.0;
  if ( layer != layerTriangleDirection(layer,triangle0,direction0))return -1.0;
  if ( layer != layerTriangleDirection(layer,triangle1,direction1))return -2.0;
  if ( layer != layerCommonEdge(layer,triangle0,triangle1,commonEdge) ) 
    return -4.0;

  dot = gridDotProduct(direction0, direction1);
  radian = acos(dot);  

  gridNodeXYZ( layerGrid(layer), commonEdge[0], xyzstart );
  gridNodeXYZ( layerGrid(layer), commonEdge[1], xyzend );

  gridSubtractVector( xyzstart, xyzend, edge );
  gridCrossProduct( direction0, direction1, cross);
  
  if ( gridDotProduct( edge, cross ) < 0.0 ) {
    radian += PI;
  }else{
    radian = PI - radian;
  }
  return ConvertRadianToDegree(radian);
}

Layer *layerNormalDirection(Layer *layer, int normal, double *direction )
{
  int i;
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  for ( i=0;i<3;i++) direction[i] = layer->normal[normal].direction[i];

  return layer;
}
Layer *layerAssignPolynomialNormalHeight(Layer *layer, 
                                         double constant, double slope,
					 double exponent,
                                         double *origin, double *direction)
{
  int normal;
  double distance;
  double distanceVector[3], normalOrigin[3];
  
  for(normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(layerGrid(layer), layerNormalRoot(layer,normal), normalOrigin);
    gridSubtractVector(normalOrigin,origin,distanceVector);
    distance = gridDotProduct(distanceVector,direction);
    if (distance >= 0.0 ) 
      layerSetNormalHeight( layer, normal, constant + slope*pow(distance,exponent) );
  }
  return layer;	
}


Layer *layerSetHeightOfAllNormals(Layer *layer, double height )
{
  int normal;

  for(normal=0;normal<layerNNormal(layer);normal++)
    layerSetNormalHeight( layer, normal, height );

  return layer;
}

Layer *layerLaminarInitialHeight(Layer *layer, double Re, double xStart)
{
  int normal;
  double xyz[3];
  double totalHeight, initialHeight;

  for(normal=0;normal<layerNNormal(layer); normal++){
    gridNodeXYZ(layerGrid(layer),layerNormalRoot(layer,normal),xyz);
    totalHeight = 5.2 * sqrt(ABS(xyz[0]-xStart)) / sqrt(ABS(Re));
    totalHeight = MIN(ABS(xyz[0]-xStart),totalHeight);
    initialHeight = totalHeight / 50.0;
    layerSetNormalHeight(layer,normal,initialHeight);
  }

  return layer;
}

Layer *layerLaminarInitialHeightNegZ(Layer *layer)
{
  int normal;
  double xyz[3];
  double initialHeight;

  for(normal=0;normal<layerNNormal(layer); normal++){
    gridNodeXYZ(layerGrid(layer),layerNormalRoot(layer,normal),xyz);
    initialHeight = 0.0004 - 1.2e-6 * xyz[2];
    if (xyz[0]>20.0) initialHeight += (xyz[0]-20.0)*0.00005;
    layerSetNormalHeight(layer,normal,initialHeight);
  }

  return layer;
}

Layer *layerSetNormalHeightOfFace(Layer *layer, int faceId, double height )
{
  int face, nodes[3], id;
  int i, normal;

  Grid *grid; grid = layerGrid(layer);

  for (face=0;face<gridMaxFace(grid);face++){
    if (grid == gridFace(grid,face,nodes,&id) && id == faceId ){
      for (i=0;i<3;i++){
	normal = layer->globalNode2Normal[nodes[i]];
	layerSetNormalHeight( layer, normal, height );
      }
    }
  }

  return layer;
}

Layer *layerSetNormalHeight(Layer *layer, int normal, double height)
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  layer->normal[normal].height=height;
  if (height<=0.0) layerTerminateNormal(layer,normal);

  return layer;
}

Layer *layerGetNormalHeight(Layer *layer, int normal, double *height)
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  *height = layer->normal[normal].height;

  return layer;
}

Layer *layerScaleNormalHeight(Layer *layer, double scale)
{
  int normal;

  for(normal=0;normal<layerNNormal(layer); normal++)
    layer->normal[normal].height=scale*layer->normal[normal].height;

  return layer;
}

Layer *layerNormalMinDot(Layer *layer, int normal,
			 double *mindot, double *mindir,
			 int *minTriangle )
{
  int triangle;
  double dir[3], norm[3], dot;
  AdjIterator it;

  if (layer != layerNormalDirection(layer,normal,dir)) return NULL;

  *mindot = 2.0;
  mindir[0]=0.0;
  mindir[1]=0.0;
  mindir[2]=0.0;
  *minTriangle = EMPTY;

  for ( it = adjFirst(layer->adj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
    triangle = adjItem(it);
    layerTriangleDirection(layer,triangle,norm);
    dot = norm[0]*dir[0] + norm[1]*dir[1] + norm[2]*dir[2];
    if (dot<*mindot) {
      *mindot = dot;
      mindir[0]=norm[0];
      mindir[1]=norm[1];
      mindir[2]=norm[2];
      *minTriangle = triangle;
    }
  }
  
  return layer;
}

Layer *layerVisibleNormals(Layer *layer, double dotLimit, double radianLimit )
{
  int normal, iter, i;
  double *dir, mindir[3], dot, mindot, radian, length; 
  int minTriangle, lastTriangle;

  if (layerNNormal(layer) == 0 ) return NULL;
  if (dotLimit < 0) dotLimit = 0.90;
  if (radianLimit < 0) radianLimit = 1.0e-10;

  layerProjectNormalsToConstraints(layer);

  for (normal=0;normal<layerNNormal(layer);normal++){
    if ( 0 != layerNormalDeg(layer, normal ) ) {
      dir = layer->normal[normal].direction;
      lastTriangle = EMPTY;
      radian = 0.01;
      layerNormalMinDot(layer, normal, &mindot, mindir, &minTriangle );
      for (iter=0;iter<1000 && radian > radianLimit && mindot < dotLimit;iter++){
	if (minTriangle != lastTriangle) {
	  radian = radian * 0.5;
	  lastTriangle = minTriangle;
	  //printf("normal %d, dot %f rad %e triangle %d\n",normal,mindot,radian,minTriangle);
	}
	dir[0] += radian*mindir[0];
	dir[1] += radian*mindir[1];
	dir[2] += radian*mindir[2];
	length = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
	if (length > 0.0) {
	  for ( i=0;i<3;i++) dir[i] = dir[i]/length;
	}else{
	  for ( i=0;i<3;i++) dir[i] = 0.0;
	}
	layerNormalMinDot(layer, normal, &mindot, mindir, &minTriangle );
      }
      if (mindot <= 0.0 ) {
	double xyz[3];
	gridNodeXYZ(layerGrid(layer),layerNormalRoot(layer,normal),xyz);
	printf("ERROR: %s, %d, Invisible norm %d dot%9.5f X%9.5f Y%9.5f Z%9.5f\n",
	       __FILE__, __LINE__, normal, mindot,xyz[0],xyz[1],xyz[2]);
      }
    }
  }

  return layer;
}

Layer *layerSmoothNormalDirection(Layer *layer)
{
  int normal, iter, triangle, normals[3], total, i;
  double norm[3], avgdir[3], denom; 
  AdjIterator it;
  int minTriangle, lastTriangle;

  if (layerNNormal(layer) == 0 ) return NULL;
  if (layerNBlend(layer) != 0 ) return NULL;

  layerProjectNormalsToConstraints(layer);

  for (iter=0;iter<1;iter++){
    for (normal=0;normal<layerNNormal(layer);normal++){
      if ( 0 < layerConstrained(layer,normal) ){
	total = 0;
	avgdir[0]=0.0;
	avgdir[1]=0.0;
	avgdir[2]=0.0;
	for ( it = adjFirst(layer->adj,normal); 
	      adjValid(it); 
	      it = adjNext(it) ){
	  triangle = adjItem(it);
	  layerTriangleNormals(layer,triangle,normals);
	  for (i=0;i<3;i++){
	    if (normal != normals[i] && 0 != layerConstrained(layer,normal) ){
	      layerNormalDirection(layer,normals[i],norm);
	      avgdir[0] += norm[0];
	      avgdir[1] += norm[1];
	      avgdir[2] += norm[2];
	      total++;
	    }
	  }
	}
	denom = 1.0 / (double)total;
	layer->normal[normal].direction[0] = avgdir[0] * denom;
	layer->normal[normal].direction[1] = avgdir[1] * denom;
	layer->normal[normal].direction[2] = avgdir[2] * denom;
      }
    }
    layerVisibleNormals(layer,0.5,1.0e-5);
    for (normal=0;normal<layerNNormal(layer);normal++){
      if ( 0 == layerConstrained(layer,normal) ){
	total = 0;
	avgdir[0]=0.0;
	avgdir[1]=0.0;
	avgdir[2]=0.0;
	for ( it = adjFirst(layer->adj,normal); 
	      adjValid(it); 
	      it = adjNext(it) ){
	  triangle = adjItem(it);
	  layerTriangleNormals(layer,triangle,normals);
	  for (i=0;i<3;i++){
	    if (normal != normals[i]){
	      layerNormalDirection(layer,normals[i],norm);
	      avgdir[0] += norm[0];
	      avgdir[1] += norm[1];
	      avgdir[2] += norm[2];
	      total++;
	    }
	  }
	}
	denom = 1.0 / (double)total;
	layer->normal[normal].direction[0] = avgdir[0] * denom;
	layer->normal[normal].direction[1] = avgdir[1] * denom;
	layer->normal[normal].direction[2] = avgdir[2] * denom;
      }
    }
    layerVisibleNormals(layer,0.5,1.0e-5);
  }

  return layer;
}

Layer *layerProjectNormalsToConstraints(Layer *layer)
{
  int normal, tipnode, edgeId, faceId;
  double xyzroot[3], xyztip[3], direction[3], height;
  Grid *grid;

  grid = layerGrid(layer);

  tipnode = gridAddNode(grid,0,0,0);

  if (EMPTY == tipnode) return NULL;

  for (normal=0;normal<layerNNormal(layer);normal++){
    faceId = layerConstrained(layer,normal);
    if (faceId != 0 ) {
      gridNodeXYZ(grid, layerNormalRoot(layer,normal), xyzroot );
      layerNormalDirection(layer,normal,direction);
      layerGetNormalHeight(layer,normal,&height);
      direction[0] *= height;
      direction[1] *= height;
      direction[2] *= height;
      xyztip[0] = xyzroot[0] + direction[0];
      xyztip[1] = xyzroot[1] + direction[1];
      xyztip[2] = xyzroot[2] + direction[2];
      gridSetNodeXYZ(grid,tipnode,xyztip);
      if (faceId >0) {
	gridForceNodeToFace(grid, tipnode, faceId );
      }else{
	edgeId = -faceId;
	gridForceNodeToEdge(grid, tipnode, edgeId );
      }
      gridNodeXYZ(grid, tipnode, xyztip);
      gridSubtractVector(xyztip,xyzroot,direction);
      height = sqrt(gridDotProduct(direction,direction));
      layer->normal[normal].direction[0] = direction[0] / height;
      layer->normal[normal].direction[1] = direction[1] / height;
      layer->normal[normal].direction[2] = direction[2] / height;
    }
    
  }

  gridRemoveNode(grid,tipnode);

  return layer;
}

Layer *layerConstrainNormal(Layer *layer, int edgeface )
{
  int edgeId, faceId, face, nodes[3], id, i, normal;
  int n0, n1, normal0, normal1;
  int nCurveNode, *curve;
  int exist;

  if (layerNNormal(layer) == 0 ) return NULL;
  if (edgeface == 0 ) return NULL;
  
  if (edgeface > 0) {
    faceId = edgeface;
    for(face=0;face<gridMaxFace(layer->grid);face++){
      if (layer->grid == gridFace(layer->grid,face,nodes,&id) && id==faceId ){
	for(i=0;i<3;i++){
	  normal = layer->globalNode2Normal[nodes[i]];
	  if (normal != EMPTY) {
	    if ( layer->normal[normal].constrained >= 0) {
	      layer->normal[normal].constrained=faceId;
	    }
	  }
	}
	for(n0=0;n0<3;n0++){
	  n1 = n0 + 1; if ( n1 > 2 ) n1 = 0;
	  normal0 = layer->globalNode2Normal[nodes[n0]];
	  normal1 = layer->globalNode2Normal[nodes[n1]];
	  layerConstrainTriangleSide( layer, normal0, normal1, faceId );
	}
      }
    }
  }else{
    edgeId = -edgeface;
    nCurveNode = gridGeomEdgeSize( layer->grid, edgeId );
    curve = malloc( nCurveNode * sizeof(int) );
    gridGeomEdge( layer->grid, edgeId, curve );
    for ( i=0; i<nCurveNode; i++){ 
      normal = layer->globalNode2Normal[curve[i]];
      if (normal != EMPTY) layer->normal[normal].constrained=-edgeId;
    }
    free(curve);
  }
  
  if ( !layerConstrainingGeometry(layer, edgeface) ) {
    if ( layer->constrainingGeometry == NULL ) {
      layer->nConstrainingGeometry = 1;
      layer->constrainingGeometry = 
	malloc( layer->nConstrainingGeometry * sizeof(int) );
      layer->constrainingGeometry[0] = edgeface;
    }else{
      layer->nConstrainingGeometry++;
      layer->constrainingGeometry = 
	realloc( layer->constrainingGeometry, 
		 layer->nConstrainingGeometry * sizeof(int) );
      layer->constrainingGeometry[layer->nConstrainingGeometry-1] = edgeface;
    }
  }

  return layer;
}

bool layerConstrainingGeometry(Layer *layer, int edgeface )
{
  int i;

  if ( edgeface == 0 ) return FALSE;

  for (i=0;i<layer->nConstrainingGeometry;i++)
    if (layer->constrainingGeometry[i]==edgeface) return TRUE;
  
  return FALSE;
}

int layerConstrained(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return 0;
  return layer->normal[normal].constrained;
}

Layer *layerConstrainTriangleSide(Layer *layer, int normal0, int normal1, int bc )
{
  AdjIterator it;
  int i0, i1, triangle, side;
 
  if (normal0 < 0 || normal0 >= layerNNormal(layer) ) return NULL;
  if (normal1 < 0 || normal1 >= layerNNormal(layer) ) return NULL;

  for ( it = adjFirst(layer->adj,normal0); 
	adjValid(it); 
	it = adjNext(it) ){
    triangle = adjItem(it);
    for (i1=0;i1<3;i1++) {
      if (normal1 == layer->triangle[triangle].normal[i1]) {
	for (i0=0;i0<3;i0++) {
	  if (normal0 == layer->triangle[triangle].normal[i0]) {
	    side = MIN(i0,i1);
	    if ( side == 0 && 2 == MAX(i0,i1) ) side = 2;
	    layer->triangle[triangle].constrainedSide[side]=bc;
	  }
	}
      }   
    }
  }
  return layer;
}

int layerConstrainedSide(Layer *layer, int triangle, int side )
{
  if (triangle < 0 || triangle >= layerNTriangle(layer) ) return 0;

  if (side < 0 || side > 2 ) return 0;

  return layer->triangle[triangle].constrainedSide[side];
}

int layerNConstrainedSides(Layer *layer, int faceId )
{
  int triangle, i, nside;

  if (faceId==0) return 0;
  nside = 0;
  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    for(i=0;i<3;i++){
      if (layerConstrainedSide(layer,triangle,i)==faceId) nside++;
    }
  }
  return nside;
}

Layer *layerFindParentGeomEdges(Layer *layer)
{
  int triangle, nodes[3], side, n0, n1, edgeId;

  if (layerNTriangle(layer) == 0 ) return NULL;
  if (layerNNormal(layer) == 0 ) return NULL;

  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangle(layer,triangle,nodes);
    for(side=0;side<3;side++){
      n0 = side;
      n1 = side+1; if (n1>2) n1 = 0;
      n0 = nodes[n0];
      n1 = nodes[n1];
      edgeId = gridEdgeId(layer->grid,n0,n1);
      if (EMPTY != edgeId) layer->triangle[triangle].parentGeomEdge[side]=edgeId;
    }
  }
  return layer;
}

Layer *layerSetParentGeomEdge(Layer *layer, int normal0, int normal1, int edgeId )
{
  AdjIterator it;
  int i0, i1, triangle, side;
 
  if (normal0 < 0 || normal0 >= layerNNormal(layer) ) return NULL;
  if (normal1 < 0 || normal1 >= layerNNormal(layer) ) return NULL;

  for ( it = adjFirst(layer->adj,normal0); 
	adjValid(it); 
	it = adjNext(it) ){
    triangle = adjItem(it);
    for (i1=0;i1<3;i1++) {
      if (normal1 == layer->triangle[triangle].normal[i1]) {
	for (i0=0;i0<3;i0++) {
	  if (normal0 == layer->triangle[triangle].normal[i0]) {
	    side = MIN(i0,i1);
	    if ( side == 0 && 2 == MAX(i0,i1) ) side = 2;
	    layer->triangle[triangle].parentGeomEdge[side]=edgeId;
	  }
	}
      }   
    }
  }
  return layer;
}

int layerParentGeomEdge(Layer *layer, int triangle, int side )
{
  if (triangle < 0 || triangle >= layerNTriangle(layer) ) return 0;

  if (side < 0 || side > 2 ) return 0;

  return layer->triangle[triangle].parentGeomEdge[side];
}

int layerNParentGeomEdgeSegments(Layer *layer, int edgeId )
{
  int triangle, i, nSegments;

  if (edgeId==0) return 0;
  nSegments = 0;
  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    for(i=0;i<3;i++){
      if (layerParentGeomEdge(layer,triangle,i)==edgeId) nSegments++;
    }
  }
  return nSegments;
}

Layer *layerTerminateNormal(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  layer->normal[normal].terminated = TRUE;

  return layer;
}

bool layerNormalTerminated(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return FALSE;

  return layer->normal[normal].terminated;
}

Layer *layerTerminateFaceNormals(Layer *layer, int faceId)
{
  int face, nodes[3], id;
  int i, normal;

  Grid *grid; grid = layerGrid(layer);

  for (face=0;face<gridMaxFace(grid);face++){
    if (grid == gridFace(grid,face,nodes,&id) && id == faceId ){
      for (i=0;i<3;i++){
	normal = layer->globalNode2Normal[nodes[i]];
	layerTerminateNormal( layer, normal );
      }
    }
  }

  return layer;
}


int layerNActiveNormal(Layer *layer )
{
  int normal, nActive;

  nActive=0;

  for ( normal=0 ; normal < layerNNormal(layer); normal++ ) 
    if ( !layer->normal[normal].terminated ) nActive++;

  return nActive;
}

bool layerCellInLayer(Layer *layer, int cell)
{
  if (cell < 0 || cell >= gridMaxCell(layerGrid(layer)) ) return FALSE;
  return layer->cellInLayer[cell];
}

bool layerFaceInLayer(Layer *layer, int face)
{
  if (face < 0 || face >= gridMaxFace(layerGrid(layer)) ) return FALSE;
  return layer->faceInLayer[face];
}

bool layerEdgeInLayer(Layer *layer, int edge)
{
  int nodes[2], edgeId;
  Grid *grid;
  grid = layerGrid(layer);

  if ( grid != gridEdge(grid, edge, nodes, &edgeId) ) return FALSE;

  return (gridNodeFrozen(grid,nodes[0]) && gridNodeFrozen(grid,nodes[1]) );
}

Layer *layerReconnectCellUnlessInLayer(Layer *layer, int oldNode, int newNode )
{
  AdjIterator it;
  int cell, nodes[4], i;
  Grid *grid;

  if (oldNode < 0 || oldNode >= layerMaxNode(layer) ) return NULL;
  if (newNode < 0 || newNode >= layerMaxNode(layer) ) return NULL;
  if (newNode == oldNode) return layer;
  
  grid = layerGrid(layer);

  it = adjFirst(gridCellAdj(grid),oldNode);
  while (adjValid(it)){
    cell = adjItem(it);
    if (!layerCellInLayer(layer,cell) ) {
      gridCell(grid, cell, nodes);
      gridRemoveCell(grid,cell);
      for (i=0;i<4;i++) if (oldNode == nodes[i]) nodes[i] = newNode;
      gridAddCell(grid, nodes[0],  nodes[1],  nodes[2],  nodes[3]);
      it = adjFirst(gridCellAdj(grid),oldNode);
    }else{
      it = adjNext(it);
    }      
  }
  
  return layer;
}

Layer *layerReconnectFaceUnlessInLayer(Layer *layer, int faceId,
				       int oldNode, int newNode )
{
  AdjIterator it;
  int face, nodes[3], id;
  int i;
  double uv[6];
  Grid *grid;

  if (oldNode < 0 || oldNode >= layerMaxNode(layer) ) return NULL;
  if (newNode < 0 || newNode >= layerMaxNode(layer) ) return NULL;
  if (newNode == oldNode) return layer;
  
  grid = layerGrid(layer);

  it = adjFirst(gridFaceAdj(grid),oldNode);
  while (adjValid(it)){
    face = adjItem(it);
    gridFace(grid, face, nodes, &id);
    if (id == faceId && !layerFaceInLayer(layer,face) ) {
      gridNodeUV(grid,nodes[0],faceId,&uv[0]);
      gridNodeUV(grid,nodes[1],faceId,&uv[2]);
      gridNodeUV(grid,nodes[2],faceId,&uv[4]);
      gridRemoveFace(grid,face);
      for (i=0;i<3;i++) if (oldNode == nodes[i]) nodes[i] = newNode;
      gridAddFaceUV( grid, 
		     nodes[0], uv[0], uv[1],
		     nodes[1], uv[2], uv[3],
		     nodes[2], uv[4], uv[5],
		     faceId );
      it = adjFirst(gridFaceAdj(grid),oldNode);
    }else{
      it = adjNext(it);
    }      
  }
  
  return layer;
}

Layer *layerReconnectEdgeUnlessInLayer(Layer *layer, int edgeId,
				       int oldNode, int newNode )
{
  AdjIterator it;
  int edge, nodes[2], id;
  int i;
  double t[2];
  Grid *grid;

  if (oldNode < 0 || oldNode >= layerMaxNode(layer) ) return NULL;
  if (newNode < 0 || newNode >= layerMaxNode(layer) ) return NULL;
  if (newNode == oldNode) return layer;
  
  grid = layerGrid(layer);

  it = adjFirst(gridEdgeAdj(grid),oldNode);
  while (adjValid(it)){
    edge = adjItem(it);
    gridEdge(grid, edge, nodes, &id);
    if (id == edgeId && !layerEdgeInLayer(layer,edge) ) {
      gridNodeT(grid,nodes[0],edgeId,&t[0]);
      gridNodeT(grid,nodes[1],edgeId,&t[1]);
      gridRemoveEdge(grid,edge);
      for (i=0;i<2;i++) if (oldNode == nodes[i]) nodes[i] = newNode;
      gridAddEdge(grid, nodes[0], nodes[1], edgeId, t[0], t[1]);
      it = adjFirst(gridEdgeAdj(grid),oldNode);
    }else{
      it = adjNext(it);
    }      
  }
  
  return layer;
}

Layer *layerAdvanceConstantHeight(Layer *layer, double height )
{
  layerSetHeightOfAllNormals(layer, height );
  return layerAdvance(layer);
}

Layer *layerAdvance(Layer *layer)
{
  Grid *grid = layer->grid;
  int normal, normal0, normal1, root, tip, faceId, edgeId, i;
  int cell, node;
  int triangle, normals[3], nodes[3], n[6];
  int side[2], sidenode[2], sidenormal[2];
  double xyz[3];
  int nterminated;

  if (layerNNormal(layer) == 0 ) return NULL;

  for (normal=0;normal<layerNNormal(layer);normal++){
    root = layer->normal[normal].root;
    gridFreezeNode( grid, root );
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    if (layerNormalTerminated(layer,normal)){
      layer->normal[normal].tip = layer->normal[normal].root;
    }else{
      root = layer->normal[normal].root;
      gridNodeXYZ(grid,root,xyz);
      for(i=0;i<3;i++)
	xyz[i] = xyz[i] 
	  + layer->normal[normal].height
          * layer->normal[normal].direction[i];
      tip = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);
      if ( EMPTY == tip) return NULL;
      layer->normal[normal].tip = tip;
      layerReconnectCellUnlessInLayer(layer, root, tip);
      faceId = layerConstrained(layer,normal);
      if (0 > faceId) {
	edgeId = -faceId;
	layerReconnectEdgeUnlessInLayer(layer, edgeId, root, tip);
      }
      gridCopySpacing(grid, root, tip );
      gridFreezeNode( grid, tip );
    }
  }

  /* reconnect faces for constrained triangleside */
  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    for (i=0;i<3;i++){
      faceId = layerConstrainedSide(layer, triangle, i);
      if (faceId > 0) {
	normal0 = i;
	normal1 = i+1; if (normal1>2) normal1 = 0;
	normal0 = layer->triangle[triangle].normal[normal0];
	normal1 = layer->triangle[triangle].normal[normal1];
	layerReconnectFaceUnlessInLayer(layer, faceId, 
					layer->normal[normal0].root, 
					layer->normal[normal0].tip);
	layerReconnectFaceUnlessInLayer(layer, faceId, 
					layer->normal[normal1].root, 
					layer->normal[normal1].tip);
      }
    }    
  }

  /* advance edges */
  for (normal=0;normal<layerNNormal(layer);normal++) {
    edgeId = -layerConstrained(layer,normal);
    if (edgeId > 0) {
      root = layer->normal[normal].root;
      tip  = layer->normal[normal].tip;
      /* note that tip has been set to root on terminated normals */
      if (root != tip) gridAddEdge(grid,root,tip,edgeId,DBL_MAX,DBL_MAX);
    }

  }

  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangleNormals(layer, triangle, normals);
    layerTriangle(layer, triangle, nodes);

    /* note that tip has been set to root on terminated normals */
    /* the if (n[0]!=n[3]) checks are for layer termiantion */

    // advance faces
    for (i=0;i<3;i++){
      faceId = layerConstrainedSide(layer, triangle, i);
      if (faceId > 0) {
	side[1] = i;
	side[0] = i+1; if (side[0]>2) side[0] = 0;
	sidenormal[0] = normals[side[0]];
	sidenormal[1] = normals[side[1]];
	sidenode[0]   = nodes[side[0]];
	sidenode[1]   = nodes[side[1]];
	n[0] = layer->normal[sidenormal[0]].root;
	n[1] = layer->normal[sidenormal[1]].root;
	n[2] = layer->normal[sidenormal[0]].tip;
	n[3] = layer->normal[sidenormal[1]].tip;
	if (layerTetrahedraOnly(layer) || n[0]==n[2] || n[1]==n[3]){
	  if (sidenode[0]<sidenode[1]){
	    if (n[1]!=n[3]) 
	      layer->faceInLayer[gridAddFace(grid,n[0],n[1],n[3],faceId)]=TRUE;
	    if (n[0]!=n[2])  
	      layer->faceInLayer[gridAddFace(grid,n[0],n[3],n[2],faceId)]=TRUE;
	  }else{
	    if (n[0]!=n[2])  
	      layer->faceInLayer[gridAddFace(grid,n[0],n[1],n[2],faceId)]=TRUE;
	    if (n[1]!=n[3])  
	      layer->faceInLayer[gridAddFace(grid,n[2],n[1],n[3],faceId)]=TRUE;
	  }
	}else{
	  gridAddQuad(grid,n[0],n[1],n[3],n[2],faceId);
	}
      }
    }

    // advance cells
    /* pg. 82-85 of Garimella Thesis*/
    /* sort so that normals[0] is the smallest normal id*/
    if (nodes[1]<nodes[0] && nodes[1]<nodes[2]){
      node = nodes[1];
      nodes[1] = nodes[2];
      nodes[2] = nodes[0];
      nodes[0] = node;
      normal = normals[1];
      normals[1] = normals[2];
      normals[2] = normals[0];
      normals[0] = normal;
    }
    if (nodes[2]<nodes[0] && nodes[2]<nodes[1]){
      node = nodes[2];
      nodes[2] = nodes[1];
      nodes[1] = nodes[0];
      nodes[0] = node;
      normal = normals[2];
      normals[2] = normals[1];
      normals[1] = normals[0];
      normals[0] = normal;
    }

    nterminated = 0;
    for (i=0;i<3;i++){
      n[i]   = layer->normal[normals[i]].root;
      n[i+3] = layer->normal[normals[i]].tip;
      if (layerNormalTerminated(layer,normals[i])) nterminated++;
    }

    
    if (layerTetrahedraOnly(layer) || nterminated >= 2){
      if (nodes[2]<nodes[1]){
	if (n[0]!=n[3]) 
	  layer->cellInLayer[gridAddCell(grid, n[0], n[4], n[5], n[3])]=TRUE;
	if (n[2]!=n[5]) 
	  layer->cellInLayer[gridAddCell(grid, n[2], n[0], n[4], n[5])]=TRUE;
	if (n[1]!=n[4]) 
	  layer->cellInLayer[gridAddCell(grid, n[2], n[0], n[1], n[4])]=TRUE;
      }else{
	if (n[0]!=n[3]) 
	  layer->cellInLayer[gridAddCell(grid, n[0], n[4], n[5], n[3])]=TRUE;
	if (n[1]!=n[4]) 
	  layer->cellInLayer[gridAddCell(grid, n[0], n[1], n[5], n[4])]=TRUE;
	if (n[2]!=n[5]) 
	  layer->cellInLayer[gridAddCell(grid, n[2], n[0], n[1], n[5])]=TRUE;
      }
    }else{
      if ( nterminated == 1 ) {
	if (n[0]==n[3]) gridAddPyramid(grid,n[1],n[2],n[0],n[4],n[5]);
	if (n[1]==n[4]) gridAddPyramid(grid,n[2],n[0],n[1],n[5],n[3]);
	if (n[2]==n[5]) gridAddPyramid(grid,n[0],n[1],n[2],n[3],n[4]);
      }else{
	gridAddPrism(grid,n[0],n[1],n[2],n[3],n[4],n[5]);
      }
    }
  }


  if (layerNBlend(layer) > 0){
    int blend, blendnormals[4];
    int triangle0, triangle1;
    for (blend=0;blend<layerNBlend(layer);blend++){
      layerBlendNormals(layer, blend, blendnormals );
      
      triangle0 = layerForceTriangle(layer,blendnormals[0],
				     blendnormals[1],blendnormals[2]);
      triangle1 = layerForceTriangle(layer,blendnormals[1],
				     blendnormals[3],blendnormals[2]);

      faceId = layerConstrained(layer,blendnormals[0]);
      if (faceId>0){
	n[0] = layerNormalRoot(layer,blendnormals[0]);
	n[1] = layer->normal[blendnormals[0]].tip;
	n[2] = layer->normal[blendnormals[1]].tip;
	layer->faceInLayer[gridAddFace(grid,n[0],n[1],n[2],faceId)]=TRUE;
	layer->triangle[triangle0].constrainedSide[0]=faceId;
	layer->triangle[triangle0].parentGeomEdge[0] =
	  layer->blend[blend].edgeId[0];
      }
      faceId = layerConstrained(layer,blendnormals[2]);
      if (faceId>0){
	n[0] = layerNormalRoot(layer,blendnormals[2]);
	n[1] = layer->normal[blendnormals[3]].tip;
	n[2] = layer->normal[blendnormals[2]].tip;
	layer->faceInLayer[gridAddFace(grid,n[0],n[1],n[2],faceId)]=TRUE;
	layer->triangle[triangle1].constrainedSide[1]=faceId;
	layer->triangle[triangle1].parentGeomEdge[1] =
	  layer->blend[blend].edgeId[1];
      }

      if ( layerNormalRoot(layer,blendnormals[0]) > 
	   layerNormalRoot(layer,blendnormals[2]) ) {
	normal = blendnormals[0];
	blendnormals[0] = blendnormals[3];
	blendnormals[3] = normal;
	normal = blendnormals[1];
	blendnormals[1] = blendnormals[2];
	blendnormals[2] = normal;
      }

      n[0] = layerNormalRoot(layer,blendnormals[0]);
      n[1] = layer->normal[blendnormals[0]].tip;
      n[2] = layer->normal[blendnormals[1]].tip;
      n[3] = layerNormalRoot(layer,blendnormals[2]);
      n[4] = layer->normal[blendnormals[2]].tip;
      n[5] = layer->normal[blendnormals[3]].tip;

      layer->cellInLayer[gridAddCell(grid, n[0], n[4], n[5], n[3])]=TRUE;
      layer->cellInLayer[gridAddCell(grid, n[2], n[0], n[4], n[5])]=TRUE;
      layer->cellInLayer[gridAddCell(grid, n[2], n[0], n[1], n[4])]=TRUE;

    }
    layerBuildNormalTriangleAdjacency(layer);
    layer->nblend=0;
  }


  for (normal=0;normal<layerNNormal(layer);normal++){
    faceId = layerConstrained(layer,normal);
    if (0 > faceId) {
      edgeId = -faceId;
      gridProjectNodeToEdge(grid, layer->normal[normal].root, edgeId );
      // need to do faces to get the UV vals but have xyz.
    }
    layer->normal[normal].root = layer->normal[normal].tip;
    layer->normal[normal].tip = EMPTY;
    if (0 < faceId) {
      gridProjectNodeToFace(grid, layer->normal[normal].root, faceId );
    }
    if (0 > faceId) {
      edgeId = -faceId;
      gridProjectNodeToEdge(grid, layer->normal[normal].root, edgeId );
      // need to do faces to get the UV vals but have xyz.
    }
    gridFreezeNode(grid,layer->normal[normal].root);
  }

  return layer;
}

Layer *layerWiggle(Layer *layer, double height )
{
  Grid *grid = layer->grid;
  int normal, root, faceId, edgeId, i;
  double xyz[3];

  if (layerNNormal(layer) == 0 ) return NULL;

  for (normal=0;normal<layerNNormal(layer);normal++){
    if ( !layerNormalTerminated(layer,normal) ) {
      root = layer->normal[normal].root;
      gridNodeXYZ(grid,root,xyz);
      for(i=0;i<3;i++)xyz[i]=xyz[i]+height*layer->normal[normal].direction[i];
      gridSetNodeXYZ(grid, root, xyz);
      faceId = layerConstrained(layer,normal);
      if (0 < faceId) {
	gridProjectNodeToFace(grid, layer->normal[normal].root, faceId );
      }
      if (0 > faceId) {
	edgeId = -faceId;
	gridProjectNodeToEdge(grid, layer->normal[normal].root, edgeId );
      }
    }
  }

  return layer;
}

Layer *layerSmoothLayerNeighbors(Layer *layer )
{
  int normal;

  for (normal=0;normal<layerNNormal(layer);normal++)
    gridSmoothNearNode(layer->grid,layer->normal[normal].root);

  return layer;
}

Layer *layerTerminateNormalWithSpacing(Layer *layer, double spacing)
{
  int normal, nterm;;
  if (layerNNormal(layer) == 0 ) return NULL;

  nterm = 0;
  for (normal=0;normal<layerNNormal(layer);normal++){
    if (gridSpacing(layer->grid, layer->normal[normal].root ) < spacing ) {
      nterm++;
      layerTerminateNormal(layer, normal);
    }
  }
  printf("normals %d of %d terminated. %d active.\n",
	 nterm,layerNNormal(layer),layerNActiveNormal(layer) );
  return layer;
}

Layer *layerTerminateNormalWithX(Layer *layer, int direction, double x)
{
  int normal, nterm;
  double xyz[3];
  
  if (layerNNormal(layer) == 0 ) return NULL;

  nterm = 0;
  for (normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(layer->grid, layer->normal[normal].root, xyz);
    if (direction > 0 ) {
      if (xyz[0]>x) { layerTerminateNormal(layer, normal); nterm++; }
    }else{
      if (xyz[0]<x) { layerTerminateNormal(layer, normal); nterm++; }     
    }
  }
  printf("normals %d of %d terminated\n",nterm,layerNNormal(layer) );
  return layer;
}

Layer *layerInsertPhantomTriangle(Layer *layer, double dz )
{
  Grid *grid = layer->grid;
  int normal, faceId, edgeId, newnode;
  double xyz[3], spacing, m;

  layerVisibleNormals(layer,-1.0,-1.0);

  for (normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(grid,layer->normal[normal].root,xyz);
    spacing = gridSpacing(grid, layer->normal[normal].root );
    layer->normal[normal].root = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);
    m = 1.0/spacing/spacing;
    gridSetMap(grid,layer->normal[normal].root,m,0.0,0.0,m,0.0,m); 
    layer->normal[normal].tip = EMPTY;
  }

  layerWiggle(layer, dz );

  for (normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(grid,layer->normal[normal].root,xyz);
    spacing = gridSpacing(grid, layer->normal[normal].root );
    if (0 == layerConstrained(layer,normal)){
      newnode = gridInsertInToVolume(grid, xyz[0], xyz[1], xyz[2]);
    }else{
      faceId = layerConstrained(layer,normal);
      if (0 < faceId) {
	gridForceNodeToFace(grid, layer->normal[normal].root, faceId );
      }
      if (0 > faceId) {
	edgeId = -faceId;
	gridForceNodeToEdge(grid, layer->normal[normal].root, edgeId );
      }
      newnode = gridInsertInToGeomFace(grid, xyz[0], xyz[1], xyz[2]);
    }
    if (EMPTY == newnode) 
      printf("Could not insert node norm %8d x %10.5f y %10.5f z %10.5f\n",
	     normal, xyz[0], xyz[1], xyz[2]);
    layer->normal[normal].tip = newnode;
    gridFreezeNode( grid, newnode );
    m = 1.0/spacing/spacing;
    gridSetMap(grid,newnode,m,0.0,0.0,m,0.0,m); 
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    gridRemoveNode(grid,layer->normal[normal].root);
    layer->normal[normal].root = EMPTY;
  }

  return layer;
}

Layer *layerVerifyPhantomEdges(Layer *layer)
{
  Grid *grid = layer->grid;
  int triangle, normals[3], nline, ngot, n0, n1;

  nline = 0;
  ngot  = 0;
  for(triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangleNormals(layer, triangle, normals );
    n0 = layer->normal[normals[0]].tip;
    n1 = layer->normal[normals[1]].tip;
    if ( n0 != EMPTY && n1 != EMPTY ) {
      nline++;
      if (grid == gridVerifyEdgeExists(grid, n0, n1) ){
	ngot++;
      }else{
	printf("line %d<->%d is missing.\n",n0,n1);
      }
    }
    n0 = layer->normal[normals[1]].tip;
    n1 = layer->normal[normals[2]].tip;
    if ( n0 != EMPTY && n1 != EMPTY ) {
      nline++;
      if (grid == gridVerifyEdgeExists(grid, n0, n1) ){
	ngot++;
      }else{
	printf("line %d<->%d is missing.\n",n0,n1);
      }
    }
    n0 = layer->normal[normals[2]].tip;
    n1 = layer->normal[normals[0]].tip;
    if ( n0 != EMPTY && n1 != EMPTY ) {
      nline++;
      if (grid == gridVerifyEdgeExists(grid, n0, n1) ){
	ngot++;
      }else{
	printf("line %d<->%d is missing.\n",n0,n1);
      }
    }
  }
  printf("lines %d of %d verified.\n",ngot,nline);
 
  return layer;
}

Layer *layerVerifyPhantomFaces(Layer *layer)
{
  Grid *grid = layer->grid;
  int triangle, normals[3], nface, ngot, n0, n1, n2;

  nface = 0;
  ngot  = 0;
  for(triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangleNormals(layer, triangle, normals );
    n0 = layer->normal[normals[0]].tip;
    n1 = layer->normal[normals[1]].tip;
    n2 = layer->normal[normals[2]].tip;
    if ( n0 != EMPTY && n1 != EMPTY && n2 != EMPTY ) {
      nface++;
      if (grid == gridVerifyFaceExists(grid, n0, n1, n2) ){
	ngot++;
      }else{
	printf("face %d %d %d is missing.\n",n0,n1,n2);
      }
    }
  }
  printf("faces %d of %d verified.\n",ngot,nface);
 
  return layer;
}

Layer *layerThaw(Layer *layer)
{
  int normal;
  for(normal=0;normal<layerNNormal(layer);normal++)
    gridThawNode(layerGrid(layer),layerNormalRoot(layer,normal));

  return layer;
}


bool layerTetrahedraOnly(Layer *layer)
{
  return !layer->mixedElementMode;
}

Layer *layerToggleMixedElementMode(Layer *layer)
{
  layer->mixedElementMode = !layer->mixedElementMode;
  return layer;
}

Adj *layerBuildNormalBlendAdjacency(Layer *layer)
{
  Adj *adj;
  int blend, i, normals[4];

  adj = adjCreate( layerNNormal(layer), layerNBlend(layer)*4  );

  for (blend=0;blend<layerNBlend(layer);blend++){
    layerBlendNormals(layer,blend,normals);
    for (i=0;i<4;i++) adjRegister( adj, normals[i], blend );
  }

  return adj;
}

Layer *layerSplitBlend(Layer *layer)
{
  int blend, origblend;
  int normals[4];
  int i, nnode, *newnormal;
  int node0, node1, normal;
  double xyz0[3], xyz1[3], edge[3], n0[3], n1[3], c0[3], c1[3], new[3];
  double length;
  Grid *grid;

  int triangle;
  double direction[3],*norm;

  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    for(i=0;i<3;i++){
      normal = layer->triangle[triangle].normal[i];
      layer->normal[normal].direction[0]=0.0;
      layer->normal[normal].direction[1]=0.0;
      layer->normal[normal].direction[2]=0.0;
    }
  }

  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    for(i=0;i<3;i++){
      normal = layer->triangle[triangle].normal[i];
      if (layer != layerTriangleDirection(layer,triangle,direction))
	printf("Error: layerSplitBlend: %s: %d: %s\n",
	       __FILE__,__LINE__,"triangle direction");
      layer->normal[normal].direction[0] += direction[0];
      layer->normal[normal].direction[1] += direction[1];
      layer->normal[normal].direction[2] += direction[2];
    }
  }

  for (normal=0;normal<layerNNormal(layer);normal++){
    norm = layer->normal[normal].direction;
    length = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    if (length > 0.0) {
      for ( i=0;i<3;i++) norm[i] = norm[i]/length;
    }else{
      printf("zero length normal\n");
      for ( i=0;i<3;i++) norm[i] = 0.0;
    }
  }

  layerVisibleNormals(layer, 0.99, 1.0e-15 );

  origblend = layerNBlend(layer);

  grid = layerGrid(layer);
  nnode = gridNNode(grid);
  newnormal = malloc(nnode*sizeof(int));
  for (i=0; i<nnode; i++) newnormal[i]=EMPTY;

  for (blend=0; blend<layerNBlend(layer); blend++){

    layerBlendNormals(layer,blend,normals);

    node0 = layerNormalRoot(layer,normals[0]);
    node1 = layerNormalRoot(layer,normals[2]);
    if (newnormal[node0]==EMPTY){
      newnormal[node0] = layerDuplicateNormal(layer, normals[0] );
      gridNodeXYZ(grid,node0,xyz0);
      gridNodeXYZ(grid,node1,xyz1);
      gridSubtractVector(xyz1,xyz0,edge);
      layerNormalDirection(layer,normals[0],n0);
      layerNormalDirection(layer,normals[1],n1);
      gridCrossProduct(edge,n0,c0);
      gridCrossProduct(n1,edge,c1);
      for (i=0; i<3; i++) new[i] = 0.5*(c0[i]+c1[i]);
      length = sqrt(gridDotProduct(new,new));
      for (i=0; i<3; i++) 
	layer->normal[newnormal[node0]].direction[i]=new[i]/length;
    }

    normal = normals[0];
    normals[0] = normals[3];
    normals[3] = normal;
    normal = normals[1];
    normals[1] = normals[2];
    normals[2] = normal;

    node0 = layerNormalRoot(layer,normals[0]);
    node1 = layerNormalRoot(layer,normals[2]);
    if (newnormal[node0]==EMPTY){
      newnormal[node0] = layerDuplicateNormal(layer, normals[0] );
      gridNodeXYZ(grid,node0,xyz0);
      gridNodeXYZ(grid,node1,xyz1);
      gridSubtractVector(xyz1,xyz0,edge);
      layerNormalDirection(layer,normals[0],n0);
      layerNormalDirection(layer,normals[1],n1);
      gridCrossProduct(edge,n0,c0);
      gridCrossProduct(n1,edge,c1);
      for (i=0; i<3; i++) new[i] = 0.5*(c0[i]+c1[i]);
      length = sqrt(gridDotProduct(new,new));
      for (i=0; i<3; i++) 
	layer->normal[newnormal[node0]].direction[i]=new[i]/length;
    }
  }

  layerDuplicateAllBlend(layer);
  for (blend=0; blend<origblend; blend++){
    layerBlendNormals(layer,blend,normals);
    node0 = layerNormalRoot(layer,normals[0]);
    node1 = layerNormalRoot(layer,normals[2]);
    if (newnormal[node0]!=EMPTY) layer->blend[blend].normal[0]=newnormal[node0];
    if (newnormal[node1]!=EMPTY) layer->blend[blend].normal[2]=newnormal[node1];
  }  
  for (blend=origblend; blend<layerNBlend(layer); blend++){
    layerBlendNormals(layer,blend,normals);
    node0 = layerNormalRoot(layer,normals[1]);
    node1 = layerNormalRoot(layer,normals[3]);
    if (newnormal[node0]!=EMPTY) layer->blend[blend].normal[1]=newnormal[node0];
    if (newnormal[node1]!=EMPTY) layer->blend[blend].normal[3]=newnormal[node1];
  }

  free(newnormal);
  return layer;
}

Layer *layerBlend(Layer *layer)
{
  int normal, originalNormals;
  AdjIterator it;
  int triangle, splitTriangle, nextTriangle, previousTriangle;
  double edgeAngle, largestEdgeAngle, angleLimit;
  int commonEdge[2];

  int newNormal, i;
  bool done;

  angleLimit = 250; /* deg */

  originalNormals = layerNNormal(layer);
  for ( normal = 0 ; normal < originalNormals ; normal++ ) {
    
    largestEdgeAngle = 0;
    splitTriangle = EMPTY;
    for ( it = adjFirst(layer->adj,normal); 
	  adjValid(it); 
	  it = adjNext(it) ){
      triangle = adjItem(it);
      nextTriangle = layerNextTriangle(layer, normal, triangle);
      edgeAngle = layerEdgeAngle(layer,triangle,nextTriangle);
      if (largestEdgeAngle <= edgeAngle){
	largestEdgeAngle = edgeAngle;
	splitTriangle = triangle;
      }
      previousTriangle = layerPreviousTriangle(layer, normal, triangle);
      edgeAngle = layerEdgeAngle(layer,previousTriangle,triangle);
      if (largestEdgeAngle <= edgeAngle){
	largestEdgeAngle = edgeAngle;
	splitTriangle = previousTriangle;
      }
    }
   
    if (splitTriangle != EMPTY && largestEdgeAngle > angleLimit){
      newNormal = layerDuplicateNormal(layer, normal );
      nextTriangle = layerNextTriangle(layer, normal, splitTriangle);
      layerCommonEdge(layer, splitTriangle, nextTriangle, commonEdge);
      if (layerNormalRoot(layer,normal) == commonEdge[0] ) {
	layerAddBlend(layer,newNormal,normal,commonEdge[1]);
      }else{
	layerAddBlend(layer,newNormal,normal,commonEdge[0]);
      }
      triangle = nextTriangle;
      done = FALSE;
      while (!done) {
	for (i=0;i<3;i++)
	  if (layer->triangle[triangle].normal[i] == normal) 
	    layer->triangle[triangle].normal[i] = newNormal;
	nextTriangle = layerNextTriangle(layer, normal, triangle);
	done = (EMPTY == nextTriangle);
	if (!done && angleLimit < layerEdgeAngle(layer,triangle,nextTriangle)){
	  layerCommonEdge(layer, triangle, nextTriangle, commonEdge);
	  if (layerNormalRoot(layer,normal) == commonEdge[0] ) {
	    layerAddBlend(layer,normal,newNormal,commonEdge[1]);
	  }else{
	    layerAddBlend(layer,normal,newNormal,commonEdge[0]);
	  }
	  done = TRUE;
	}
	triangle = nextTriangle;
      }
    }
 
  }

  layerBuildNormalTriangleAdjacency(layer);
  layerVisibleNormals(layer , 
		      sin(ConvertDegreeToRadian(largestEdgeAngle*0.333)), 
		      1.0e-8 );
  return layer;
}

Layer *layerDuplicateAllBlend(Layer *layer)
{
  int blend, maxblend;
  
  maxblend = layerNBlend(layer);
  for (blend=0; blend<maxblend; blend++){
    
    if (layer->nblend >= layer->maxblend) {
      layer->maxblend += 5000;
      if (layer->blend == NULL) {
	layer->blend = malloc(layer->maxblend*sizeof(Blend));
      }else{
	layer->blend = realloc(layer->blend,layer->maxblend*sizeof(Blend));
      }
    }

    layer->blend[layer->nblend].nodes[0] = 
      layer->blend[blend].nodes[0];
    layer->blend[layer->nblend].nodes[1] = 
      layer->blend[blend].nodes[1];

    layer->blend[layer->nblend].normal[0] = 
      layer->blend[blend].normal[0];
    layer->blend[layer->nblend].normal[1] = 
      layer->blend[blend].normal[1];
    layer->blend[layer->nblend].normal[2] = 
      layer->blend[blend].normal[2];
    layer->blend[layer->nblend].normal[3] = 
      layer->blend[blend].normal[3];
    
    layer->blend[layer->nblend].edgeId[0] = 
      layer->blend[blend].edgeId[0];
    layer->blend[layer->nblend].edgeId[1] = 
      layer->blend[blend].edgeId[1];

    layer->nblend++;
  }

  return layer;
}

Layer *layerAddBlend(Layer *layer, int normal0, int normal1, int otherNode )
{
  int i, node0, node1, n0, n1;
  bool newEdge;
  int blend;
  int triangle, edge, nodes[2], excludeId, edgeId, side;
  AdjIterator it;
  Grid *grid;
  grid=layerGrid(layer);

  if (layer->nblend >= layer->maxblend) {
    layer->maxblend += 5000;
    if (layer->blend == NULL) {
      layer->blend = malloc(layer->maxblend*sizeof(Blend));
    }else{
      layer->blend = realloc(layer->blend,layer->maxblend*sizeof(Blend));
    }
  }
  
  n0 = layerNormalRoot(layer,normal0);
  n1 = otherNode;

  newEdge=TRUE;
  for (i=0;i<layerNBlend(layer);i++){
    node0 = layer->blend[i].nodes[0];
    node1 = layer->blend[i].nodes[1];

    if ( ( node0 == n0 && node1 == n1 ) ||
	 ( node0 == n1 && node1 == n0 ) ){
      newEdge=FALSE;
      blend = i;
      layer->blend[blend].normal[2] = normal1;
      layer->blend[blend].normal[3] = normal0;
    }
  }

  if (newEdge){
    blend = layer->nblend;
    layer->blend[blend].nodes[0] = n0;
    layer->blend[blend].nodes[1] = n1;
    layer->blend[blend].normal[0] = normal0;
    layer->blend[blend].normal[1] = normal1;
    layer->blend[blend].edgeId[0] = EMPTY;
    layer->blend[blend].edgeId[1] = EMPTY;
    layer->nblend++;
  }

  edge = gridFindEdge(grid,n0,n1);
  excludeId=0;
  gridEdge(grid,edge,nodes,&excludeId);
  
  side = newEdge?0:1;

  for ( it = adjFirst(layer->adj,normal0); adjValid(it); it=adjNext(it) ){
    triangle = adjItem(it);
    for (i=0;i<3;i++) {
      edgeId = layer->triangle[triangle].parentGeomEdge[i];
      if (edgeId > 0 && edgeId != excludeId) 
	layer->blend[blend].edgeId[side] = edgeId;
    }
  }
  for ( it = adjFirst(layer->adj,normal1); adjValid(it); it=adjNext(it) ){
    triangle = adjItem(it);
    for (i=0;i<3;i++) {
      edgeId = layer->triangle[triangle].parentGeomEdge[i];
      if (edgeId > 0 && edgeId != excludeId) 
	layer->blend[blend].edgeId[side] = edgeId;
    }
  }

  return layer;
}

Layer *layerBlendNormals(Layer *layer, int blend, int *normals )
{
  int i;
  if (blend < 0 || blend >= layerNBlend(layer)) return NULL;
  for(i=0;i<4;i++) normals[i] = layer->blend[blend].normal[i];
  return layer;
}

Layer *layerBlendExtend(Layer *layer, double dx, double dy, double dz )
{

  int blend, i, fix, fixblend, node, newnode, normal, newnormal;
  int normals[4];
  double xyz[3];
  Adj *adj;
  AdjIterator it;
  Grid *grid;

  grid = layerGrid(layer);
  adj = layerBuildNormalBlendAdjacency(layer);

  for (blend=0; blend < layerNBlend(layer); blend++){
    for(i=0;i<4;i++){
      layer->blend[blend].oldnormal[i]=layer->blend[blend].normal[i];
      layer->blend[blend].normal[i]=EMPTY;
    }
  }

  for (blend=0; blend < layerNBlend(layer); blend++){
    for(i=0;i<4;i++){
      if (layer->blend[blend].normal[i]==EMPTY){
	normal = layer->blend[blend].oldnormal[i];
	node = layerNormalRoot(layer,normal);
	gridNodeXYZ(grid,node,xyz);
	xyz[0] += dx; xyz[1] += dy; xyz[2] += dz;
	newnode = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);
	newnormal = layerDuplicateNormal(layer, normal);
	layer->normal[newnormal].root = newnode;
	for ( it = adjFirst(adj,normal); adjValid(it); it=adjNext(it) ){
	  fixblend = adjItem(it);
	  for (fix=0;fix<4;fix++)
	    if ( layer->blend[fixblend].oldnormal[fix] == normal )
	      layer->blend[fixblend].normal[fix] = newnormal;
	}
      }
    }
  }

  for (blend=0; blend < layerNBlend(layer); blend++){
    layer->blend[blend].nodes[0] = 
      layerNormalRoot(layer,layer->blend[blend].normal[0]);
    layer->blend[blend].nodes[1] = 
      layerNormalRoot(layer,layer->blend[blend].normal[2]);
    normals[0] = layer->blend[blend].oldnormal[0];
    normals[1] = layer->blend[blend].oldnormal[2];
    normals[2] = layer->blend[blend].normal[0];
    normals[3] = layer->blend[blend].normal[2];
    printf("norm %d %d %d %d\n",normals[0],normals[1],normals[2],normals[3]);
    layerForceTriangle(layer,normals[0],normals[1],normals[2]);
    layerForceTriangle(layer,normals[1],normals[3],normals[2]);
    normals[0] = layer->blend[blend].oldnormal[3];
    normals[1] = layer->blend[blend].oldnormal[1];
    normals[2] = layer->blend[blend].normal[3];
    normals[3] = layer->blend[blend].normal[1];
    printf("norm %d %d %d %d\n",normals[0],normals[1],normals[2],normals[3]);
    layerForceTriangle(layer,normals[0],normals[1],normals[2]);
    layerForceTriangle(layer,normals[1],normals[3],normals[2]);
  }

  layerBuildNormalTriangleAdjacency(layer);

  return layer;
}
