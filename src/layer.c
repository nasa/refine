
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

struct Layer {
  Grid *grid;
  int maxtriangle, ntriangle;
  Triangle *triangle;
  int maxParentGeomFace, nParentGeomFace, *ParentGeomFace;
  int maxnormal, nnormal;
  Normal *normal;
  int *globalNode2Normal;
  int nConstrainingGeometry, *constrainingGeometry;
  Adj *adj;
  bool mixedElementMode;
};

Layer *layerCreate( Grid *grid )
{
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
  layer->maxnormal=0;
  layer->nnormal=0;
  layer->normal=NULL;
  layer->globalNode2Normal=NULL;
  layer->nConstrainingGeometry=0;
  layer->constrainingGeometry=NULL;
  layer->adj=NULL;
  layer->mixedElementMode=FALSE;
  return layer;
}

Layer *formAdvancingFront( Grid *grid, char *project )
{
  Layer *layer;
  int i, nbc, bc[4];
  bool box, plate, om6, n12, quarter, spherecone;
  
  box = (NULL != strstr( project, "box"));
  plate = (NULL != strstr( project, "plate"));
  om6 = (NULL != strstr( project, "om6"));
  n12 = (NULL != strstr( project, "n12"));
  quarter = (NULL != strstr( project, "quarter"));
  spherecone = (NULL != strstr( project, "Mach6sphere"));

  if (box) printf("string %s has box.\n",project);
  if (plate) printf("string %s has plate.\n",project);
  if (om6) printf("string %s has om6.\n",project);
  if (n12) printf("string %s has n12.\n",project);
  if (quarter) printf("string %s has quarter.\n",project);
  if (spherecone) printf("string %s has spherecone (Mach6sphere).\n",project);

  bc[0]=1;
  bc[1]=2;
  bc[2]=3;
  if(box) nbc = 1;
  if(om6) nbc = 2;
  if(plate) nbc = 3;
  if(n12){
    nbc=2;
    bc[0]=5;
    bc[1]=6;
  }
  if(quarter) nbc = 1;
  if(spherecone){
    nbc=4;
    bc[0]=1;
    bc[1]=3;
    bc[2]=4;
    bc[3]=5;
  }
  /* skip freezing
  printf("freezing distant volume nodes.\n");
  gridFreezeAll(grid);
  for (i=0;i<nbc;i++){  
    printf("thaw bc %d.\n",bc[i]);
    gridThawNearBC(grid,0.5,bc[i]);
    gridFreezeBCFace(grid,bc[i]);
  }
  */
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
  if (plate) {
    layerTerminateNormalWithX(layer,-1, 0.0001);
    layerTerminateNormalWithX(layer, 1, 1.0001);

    /*
      face and four edges
      04 13 12 07 11
      05 15 11 02 14
      06 17 14 10 16

      14 20 26 16 08

      07 20 09 19 18
      08 19 04 22 21
      09 24 23 22 05

      13 12 28 24 06
    */
    layerConstrainNormal(layer,-11);
    layerConstrainNormal(layer,-14);
    layerConstrainNormal(layer,-16);

    layerConstrainNormal(layer,-20);

    layerConstrainNormal(layer,-19);
    layerConstrainNormal(layer,-22);
    layerConstrainNormal(layer,-24);

    layerConstrainNormal(layer,-12);
    /* y1 */
    layerConstrainNormal(layer,4);
    layerConstrainNormal(layer,5);
    layerConstrainNormal(layer,6);
    /* y0 */
    layerConstrainNormal(layer,7);
    layerConstrainNormal(layer,8);
    layerConstrainNormal(layer,9);
    /* inout */
    layerConstrainNormal(layer,13);
    layerConstrainNormal(layer,14);
  }
  if (om6) {
    layerConstrainNormal(layer,5);
  }
  if(n12){
    layerConstrainNormal(layer,1);
    layerConstrainNormal(layer,2);
  }
  if(quarter){
    layerConstrainNormal(layer,2);
    layerConstrainNormal(layer,4);
    layerConstrainNormal(layer,-5);
    layerConstrainNormal(layer,-6);
  }
  if(spherecone){
    layerConstrainNormal(layer,2);
  }
  printf("make advancing layer triangle normals visible to triangle.\n");
  layerVisibleNormals(layer,-1.0,-1.0);
  return layer;
}

Grid *layerGrid(Layer *layer)
{
  return layer->grid;
}

void layerFree(Layer *layer)
{
  gridDetachNodeSorter( layer->grid );
  if (layer->adj != NULL) adjFree(layer->adj);
  if (layer->constrainingGeometry != NULL) free(layer->constrainingGeometry);
  if (layer->globalNode2Normal != NULL) free(layer->globalNode2Normal);
  if (layer->normal != NULL) free(layer->normal);
  if (layer->ParentGeomFace != NULL) free(layer->ParentGeomFace);
  if (layer->triangle != NULL) free(layer->triangle);
  free(layer);
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
      }
    }
  }

  layer->adj = adjCreate( layer->nnormal,layerNTriangle(layer)*3  );
  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    for(i=0;i<3;i++){
      normal = layer->triangle[triangle].normal[i];
      adjRegister( layer->adj, normal, triangle );
      layerTriangleDirection(layer,triangle,direction);
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

  layer->ntriangle++;

  return layer;
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
  if (normal < 0 || normal >= layerNNormal(layer) ) return 0;
  return layer->normal[normal].root;
}

int layerNormalDeg(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return 0;
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

Layer *layerNormalDirection(Layer *layer, int normal, double *direction )
{
  int i;
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  for ( i=0;i<3;i++) direction[i] = layer->normal[normal].direction[i];

  return layer;
}

Layer *layerSetHeightOfAllNormals(Layer *layer, double height )
{
  int normal;

  for(normal=0;normal<layerNNormal(layer);normal++)
    layerSetNormalHeight( layer, normal, height );

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
    if (mindot <= 0.0 ) 
      printf("ERROR: %s, %d, Invisible normal %d, min dot product %f\n",
	     __FILE__, __LINE__, normal, mindot);
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
  int triangle, normals[3], n[6], side[2];
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
      gridReconnectCellUnlessFrozen(grid, root, tip);
      faceId = layerConstrained(layer,normal);
      if (0 > faceId) {
	edgeId = -faceId;
	gridReconnectEdgeUnlessFrozen(grid, edgeId, root, tip);
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
	gridReconnectFaceUnlessFrozen(grid, faceId, 
				      layer->normal[normal0].root, 
				      layer->normal[normal0].tip);
	gridReconnectFaceUnlessFrozen(grid, faceId, 
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

    /* note that tip has been set to root on terminated normals */
    /* the if (n[0]!=n[3]) checks are for layer termiantion */

    // advance faces
    for (i=0;i<3;i++){
      faceId = layerConstrainedSide(layer, triangle, i);
      if (faceId > 0) {
	side[1] = i;
	side[0] = i+1; if (side[0]>2) side[0] = 0;
	side[0] = normals[side[0]];
	side[1] = normals[side[1]];
	n[0] = layer->normal[side[0]].root;
	n[1] = layer->normal[side[1]].root;
	n[2] = layer->normal[side[0]].tip;
	n[3] = layer->normal[side[1]].tip;
	if (layerTetrahedraOnly(layer) || n[0]==n[2] || n[1]==n[3]){
	  if (side[0]<side[1]){
	    if (n[1]!=n[3]) gridAddFace(grid,n[0],n[1],n[3],faceId);
	    if (n[0]!=n[2]) gridAddFace(grid,n[0],n[3],n[2],faceId);
	  }else{
	    if (n[0]!=n[2]) gridAddFace(grid,n[0],n[1],n[2],faceId);
	    if (n[1]!=n[3]) gridAddFace(grid,n[2],n[1],n[3],faceId);
	  }
	}else{
	  gridAddQuad(grid,n[0],n[1],n[3],n[2],faceId);
	}
      }
    }

    // advance cells
    /* pg. 82-85 of Garimella Thesis*/
    /* sort so that normals[0] is the smallest normal id*/
    if (normals[1]<normals[0] && normals[1]<normals[2]){
      normal = normals[1];
      normals[1] = normals[2];
      normals[2] = normals[0];
      normals[0] = normal;
    }
    if (normals[2]<normals[0] && normals[2]<normals[1]){
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
      if (normals[2]<normals[1]){
	if (n[0]!=n[3]) gridAddCell(grid, n[0], n[4], n[5], n[3]);
	if (n[2]!=n[5]) gridAddCell(grid, n[2], n[0], n[4], n[5]);
	if (n[1]!=n[4]) gridAddCell(grid, n[2], n[0], n[1], n[4]);
      }else{
	if (n[0]!=n[3]) gridAddCell(grid, n[0], n[4], n[5], n[3]);
	if (n[1]!=n[4]) gridAddCell(grid, n[0], n[1], n[5], n[4]);
	if (n[2]!=n[5]) gridAddCell(grid, n[2], n[0], n[1], n[5]);
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
  printf("normals %d of %d terminated\n",nterm,layerNNormal(layer) );
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

bool layerTetrahedraOnly(Layer *layer)
{
  return !layer->mixedElementMode;
}

Layer *layerToggleMixedElementMode(Layer *layer)
{
  layer->mixedElementMode = !layer->mixedElementMode;
  return layer;
}
