/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <values.h>
#include <string.h>
#include "layer.h"
#include "gridmath.h"
#include "gridmetric.h"
#include "gridcad.h"
#include "gridinsert.h"
#include "intersect.h"

Layer *layerCreate( Grid *grid )
{
  int i;
  Layer *layer;
  layer = malloc(sizeof(Layer));
  layer->grid = grid;
  gridAttachPacker( grid, layerPack, (void *)layer );
  gridAttachNodeSorter( grid, layerSortNodes, (void *)layer );
  gridAttachReallocator( grid, layerReallocator, (void *)layer );
  gridAttachFreeNotifier( grid, layerGridHasBeenFreed, (void *)layer );
  layer->maxtriangle=0;
  layer->ntriangle=0;
  layer->triangle=NULL;
  layer->triangleAdj=NULL;
  layer->maxParentGeomFace=0;
  layer->nParentGeomFace=0;
  layer->ParentGeomFace=NULL;
  layer->maxblend=0;
  layer->nblend=0;
  layer->blend=NULL;
  layer->blendAdj=NULL;
  layer->maxnormal=0;
  layer->nnormal=0;
  layer->originalnormal=0;
  layer->normal=NULL;
  layer->globalNode2Normal=NULL;
  layer->vertexNormal=NULL;
  layer->nConstrainingGeometry=0;
  layer->constrainingGeometry=NULL;
  layer->nearTree=NULL;
  layer->mixedElementMode=FALSE;

  layer->normalTriangleHub = EMPTY;
  layer->normalTriangleDegree = EMPTY;

  layer->cellInLayer = malloc(gridMaxCell(grid)*sizeof(GridBool));
  for (i=0;i<gridMaxCell(grid);i++) layer->cellInLayer[i] = FALSE;
  layer->faceInLayer = malloc(gridMaxFace(grid)*sizeof(GridBool));
  for (i=0;i<gridMaxFace(grid);i++) layer->faceInLayer[i] = FALSE;
  layer->edgeInLayer = malloc(gridMaxEdge(grid)*sizeof(GridBool));
  for (i=0;i<gridMaxEdge(grid);i++) layer->edgeInLayer[i] = FALSE;

  layer->tecplotFile = NULL;

  return layer;
}

Grid *layerGrid(Layer *layer)
{
  return layer->grid;
}

void layerFree(Layer *layer)
{
  if ( layer->tecplotFile != NULL ) fclose(layer->tecplotFile);
  free(layer->edgeInLayer);
  free(layer->faceInLayer);
  free(layer->cellInLayer);
  if (NULL != layer->grid) {
    gridDetachPacker( layer->grid );
    gridDetachNodeSorter( layer->grid );
    gridDetachReallocator( layer->grid );
    gridDetachFreeNotifier( layer->grid );
  }
  if (layer->nearTree != NULL) free(layer->nearTree);
  if (layer->constrainingGeometry != NULL) free(layer->constrainingGeometry);
  if (layer->vertexNormal != NULL) free(layer->vertexNormal);
  if (layer->globalNode2Normal != NULL) free(layer->globalNode2Normal);
  if (layer->normal != NULL) free(layer->normal);
  if (layer->ParentGeomFace != NULL) free(layer->ParentGeomFace);
  if (layer->blendAdj != NULL) adjFree(layer->blendAdj);
  if (layer->blend != NULL) free(layer->blend);
  if (layer->triangleAdj != NULL) adjFree(layer->triangleAdj);
  if (layer->triangle != NULL) free(layer->triangle);
  free(layer);
}

Layer *formAdvancingFront( Grid *grid, char *project )
{
  Layer *layer;
  int nbc, bc[2];
  GridBool box, om6, n12;
  
  box = (NULL != strstr( project, "box"));
  om6 = (NULL != strstr( project, "om6"));
  n12 = (NULL != strstr( project, "n12"));

  if (box) printf("string %s has box.\n",project);
  if (om6) printf("string %s has om6.\n",project);
  if (n12) printf("string %s has n12.\n",project);

  bc[0]=1;
  bc[1]=2;
  nbc=0;
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

void layerPack(void *voidLayer, 
	       int nnode, int maxnode, int *nodeo2n,
	       int ncell, int maxcell, int *cello2n,
	       int nface, int maxface, int *faceo2n,
	       int nedge, int maxedge, int *edgeo2n)
{
  Layer *layer = (Layer *)voidLayer;
  int normal;
  int packcell, origcell;
  int packface, origface;
  int packedge, origedge;

  for (normal = 0 ; normal < layerNNormal(layer) ; normal++ ) {
    if (EMPTY != layer->normal[normal].root)
      layer->normal[normal].root = nodeo2n[layer->normal[normal].root];
    if (EMPTY != layer->normal[normal].tip)
      layer->normal[normal].tip = nodeo2n[layer->normal[normal].tip];
  }

  for (origcell=0; origcell<maxcell; origcell++) {
    packcell = cello2n[origcell];
    if (packcell != EMPTY)
      layer->cellInLayer[packcell]=layer->cellInLayer[origcell];
  }
  for ( packcell=ncell ; packcell < maxcell ; packcell++ ){ 
    layer->cellInLayer[packcell] = FALSE;
  }

  for (origface=0; origface<maxface; origface++) {
    packface = faceo2n[origface];
    if (packface != EMPTY)
      layer->faceInLayer[packface]=layer->faceInLayer[origface];
  }
  for ( packface=nface ; packface < maxface ; packface++ ){ 
    layer->faceInLayer[packface] = FALSE;
  }

  for (origedge=0; origedge<maxedge; origedge++) {
    packedge = edgeo2n[origedge];
    if (packedge != EMPTY)
      layer->edgeInLayer[packedge]=layer->edgeInLayer[origedge];
  }
  for ( packedge=nedge ; packedge < maxedge ; packedge++ ){
    layer->edgeInLayer[packedge] = FALSE;
  }
}

void layerSortNodes(void *voidLayer, int maxnode, int *o2n)
{
  Layer *layer = (Layer *)voidLayer;
  int normal;

  for (normal = 0 ; normal < layerNNormal(layer) ; normal++ ) {
    if (EMPTY != layer->normal[normal].root)
      layer->normal[normal].root = o2n[layer->normal[normal].root];
    if (EMPTY != layer->normal[normal].tip)
      layer->normal[normal].tip = o2n[layer->normal[normal].tip];
  }

}

void layerReallocator(void *voidLayer, int reallocType, 
		      int lastSize, int newSize)
{
  Layer *layer = (Layer *)voidLayer;
  int i;

  switch (reallocType) {
  case gridREALLOC_EDGE:
    layer->edgeInLayer = realloc(layer->edgeInLayer, newSize*sizeof(GridBool));
    for (i=lastSize;i<newSize;i++) layer->edgeInLayer[i] = FALSE;
    break;
  case gridREALLOC_FACE:
    layer->faceInLayer = realloc(layer->faceInLayer, newSize*sizeof(GridBool));
    for (i=lastSize;i<newSize;i++) layer->faceInLayer[i] = FALSE;
    break;
  case gridREALLOC_CELL:
    layer->cellInLayer = realloc(layer->cellInLayer, newSize*sizeof(GridBool));
    for (i=lastSize;i<newSize;i++) layer->cellInLayer[i] = FALSE;
    break;
  }

}

void layerGridHasBeenFreed(void *voidLayer )
{
  Layer *layer = (Layer *)voidLayer;
  layer->grid = NULL;
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
  int ibc, face, id, nodes[3];
  Grid *grid;

  grid = layerGrid(layer);

  for(ibc=0;ibc<nbc;ibc++) layerAddParentGeomFace(layer,bc[ibc]);

  for (ibc=0;ibc<nbc;ibc++){
    for(face=0;face<gridMaxFace(grid);face++){
      if (grid == gridFace(grid,face,nodes,&id) &&
	  id==bc[ibc] ) {
	layerAddTriangle(layer,nodes[0],nodes[1],nodes[2]);
	layer->faceInLayer[face]=TRUE;
      }
    }
  }

  layerBuildNormalTriangleAdjacency(layer);

  layerInitializeTriangleNormalDirection(layer);

  return layer;

}

Layer *layerInitializeTriangleNormalDirection(Layer *layer)
{
  int triangle, normal, i;
  double direction[3];

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
	printf("Error: layerInitializeTriangleNormalDirection: %s: %d: %s\n",
	       __FILE__,__LINE__,"NULL layerTriangleDirection");
      layer->normal[normal].direction[0] += direction[0];
      layer->normal[normal].direction[1] += direction[1];
      layer->normal[normal].direction[2] += direction[2];
    }
  }

  for (normal=0;normal<layerNNormal(layer);normal++)
    gridVectorNormalize(layer->normal[normal].direction);
  
  return layer;
}

Layer *layerBuildNormalTriangleAdjacency(Layer *layer)
{
  int triangle, i, normals[3];

  layer->normalTriangleHub = EMPTY;
  if (NULL != layer->triangleAdj) adjFree(layer->triangleAdj);
  layer->triangleAdj = adjCreate( layerNNormal(layer), 
				  layerNTriangle(layer)*3, 1000 );

  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangleNormals(layer,triangle,normals);
    for (i=0;i<3;i++) adjRegister( layer->triangleAdj, normals[i], triangle );
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

GridBool layerParentGeomFace(Layer *layer, int faceId )
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
  double edge1[3], edge2[3], norm[3]; 
  
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

  for ( i=0;i<3;i++) direction[i] = norm[i];
  gridVectorNormalize(direction);

  return layer;
}

Layer *layerTriangleArea(Layer *layer, int triangle, double *area )
{
  int i, nodes[3];
  double node0[3], node1[3], node2[3];
  double edge1[3], edge2[3], norm[3]; 
  
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

  *area = gridVectorLength(norm) * 0.5;

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

Layer *layerTriangleFourthNode(Layer *layer, int triangle, double *xyz )
{
  int i;
  double center[3], h;
  double direction[3];

  if ( layer != layerTriangleCenter(layer,triangle,center) ) return NULL;
  layerTriangleDirection(layer,triangle,direction);
  layerTriangleMaxEdgeLength(layer,triangle,&h);
  
  for (i=0;i<3;i++) xyz[i] = center[i] + 0.5*h*direction[i];

  return layer;

}

Layer *layerTriangleInviscidTet(Layer *layer, int triangle, 
				double *node0, double *node1,
				double *node2, double *node3) 
{
  int nodes[3];

  if ( layer != layerTriangle(layer,triangle,nodes) ) return NULL;

  if (layer->grid != gridNodeXYZ( layer->grid, nodes[0], node0 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[1], node1 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[2], node2 )) return NULL;

  if ( layer != layerTriangleFourthNode(layer,triangle,node3) ) return NULL;

  return layer;
}

Layer *layerTriangleMaxEdgeLength(Layer *layer, int triangle, double *length )
{
  int nodes[3];
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

Layer *layerAdvancedTriangleMaxEdgeLength(Layer *layer, 
					  int triangle, double *length )
{
  int i, nodes[3], normals[3];
  double node0[3], node1[3], node2[3];
  double edge0[3], edge1[3], edge2[3], maxLength; 
  
  if ( layer != layerTriangle(layer,triangle,nodes) ) return NULL;
  if ( layer != layerTriangleNormals(layer,triangle,normals) ) return NULL;

  if (layer->grid != gridNodeXYZ( layer->grid, nodes[0], node0 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[1], node1 )) return NULL;
  if (layer->grid != gridNodeXYZ( layer->grid, nodes[2], node2 )) return NULL;

  for(i=0;i<3;i++){
    if (!layerNormalTerminated(layer,normals[0]))
      node0[i] += ( layer->normal[normals[0]].height *
		    layer->normal[normals[0]].direction[i] );
    if (!layerNormalTerminated(layer,normals[1]))
      node1[i] += ( layer->normal[normals[1]].height *
		    layer->normal[normals[1]].direction[i] );
    if (!layerNormalTerminated(layer,normals[2]))
      node2[i] += ( layer->normal[normals[2]].height *
		    layer->normal[normals[2]].direction[i] );
  }

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

Layer *layerNormalMaxEdgeLength(Layer *layer, int normal, double *length )
{
  int triangle;
  double maxLength, currentLength;
  AdjIterator it;

  maxLength = 0;

  for ( it = adjFirst(layer->triangleAdj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
     triangle = adjItem(it);
     layerTriangleMaxEdgeLength(layer, triangle, &currentLength );
     maxLength = MAX(maxLength, currentLength);
  }

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
  layer->normal[normal].initialheight = layer->normal[normal].height;
  layer->normal[normal].length = 0.0;
  layer->normal[normal].maxlength = 1.0;
  layer->normal[normal].rate = 1.0;
  layer->normal[normal].terminated = FALSE;
  layer->normal[normal].directionFrozen = FALSE;

  return layer;
}

int layerDuplicateNormal(Layer *layer, int normal)
{
  int newone, root;

  root = layerNormalRoot(layer,normal);
  if (EMPTY == root) return EMPTY;
 
  newone = layerAddNormal(layer, root);
  if (EMPTY == newone) return EMPTY;

  layer->normal[newone].constrained =     layer->normal[normal].constrained;
  layer->normal[newone].root =            layer->normal[normal].root;
  layer->normal[newone].tip =             layer->normal[normal].tip;
  layer->normal[newone].direction[0] =    layer->normal[normal].direction[0];
  layer->normal[newone].direction[1] =    layer->normal[normal].direction[1];
  layer->normal[newone].direction[2] =    layer->normal[normal].direction[2];
  layer->normal[newone].height =          layer->normal[normal].height;
  layer->normal[newone].initialheight =   layer->normal[normal].initialheight;
  layer->normal[newone].length =          layer->normal[normal].length;
  layer->normal[newone].maxlength =       layer->normal[normal].maxlength;
  layer->normal[newone].rate =            layer->normal[normal].rate;
  layer->normal[newone].terminated =      layer->normal[normal].terminated;
  layer->normal[newone].directionFrozen = layer->normal[normal].directionFrozen;

  return newone;
}

int layerAddNormal(Layer *layer, int globalNodeId )
{
  int i;
  int oldNNormal;

  if (globalNodeId < 0 || globalNodeId >= layerMaxNode(layer) ) return EMPTY;

  if (layer->nnormal >= layer->maxnormal) {
    oldNNormal = layer->maxnormal;
    layer->maxnormal += 5000;
    if (layer->normal == NULL) {
      layer->normal = malloc(layer->maxnormal*sizeof(Normal));
      layer->globalNode2Normal = malloc(layerMaxNode(layer)*sizeof(int));
      for (i=0;i<layerMaxNode(layer);i++) layer->globalNode2Normal[i]=EMPTY;
    }else{
      layer->normal = realloc(layer->normal,layer->maxnormal*sizeof(Normal));
      if ( NULL != layer->vertexNormal ) {
	layer->vertexNormal = realloc(layer->vertexNormal, 
				      layer->maxnormal*sizeof(int));
	for (i=oldNNormal;i<layer->maxnormal;i++) layer->vertexNormal[i]=EMPTY;
      }
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

int layerNormalTip(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return EMPTY;
  return layer->normal[normal].tip;
}

int layerNormalDeg(Layer *layer, int normal )
{
  if (NULL == layer->triangleAdj) return 0;
  return adjDegree(layer->triangleAdj, normal);
}

double layerNormalAngle(Layer *layer, int normal0, int normal1)
{
  double direction0[3], direction1[3];

  if ( layer != layerNormalDirection(layer, normal0, direction0) ) return -1.0;
  if ( layer != layerNormalDirection(layer, normal1, direction1) ) return -1.0;
  return 
    gridConvertRadianToDegree(acos(gridDotProduct(direction0,direction1)));
}

Layer *layerNormalMinDot(Layer *layer, int normal,
			 double *mindot, double *mindir)
{
  int index;
  double dir[3], norm[3], dot;

  if (layer != layerNormalDirection(layer,normal,dir)) return NULL;

  if (normal != layer->normalTriangleHub)
    layerStoreNormalTriangleDirections(layer, normal);

  *mindot = 2.0;
  mindir[0]=0.0;
  mindir[1]=0.0;
  mindir[2]=0.0;

  for ( index =0; index < layer->normalTriangleDegree; index++ ){
    if (!layer->normalTriangleExclusive|| !layer->normalTriangleExclude[index]){
      layerNormalTriangleDirection(layer,index,norm);
      dot = norm[0]*dir[0] + norm[1]*dir[1] + norm[2]*dir[2];
      if (dot<*mindot) {
	*mindot = dot;
	mindir[0]=norm[0];
	mindir[1]=norm[1];
	mindir[2]=norm[2];
      }
    }
  }
  
  return layer;
}

Layer *layerNormalTriangles(Layer *layer, int normal, int ntriangle, int *triangles )
{
  int i;
  AdjIterator it;
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;
  i=0;
  for ( it = adjFirst(layer->triangleAdj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
    if (i>=ntriangle) return NULL;
    triangles[i] = adjItem(it);
    i++;
  }
  return layer;
}

Layer *layerStoreNormalTriangleDirections(Layer *layer, int normal)
{
  int tri, triangles[MAXNORMALDEG];
  int triangle, nodes[3], side0, side1, faceId;
  double xyz0[3], xyz1[3], tangent[3];
  double newxyz[3], uv[2], normalDirection[3];
  
  if (layer != layerNormalTriangles(layer, normal, MAXNORMALDEG, triangles )) {
    layer->normalTriangleDegree = EMPTY;
    return NULL;
  }
  layer->normalTriangleHub = normal;
  layer->normalTriangleDegree = layerNormalDeg(layer, normal );
  layer->normalTriangleExclusive = FALSE;
  for (tri=0;tri<layer->normalTriangleDegree;tri++) {
    triangle = triangles[tri];
    if ( 0 == layerConstrained(layer,normal) || 
	 ( 0 == layerConstrainedSide(layer, triangle, 0 ) &&
	   0 == layerConstrainedSide(layer, triangle, 1 ) &&
	   0 == layerConstrainedSide(layer, triangle, 2 ) ) ) {
      layer->normalTriangleExclude[tri] = TRUE;
      layerTriangleDirection( layer, triangle,
			      &layer->normalTriangleDirection[3*tri] );
    }else{
      layer->normalTriangleExclusive = TRUE;
      layer->normalTriangleExclude[tri] = FALSE;
      layerTriangle(layer,triangle,nodes);
      if (0<layerConstrainedSide(layer, triangle, 0 )) side0 = 0; 
      if (0<layerConstrainedSide(layer, triangle, 1 )) side0 = 1; 
      if (0<layerConstrainedSide(layer, triangle, 2 )) side0 = 2;
      side1 = side0+1; if (side1>2) side1 = 0;
      gridNodeXYZ(layerGrid(layer),nodes[side0],xyz0);
      gridNodeXYZ(layerGrid(layer),nodes[side1],xyz1);
      gridSubtractVector(xyz1,xyz0,tangent);
      gridVectorNormalize(tangent);
      
      faceId = layerConstrained(layer,normal);
      if (faceId>0) {
	gridProjectToFace(layerGrid(layer),faceId,xyz0,uv,newxyz);
	gridFaceNormalAtUV(layerGrid(layer), faceId,
			   uv, newxyz, normalDirection);
	gridCrossProduct(normalDirection,tangent,
			 &layer->normalTriangleDirection[3*tri]);
	gridVectorNormalize(&layer->normalTriangleDirection[3*tri]);
      } else {
	layerNormalDirection( layer, normal, 
			      &layer->normalTriangleDirection[3*tri] );
      }
    }
  }
  return layer;
}

Layer *layerNormalTriangleDirection(Layer *layer, int index, double *direction )
{
  if (index<0||index>layer->normalTriangleDegree) return NULL;
  direction[0] = layer->normalTriangleDirection[0+3*index];
  direction[1] = layer->normalTriangleDirection[1+3*index];
  direction[2] = layer->normalTriangleDirection[2+3*index];
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

  for ( it = adjFirst(layer->triangleAdj,normal); 
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

  for ( it = adjFirst(layer->triangleAdj,normal); 
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
    radian += gridPI;
  }else{
    radian = gridPI - radian;
  }
  return gridConvertRadianToDegree(radian);
}

Layer *layerNormalDirection(Layer *layer, int normal, double *direction )
{
  int i;
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  for ( i=0;i<3;i++) direction[i] = layer->normal[normal].direction[i];

  return layer;
}

GridBool layerNormalDirectionFrozen(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return FALSE;
  return layer->normal[normal].directionFrozen;
}

Layer *layerNormalDirectionFreeze(Layer *layer, int normal )
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;
  layer->normal[normal].directionFrozen = TRUE;
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

Layer *layerScaleNormalHeightWithPolynomial(Layer *layer, 
					    double constant, double slope,
					    double exponent, double scale,
					    double *origin, double *direction)
{
  int normal;
  double distance, rate;
  double distanceVector[3], normalOrigin[3];
  
  for(normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(layerGrid(layer), layerNormalRoot(layer,normal), normalOrigin);
    gridSubtractVector(normalOrigin,origin,distanceVector);
    distance = gridDotProduct(distanceVector,direction);
    if (distance >= 0.0 ){
      rate = constant + slope*pow(distance,exponent);
      rate = exp(scale*log(rate));
      layer->normal[normal].height=rate*layer->normal[normal].height;
    }
  }
  return layer;	
}

Layer *layerSetNormalMaxLength(Layer *layer, int normal, double maxLength)
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  layer->normal[normal].maxlength=maxLength;

  return layer;
}

Layer *layerSetPolynomialMaxHeight(Layer *layer, 
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
      layerSetNormalMaxLength( layer, normal, constant + slope*pow(distance,exponent) );
  }
  return layer;	
}

double layerNormalMaxLength(Layer *layer, int normal)
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return -1.0;

  return layer->normal[normal].maxlength;
}

Layer *layerSetNormalInitialHeight(Layer *layer, int normal, 
				   double initialHeight)
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  layer->normal[normal].initialheight = initialHeight;

  return layer;
}

Layer *layerSaveInitialNormalHeight(Layer *layer)
{
  int normal;

  for(normal=0;normal<layerNNormal(layer);normal++)
    layer->normal[normal].initialheight = layer->normal[normal].height;

  return layer;
}

double layerNormalInitialHeight(Layer *layer, int normal)
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return -1.0;

  return layer->normal[normal].initialheight;
}

Layer *layerSetNormalRate(Layer *layer, int normal, double rate)
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return NULL;

  layer->normal[normal].rate=rate;

  return layer;
}

double layerNormalRate(Layer *layer, int normal)
{
  if (normal < 0 || normal >= layerNNormal(layer) ) return -1.0;

  return layer->normal[normal].rate;
}

Layer *layerSetAllNormalRate(Layer *layer, double rate)
{
  int normal;
  for(normal=0;normal<layerNNormal(layer);normal++)
    layerSetNormalRate(layer, normal, rate);

  return layer;
}

Layer *layerSetNormalHeightWithRate(Layer *layer)
{
  int normal;
  double length;

  for(normal=0;normal<layerNNormal(layer);normal++){

    length = layer->normal[normal].length;

    layer->normal[normal].height =
      layer->normal[normal].initialheight
      + layer->normal[normal].rate*length;

  }

  return layer;
}

Layer *layerSetNormalHeightWithMaxRate(Layer *layer, double maxRate)
{
  int normal;
  double length, height, old;

  for(normal=0;normal<layerNNormal(layer);normal++){

    length = layer->normal[normal].length;
    
    height = layer->normal[normal].initialheight 
           + layer->normal[normal].rate*length;
    old = layer->normal[normal].height;
    layer->normal[normal].height = MIN(height, maxRate*old);
  }
  
  return layer;
}
Layer *layerSetNormalHeightForLayerNumber(Layer *layer, int n, double rate)
{
  int normal;
  for(normal=0;normal<layerNNormal(layer);normal++)
    layer->normal[normal].height = 
      layer->normal[normal].initialheight * pow(rate,n);

  return layer;
}

Layer *layerFeasibleNormals(Layer *layer, double dotLimit, double relaxation )
{
  int normal, iter;
  double *dir, mindir[3], mindot, worstdot; 

  if (dotLimit < 0) dotLimit = 1.0e-14;
  if (relaxation < 0) relaxation = 1.01;

  layerProjectNormalsToConstraints(layer);

  worstdot =2;
  for (normal=0;normal<layerNNormal(layer);normal++){
    if ( 0 != layerNormalDeg(layer, normal ) ) {
      dir = layer->normal[normal].direction;
      layerNormalMinDot(layer, normal, &mindot, mindir );
      for (iter=0;iter<1000 && mindot <= dotLimit; iter++){
	dir[0] -= relaxation*mindot*mindir[0];
	dir[1] -= relaxation*mindot*mindir[1];
	dir[2] -= relaxation*mindot*mindir[2];
	gridVectorNormalize(dir);
	layerProjectNormalToConstraints(layer,normal);
	layerNormalMinDot(layer, normal, &mindot, mindir );
      }
      worstdot = MIN(worstdot,mindot);
    }
  }

  if (worstdot>0.0) {
    return layer;
  } else {
    return NULL;
  }
}

Layer *layerVisibleNormals(Layer *layer, double dotLimit, double radianLimit )
{
  int normal, iter;
  double *dir, mindir[3], mindot, radian, worstdot; 
  double lastdir[3], lastdot;

  if (layerNNormal(layer) == 0 ) return NULL;
  if (dotLimit < 0) dotLimit = 0.90;
  if (radianLimit < 0) radianLimit = 1.0e-15;

  layerProjectNormalsToConstraints(layer);

  worstdot = 2.0;
  for (normal=0;normal<layerNNormal(layer);normal++){
    if ( 0 != layerNormalDeg(layer, normal ) ) {
      dir = layer->normal[normal].direction;
      radian = 0.01;
      layerNormalMinDot(layer, normal, &mindot, mindir );
      for (iter=0;
	   iter<1000 && (radian>radianLimit || mindot<0.05) && mindot<dotLimit;
	   iter++){
	lastdir[0] = dir[0]; lastdir[1] = dir[1]; lastdir[2] = dir[2];
	dir[0] += radian*mindir[0];
	dir[1] += radian*mindir[1];
	dir[2] += radian*mindir[2];
	gridVectorNormalize(dir);
	layerProjectNormalToConstraints(layer,normal);
	lastdot = mindot;
	layerNormalMinDot(layer, normal, &mindot, mindir );
	if (mindot <= lastdot) {
	  radian *= 0.8;
	  dir[0] = lastdir[0]; dir[1] = lastdir[1]; dir[2] = lastdir[2];
	  layerNormalMinDot(layer, normal, &mindot, mindir );
	}
      }
      worstdot = MIN(worstdot,mindot);
      if (mindot <= 0.0 ) {
	double xyz[3];
	gridNodeXYZ(layerGrid(layer),layerNormalRoot(layer,normal),xyz);
	printf("ERROR: %s,%5d, Invisible norm%6d dot%14.10f X%9.5f Y%9.5f Z%9.5f\n",
	       __FILE__, __LINE__, normal, mindot,xyz[0],xyz[1],xyz[2]);
      }
    }
  }

  if (worstdot>0.0) {
    return layer;
  } else {
    return NULL;
  }
}

Layer *layerSmoothNormalDirection(Layer *layer, double relax )
{
  int normal, iter, triangle, normals[3], total, i;
  double norm[3], avgdir[3], denom, relaxm1; 
  double visibility = 0.5;
  double visTol = 1.0e-10;
  AdjIterator it;

  if (layerNNormal(layer) == 0 ) return NULL;
  if (layerNBlend(layer) != 0 ) return NULL;

  relaxm1 = 1.0-relax;

  layerProjectNormalsToConstraints(layer);

  for (iter=0;iter<1;iter++){
    for (normal=0;normal<layerNNormal(layer);normal++){
      if ( 0 < layerConstrained(layer,normal) ){
	total = 0;
	avgdir[0]=0.0;
	avgdir[1]=0.0;
	avgdir[2]=0.0;
	for ( it = adjFirst(layer->triangleAdj,normal); 
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
	layer->normal[normal].direction[0] = 
	  relax*(avgdir[0] * denom) + 
	  relaxm1*layer->normal[normal].direction[0];
	layer->normal[normal].direction[1] = 
	  relax*(avgdir[1] * denom) +
	  relaxm1*layer->normal[normal].direction[1];
	layer->normal[normal].direction[2] = 
	  relax*(avgdir[2] * denom) +
	  relaxm1*layer->normal[normal].direction[2];
	gridVectorNormalize(layer->normal[normal].direction);
	layerProjectNormalToConstraints(layer,normal);
      }
    }
    layerVisibleNormals(layer,visibility,visTol);
    for (normal=0;normal<layerNNormal(layer);normal++){
      if ( 0 == layerConstrained(layer,normal) ){
	total = 0;
	avgdir[0]=0.0;
	avgdir[1]=0.0;
	avgdir[2]=0.0;
	for ( it = adjFirst(layer->triangleAdj,normal); 
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
	layer->normal[normal].direction[0] = 
	  relax*(avgdir[0] * denom) + 
	  relaxm1*layer->normal[normal].direction[0];
	layer->normal[normal].direction[1] = 
	  relax*(avgdir[1] * denom) +
	  relaxm1*layer->normal[normal].direction[1];
	layer->normal[normal].direction[2] = 
	  relax*(avgdir[2] * denom) +
	  relaxm1*layer->normal[normal].direction[2];
	gridVectorNormalize(layer->normal[normal].direction);
      }
    }
    layerVisibleNormals(layer,visibility,visTol);
  }

  return layer;
}

Layer *layerSmoothInteriorNormalDirection(Layer *layer, 
					  double relax, int iterations,
					  double visibility)
{
  int normal, iter, triangle, normals[3], total, i;
  double norm[3], avgdir[3], denom;
  double relaxm1; 
  double mindir[3], mindot;
  double visTol = 1.0e-5;
  AdjIterator it;

  if (layerNNormal(layer) == 0 ) return NULL;
  if (layerNBlend(layer) != 0 ) return NULL;

  if (relax < 0.0) relax = 0.5;
  if (iterations<0) iterations = 20;
  if (visibility<0.0) visibility = 0.1;

  relaxm1 = 1.0-relax;

  layerProjectNormalsToConstraints(layer);

  for (normal=0;normal<layerNNormal(layer);normal++){
    layerNormalMinDot(layer, normal, &mindot, mindir );
    if (mindot < 0.5) layerNormalDirectionFreeze(layer,normal);
  }

  for (iter=0;iter<iterations;iter++){
    printf("Normal smoothing iteration %d.\n",iter);
    for (normal=0;normal<layerNNormal(layer);normal++){
      if (layerNormalDirectionFrozen(layer,normal)) continue;
      if ( 0 < layerConstrained(layer,normal) ){
	total = 0;
	avgdir[0]=0.0;
	avgdir[1]=0.0;
	avgdir[2]=0.0;
	for ( it = adjFirst(layer->triangleAdj,normal); 
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
	layer->normal[normal].direction[0] = 
	  relax*(avgdir[0] * denom) + 
	  relaxm1*layer->normal[normal].direction[0];
	layer->normal[normal].direction[1] = 
	  relax*(avgdir[1] * denom) +
	  relaxm1*layer->normal[normal].direction[1];
	layer->normal[normal].direction[2] = 
	  relax*(avgdir[2] * denom) +
	  relaxm1*layer->normal[normal].direction[2];
	gridVectorNormalize(layer->normal[normal].direction);
	layerProjectNormalToConstraints(layer,normal);
      }
    }
    layerVisibleNormals(layer,visibility,visTol);
    for (normal=0;normal<layerNNormal(layer);normal++){
      if (layerNormalDirectionFrozen(layer,normal)) continue;
      if ( 0 == layerConstrained(layer,normal) ){
	total = 0;
	avgdir[0]=0.0;
	avgdir[1]=0.0;
	avgdir[2]=0.0;
	for ( it = adjFirst(layer->triangleAdj,normal); 
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
	layer->normal[normal].direction[0] = 
	  relax*(avgdir[0] * denom) + 
	  relaxm1*layer->normal[normal].direction[0];
	layer->normal[normal].direction[1] = 
	  relax*(avgdir[1] * denom) +
	  relaxm1*layer->normal[normal].direction[1];
	layer->normal[normal].direction[2] = 
	  relax*(avgdir[2] * denom) +
	  relaxm1*layer->normal[normal].direction[2];
	gridVectorNormalize(layer->normal[normal].direction);
      }
    }
    layerVisibleNormals(layer,visibility,visTol);
  }

  return layer;
}

Layer *layerProjectNormalsToConstraints(Layer *layer)
{
  int normal;
  for (normal=0;normal<layerNNormal(layer);normal++){
    if (layer!=layerProjectNormalToConstraints(layer, normal)) return NULL;
  }
  return layer;
}

Layer *layerProjectNormalToConstraints(Layer *layer, int normal)
{
  int tipnode, edgeId, faceId;
  double xyzroot[3], xyztip[3], direction[3], height;
  Grid *grid;

  faceId = layerConstrained(layer,normal);
  if (faceId==0) return layer;

  grid = layerGrid(layer);

  tipnode = gridAddNode(grid,0.0,0.0,0.0);
  if (EMPTY == tipnode) return NULL;

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

  gridRemoveNode(grid,tipnode);

  return layer;
}

Layer *layerAdjustNormalHeightToSmoothFront(Layer *layer, double maxHeight)
{
  int normal, node, nodes[3];
  double direction[3], xyz0[3], xyz1[3], tangent[3];
  double dn, height;
  AdjIterator it;
  int triangle;
  int i,n;
  Grid *grid;

  if (layerNBlend(layer) != 0 ) return NULL;  

  grid = layerGrid(layer);

  for ( normal = 0 ; normal < layerNNormal(layer) ; normal++ ) {
    layerNormalDirection(layer, normal, direction);
    node = layerNormalRoot(layer,normal);
    gridNodeXYZ(grid,node,xyz0);
    dn = 0;
    n = 0;
    for ( it = adjFirst(layer->triangleAdj,normal); 
	  adjValid(it); 
	  it = adjNext(it) ){
      triangle = adjItem(it);
      layerTriangle(layer, triangle, nodes);
      for(i=0;i<2;i++){
	if(node != nodes[i]){
	  n++;
	  gridNodeXYZ(grid,nodes[i],xyz1);
	  gridSubtractVector(xyz1,xyz0,tangent);
	  gridVectorNormalize(tangent);
	  dn += gridDotProduct(direction, tangent);
	}
      }
    }
    dn /= (double)n;
    layerGetNormalHeight(layer, normal, &height);
    height += dn*maxHeight*height;
    layerSetNormalHeight(layer, normal, height);
  }

  return layer;
}

Layer *layerConstrainNormal(Layer *layer, int edgeface )
{
  int edgeId, faceId, face, nodes[3], id, i, normal;
  int n0, n1, normal0, normal1;
  int nCurveNode, *curve;

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
    if (nCurveNode>0) {
      curve = malloc( nCurveNode * sizeof(int) );
      gridGeomEdge( layer->grid, edgeId, curve );
      for ( i=0; i<nCurveNode; i++){ 
	normal = layer->globalNode2Normal[curve[i]];
	if (normal != EMPTY) layer->normal[normal].constrained=-edgeId;
      }
      free(curve);
    }
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

GridBool layerConstrainingGeometry(Layer *layer, int edgeface )
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

  for ( it = adjFirst(layer->triangleAdj,normal0); 
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

  for ( it = adjFirst(layer->triangleAdj,normal0); 
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

GridBool layerNormalTerminated(Layer *layer, int normal )
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

Layer *layerTerminateTriangleNormals(Layer *layer, int triangle ){
  int normals[3];

  if ( NULL == layerTriangleNormals(layer,triangle,normals) ) return(NULL);
  layerTerminateNormal( layer, normals[0] );
  layerTerminateNormal( layer, normals[1] );
  layerTerminateNormal( layer, normals[2] );
	
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

GridBool layerAnyActiveNormals(Layer *layer)
{
  return (layerNActiveNormal(layer)>0);
}

GridBool layerCellInLayer(Layer *layer, int cell)
{
  if (cell < 0 || cell >= gridMaxCell(layerGrid(layer)) ) return FALSE;
  return layer->cellInLayer[cell];
}

GridBool layerFaceInLayer(Layer *layer, int face)
{
  if (face < 0 || face >= gridMaxFace(layerGrid(layer)) ) return FALSE;
  return layer->faceInLayer[face];
}

GridBool layerEdgeInLayer(Layer *layer, int edge)
{
  if (edge < 0 || edge >= gridMaxEdge(layerGrid(layer)) ) return FALSE;
  return layer->edgeInLayer[edge];
}

int layerNEdgeInLayer(Layer *layer, int edgeId)
{
  int edge, count;
  int nodes[2], currentEdgeId;
  Grid *grid;

  grid = layerGrid(layer);
  count = 0;

  for(edge=0;edge<gridMaxEdge(grid);edge++)
    if ( ( grid == gridEdge(grid, edge, nodes, &currentEdgeId ) ) &&
	 ( edgeId == currentEdgeId ) &&
	 ( layerEdgeInLayer(layer, edge) ) ) count++;
	 
  return count;
}

int layerEdgeEndPoint(Layer *layer, int edgeId, int startNode)
{
  AdjIterator it;
  int node, lastnode, edge, n1;
  int nodes[2], currentEdgeId;
  GridBool found;
  Grid *grid;

  grid = layerGrid(layer);

  node = startNode;
  lastnode = EMPTY;
  found = TRUE;
  while (found) {
    found = FALSE;
    for ( it = adjFirst(gridEdgeAdj(grid),node); 
	  adjValid(it) && !found; 
	  it = adjNext(it)) {
      edge = adjItem(it);
      gridEdge(grid, edge, nodes, &currentEdgeId);
      if ( ( currentEdgeId == edgeId ) &&
	   ( layerEdgeInLayer(layer,edge) ) ) {
	if ( node == nodes[0] ) {
	  n1 = nodes[1];
	}else{
	  n1 = nodes[0];	  
	}
	if ( n1 != lastnode ) { 
	  found = TRUE;
	  lastnode = node;
	  node = n1;
	}
      }
    }
  }

  return node;
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
  return layerAdvance(layer, TRUE);
}

#define addTet \
{ \
  layer->cellInLayer[gridAddCell(grid, tet[0], tet[1], tet[2], tet[3])]=TRUE; \
  if ( 1.0e-14 > gridVolume(grid, tet ) ) { \
    negVolume = TRUE; \
    gridNodeXYZ(grid,tet[0],xyz); \
    printf("volume%18.10e at%15.6f%15.6f%15.6f\n", \
           gridVolume(grid, tet ),xyz[0],xyz[1],xyz[2]); \
    gridWriteTecplotCellZone(grid,tet,"layerNegVolCell.t"); \
  } \
}

Layer *layerAdvanceBlends(Layer *layer)
{
  int blend, blendnormals[4];
  int subBlend;
  int triangle0, triangle1;
  int normal, faceId;
  int tet[4], n[6];
  double xyz[3];
  double minVolume0, minVolume1;
  GridBool negVolume = FALSE;
  Grid *grid = layer->grid;
  if (!layerTetrahedraOnly(layer)) {
    printf("ERROR: %s: %d: Using blends requires layerTetrahedraOnly %s.\n",
	   __FILE__,__LINE__,"(no mixed elements)");
    return NULL;
  }
  
  for (blend=0;blend<layerNBlend(layer);blend++){
    for(subBlend=0; subBlend< layerNSubBlend(layer,blend); subBlend++){
      layerSubBlendNormals(layer, blend, subBlend, blendnormals );
      
      triangle0 = EMPTY;
      triangle1 = EMPTY;
      if (blendnormals[0] != blendnormals[1]) 
	triangle0 = layerForceTriangle(layer,blendnormals[0],
				       blendnormals[1],blendnormals[2]);
      if (blendnormals[2] != blendnormals[3]) 
	triangle1 = layerForceTriangle(layer,blendnormals[1],
				       blendnormals[3],blendnormals[2]);
      
      if ( EMPTY != triangle0) {
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
      }
      if ( EMPTY != triangle1) {
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

      minVolume0 = DBL_MAX;
      if (n[4]!=n[5]) {
	tet[0] = n[0]; tet[1] = n[4]; tet[2] = n[5]; tet[3] = n[3]; 
	minVolume0=MIN(minVolume0,gridVolume(grid, tet ));
      }
      if (n[4]!=n[5]) {
	tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[4]; tet[3] = n[5];
	minVolume0=MIN(minVolume0,gridVolume(grid, tet ));
      }
      if (n[1]!=n[2]) {
	tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[1]; tet[3] = n[4];
	minVolume0=MIN(minVolume0,gridVolume(grid, tet ));
      }
      minVolume1 = DBL_MAX;
      if (n[4]!=n[5]) {
	tet[0] = n[0]; tet[1] = n[4]; tet[2] = n[5]; tet[3] = n[3];
	minVolume1=MIN(minVolume1,gridVolume(grid, tet ));
      }
      if (n[4]!=n[5]) {
	tet[0] = n[0]; tet[1] = n[1]; tet[2] = n[5]; tet[3] = n[4];
	minVolume1=MIN(minVolume1,gridVolume(grid, tet ));
      }
      if (n[1]!=n[2]) {
	tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[1]; tet[3] = n[5];
	minVolume1=MIN(minVolume1,gridVolume(grid, tet ));
      }

      if (minVolume0>minVolume1) {
	if (n[4]!=n[5]) {
	  tet[0] = n[0]; tet[1] = n[4]; tet[2] = n[5]; tet[3] = n[3]; addTet;
	}
	if (n[4]!=n[5]) {
	  tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[4]; tet[3] = n[5]; addTet;
	}
	if (n[1]!=n[2]) {
	  tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[1]; tet[3] = n[4]; addTet;
	}
      }else{
	if (n[4]!=n[5]) {
	  tet[0] = n[0]; tet[1] = n[4]; tet[2] = n[5]; tet[3] = n[3]; addTet;
	}
	if (n[4]!=n[5]) {
	  tet[0] = n[0]; tet[1] = n[1]; tet[2] = n[5]; tet[3] = n[4]; addTet;
	}
	if (n[1]!=n[2]) {
	  tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[1]; tet[3] = n[5]; addTet;
	}
      }
    }
  }
 
  for ( normal = 0 ; normal < adjNNode(layer->blendAdj) ; normal++ ) {
    int sweep;
    int *allVertexNormals, vertexNormals[3], nVertexNormals;
    switch (layerBlendDegree(layer,normal)) {
    case 0: case 1: case 2: break;
    case 3:
      nVertexNormals = layerSubNormalDegree(layer,normal);
      allVertexNormals = malloc(nVertexNormals*sizeof(int));
      layerOrderedVertexNormals( layer, normal, 
				 &nVertexNormals, allVertexNormals);
      for(sweep=0;sweep<nVertexNormals;sweep++) {
	vertexNormals[0] = layer->vertexNormal[normal];
	vertexNormals[1] = allVertexNormals[sweep];
	if (sweep < (nVertexNormals-1) ) {
	  vertexNormals[2] = allVertexNormals[sweep+1];
	}else{
	  vertexNormals[2] = allVertexNormals[0];
	}
	triangle0 = layerForceTriangle(layer,vertexNormals[0],
				       vertexNormals[1],vertexNormals[2]);
	/*face oriented opposite direction to cell*/ 
	tet[0] = layerNormalTip(layer,vertexNormals[1]);
	tet[1] = layerNormalTip(layer,vertexNormals[0]); 
	tet[2] = layerNormalTip(layer,vertexNormals[2]); 
	tet[3] = layerNormalRoot(layer,normal); 
	addTet;
      }
      free(allVertexNormals);
      break;
    default:
      printf( "ERROR: %s: %d: Cannot handle %d blends. Write more code!\n",
	      __FILE__, __LINE__, layerBlendDegree(layer,normal));
      break;
    }
  }
  
  layerBuildNormalTriangleAdjacency(layer);
  layer->nblend=0;
  layer->maxblend=0;
  free(layer->blend); layer->blend = NULL;
  adjFree(layer->blendAdj); layer->blendAdj = NULL;
  free(layer->vertexNormal); layer->vertexNormal = NULL;

  if ( negVolume ) {
    return NULL;
  } else {
    return layer;
  }
}

Layer *layerAdvance(Layer *layer, GridBool reconnect)
{
  Grid *grid = layer->grid;
  int normal, normal0, normal1, root, tip, faceId, edgeId, i;
  int node;
  int triangle, normals[3], nodes[3], n[6];
  int side[2], sidenode[2], sidenormal[2];
  double xyz[3];
  int nterminated;
  int tet[4];
  GridBool negVolume = FALSE;

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
      layer->normal[normal].length += layer->normal[normal].height;
      tip = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);
      if ( EMPTY == tip) return NULL;
      layer->normal[normal].tip = tip;
      if (reconnect) layerReconnectCellUnlessInLayer(layer, root, tip);
      faceId = layerConstrained(layer,normal);
      if (reconnect && 0 > faceId) {
	edgeId = -faceId;
	layerReconnectEdgeUnlessInLayer(layer, edgeId, root, tip);
      }
      gridCopySpacing(grid, root, tip );
      gridFreezeNode( grid, tip );
    }
    linesAddNode(gridLines(grid),normal,layer->normal[normal].root);
    linesAddNode(gridLines(grid),normal,layer->normal[normal].tip);
  }

  /* reconnect faces for constrained triangleside */
  if (reconnect) {
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
  }
  
  /* advance edges */
  for (normal=0;normal<layerNNormal(layer);normal++) {
    edgeId = -layerConstrained(layer,normal);
    if (edgeId > 0) {
      root = layer->normal[normal].root;
      tip  = layer->normal[normal].tip;
      /* note that tip has been set to root on terminated normals */
      if (root != tip) layer->edgeInLayer[gridAddEdge(grid,
						      root,tip,edgeId,
						      DBL_MAX,DBL_MAX)] = TRUE;
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

    /* advance cells */
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
	if (n[0]!=n[3]) {
	  tet[0] = n[0]; tet[1] = n[4]; tet[2] = n[5]; tet[3] = n[3]; addTet;
	}
	if (n[2]!=n[5]) {
	  tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[4]; tet[3] = n[5]; addTet;
	}
	if (n[1]!=n[4]) {
	  tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[1]; tet[3] = n[4]; addTet;
	}
      }else{
	if (n[0]!=n[3]) {
	  tet[0] = n[0]; tet[1] = n[4]; tet[2] = n[5]; tet[3] = n[3]; addTet;
	}
	if (n[1]!=n[4]) {
	  tet[0] = n[0]; tet[1] = n[1]; tet[2] = n[5]; tet[3] = n[4]; addTet;
	}
	if (n[2]!=n[5]) {
	  tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[1]; tet[3] = n[5]; addTet;
	}
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

  if (layerNBlend(layer) > 0) {
    negVolume = negVolume || (layer != layerAdvanceBlends(layer));
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

  if ( negVolume ) {
    return NULL;
  } else {
    return layer;
  }
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
  printf("normals %d of %d terminated by spacing criterion. %d active.\n",
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
      if (xyz[0]>x) { layerTerminateNormal(layer, normal); nterm++;  }
    }else{
      if (xyz[0]<x) { layerTerminateNormal(layer, normal); nterm++; }     
    }
  }
  printf("normals %d of %d terminated by X criterion.\n",nterm,layerNNormal(layer) );
  return layer;
}

int layerTerminateNormalWithLength(Layer *layer, double ratio)
{
  int normal, nterm, totalterm;
  
  if (layerNNormal(layer) == 0 ) return 0;

  nterm = 0;
  for (normal=0;normal<layerNNormal(layer);normal++){
    if ( layer->normal[normal].length + layer->normal[normal].height > 
	 (ratio * layer->normal[normal].maxlength) ) {
      if( !layerNormalTerminated(layer, normal ) ) nterm++; 
      layerTerminateNormal(layer, normal); 
    }
  }
  totalterm = layerNNormal(layer)-layerNActiveNormal(layer);
  printf("%d terminated by Length test; %d of %d normals terminted.\n",
	 nterm, totalterm,layerNNormal(layer) );
  return totalterm;
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


GridBool layerTetrahedraOnly(Layer *layer)
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

  adj = adjCreate( layerNNormal(layer), layerNBlend(layer)*4, 1000 );

  for (blend=0;blend<layerNBlend(layer);blend++){
    layerBlendNormals(layer,blend,normals);
    for (i=0;i<4;i++) adjRegister( adj, normals[i], blend );
  }

  return adj;
}

int layerFirstTriangleAfterGap(Layer *layer, int normal )
{
  int firstTriangle, triangle, previousTriangle;
  AdjIterator it;
  
  firstTriangle = adjItem(adjFirst(layer->triangleAdj,normal));
  for ( it = adjFirst(layer->triangleAdj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
    triangle = adjItem(it);
    previousTriangle = layerPreviousTriangle(layer, normal, triangle);
    if (EMPTY == previousTriangle) firstTriangle = triangle;
  }
  return firstTriangle;
}

int layerNRequiredBlends(Layer *layer, int normal, double angleLimit )
{
  int blendCount;
  int triangle, nextTriangle;
  double edgeAngle;
  AdjIterator it;

  if (angleLimit < 0.0) angleLimit = 250; /* deg */

  blendCount = 0;
  for ( it = adjFirst(layer->triangleAdj,normal); 
	adjValid(it); 
	it = adjNext(it) ){
    triangle = adjItem(it);
    nextTriangle = layerNextTriangle(layer, normal, triangle);
    if (EMPTY != nextTriangle) {
      edgeAngle = layerEdgeAngle(layer,triangle,nextTriangle);
      if (edgeAngle > angleLimit) blendCount++;
    }
  }

  return blendCount;
}

Layer *layerBlend(Layer *layer, double angleLimit )
{
  int normal;
  int triangle, firstTriangle, nextTriangle;
  double edgeAngle;
  int commonEdge[2];
  int nRequiredBlends;
  GridBool continuous;

  int count, i;
  int leftNormal, rightNormal;
  int edge;

  if (angleLimit < 0.0) angleLimit = 250; /* deg */

  if (NULL != layer->vertexNormal) free(layer->vertexNormal);
  layer->vertexNormal = malloc(layer->maxnormal*sizeof(int));
  for (i=0;i<layer->maxnormal;i++) layer->vertexNormal[i]=EMPTY;

  layer->originalnormal = layerNNormal(layer);
  if (layer->blendAdj != NULL) adjFree(layer->blendAdj);
  layer->blendAdj = adjCreate( layer->originalnormal, 
			       layer->originalnormal, 1000);

  for ( normal = 0 ; normal < layer->originalnormal ; normal++ ) {
    nRequiredBlends = layerNRequiredBlends(layer, normal, angleLimit );
    if ( nRequiredBlends > 2 ) 
      layer->vertexNormal[normal] = layerDuplicateNormal(layer,normal);
    if ( nRequiredBlends > 0 ) {
      firstTriangle = layerFirstTriangleAfterGap(layer,normal);
      continuous = ( EMPTY != 
		     layerPreviousTriangle(layer, normal, firstTriangle) );
      leftNormal = normal;
      rightNormal = leftNormal;
      triangle = firstTriangle;
      nextTriangle = layerNextTriangle(layer, normal, triangle);
      edge = 0;
      for ( count = 0; count < layerNormalDeg(layer,normal); count++ ){
	if (EMPTY != nextTriangle) {
	  edgeAngle = layerEdgeAngle(layer,triangle,nextTriangle);
	  if (edgeAngle > angleLimit) {
	    edge = edge +1;
	    if (continuous && 1 == nRequiredBlends) {
	    }else{
	      if (continuous && edge == nRequiredBlends) {
		leftNormal = normal;
	      }else{
		leftNormal = layerDuplicateNormal(layer, normal );
	      }
	    }
	    layerCommonEdge(layer, triangle, nextTriangle, commonEdge);
	    if (layerNormalRoot(layer,normal) == commonEdge[0] ) {
	      adjRegister( layerBlendAdj(layer), normal,
			   layerAddBlend( layer,leftNormal,rightNormal,
					  commonEdge[1]) );
	    }else{
	      adjRegister( layerBlendAdj(layer), normal,
			   layerAddBlend( layer,leftNormal,rightNormal,
					  commonEdge[0]) );
	    }
	  }
	  for (i=0;i<3;i++)
	    if (layer->triangle[nextTriangle].normal[i] == normal) 
	      layer->triangle[nextTriangle].normal[i] = leftNormal;
	}
	triangle = nextTriangle;
	nextTriangle = layerNextTriangle(layer, normal, triangle);      
	rightNormal = leftNormal;
      }
      
    }
 
  }

  layerBuildNormalTriangleAdjacency(layer);
  layerInitializeTriangleNormalDirection(layer);
  layerFeasibleNormals(layer, -1.0, -1.0 );
  layerVisibleNormals(layer, 0.3, 1.0e-8 );
  return layer;
}

int layerAddBlend(Layer *layer, int normal0, int normal1, int otherNode )
{
  int i, node0, node1, n0, n1;
  GridBool newEdge;
  int blend;
  int triangle, edge, nodes[2], excludeId, edgeId, side;
  AdjIterator it;
  Grid *grid;
  grid=layerGrid(layer);

  if (NULL == layer->blend) {
    layer->maxblend = 5000;
    layer->blend = malloc(layer->maxblend*sizeof(Blend));
  }

  if (layer->nblend >= layer->maxblend) {
    layer->maxblend += 5000;
    layer->blend = realloc(layer->blend,layer->maxblend*sizeof(Blend));
  }
  
  n0 = layerNormalRoot(layer,normal0);
  n1 = otherNode;

  blend   = EMPTY;
  newEdge = TRUE;
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
    layer->blend[blend].normal[2] = EMPTY;
    layer->blend[blend].normal[3] = EMPTY;
    layer->blend[blend].edgeId[0] = EMPTY;
    layer->blend[blend].edgeId[1] = EMPTY;
    layer->blend[blend].nSubNormal0 = 0;
    layer->blend[blend].nSubNormal1 = 0;
    layer->nblend++;
  }


  if ( layerConstrained(layer,normal0) != 0 ) {
    edge = gridFindEdge(grid,n0,n1);
    excludeId=0;
    gridEdge(grid,edge,nodes,&excludeId);
  
    side = newEdge?0:1;

    for ( it = adjFirst(layer->triangleAdj,normal0); 
	  adjValid(it); 
	  it=adjNext(it) ){
      triangle = adjItem(it);
      for (i=0;i<3;i++) {
	edgeId = layer->triangle[triangle].parentGeomEdge[i];
	if (edgeId > 0 && edgeId != excludeId) 
	  layer->blend[blend].edgeId[side] = edgeId;
      }
    }
    for ( it = adjFirst(layer->triangleAdj,normal1);
	  adjValid(it);
	  it=adjNext(it) ){
      triangle = adjItem(it);
      for (i=0;i<3;i++) {
	edgeId = layer->triangle[triangle].parentGeomEdge[i];
	if (edgeId > 0 && edgeId != excludeId) 
	  layer->blend[blend].edgeId[side] = edgeId;
      }
    }
  }

  return blend;
}

Layer *layerBlendNormals(Layer *layer, int blend, int *normals )
{
  int i;
  if (blend < 0 || blend >= layerNBlend(layer)) return NULL;
  for(i=0;i<4;i++) normals[i] = layer->blend[blend].normal[i];
  return layer;
}

Layer *layerSubBlendNormals(Layer *layer, int blend, int subBlend, int *normals)
{
  int subnormal;
  if (subBlend <0 || subBlend >= layerNSubBlend(layer,blend)) return NULL;
  if (layer != layerBlendNormals(layer, blend, normals)) return NULL;

  subnormal = subBlend - 1;
  if (subnormal>=0 && subnormal<layer->blend[blend].nSubNormal0 ) 
    normals[0] = layer->blend[blend].subNormal0[subnormal];
  if (subnormal>=layer->blend[blend].nSubNormal0) normals[0] = normals[1];

  subnormal = subBlend;
  if (subnormal<layer->blend[blend].nSubNormal0) 
    normals[1] = layer->blend[blend].subNormal0[subnormal];

  subnormal = subBlend - 1;
  if (subnormal>=0 && subnormal<layer->blend[blend].nSubNormal1 ) 
    normals[2] = layer->blend[blend].subNormal1[subnormal];
  if (subnormal>=layer->blend[blend].nSubNormal1) normals[2] = normals[3];

  subnormal = subBlend;
  if (subnormal<layer->blend[blend].nSubNormal1) 
    normals[3] = layer->blend[blend].subNormal1[subnormal];

  return layer;
}

Layer *layerBlendAxle(Layer *layer, int blend, double *axle)
{
  double xyz0[3], xyz1[3];
  Grid *grid = layerGrid(layer);
  
  if (blend < 0 || blend >= layerNBlend(layer)) return NULL; 

  if (grid != gridNodeXYZ(grid,layer->blend[blend].nodes[0],xyz0)) return NULL;
  if (grid != gridNodeXYZ(grid,layer->blend[blend].nodes[1],xyz1)) return NULL;

  gridSubtractVector(xyz1,xyz0,axle);
  gridVectorNormalize(axle);
  return layer;
}

Layer *layerNormalBlendAxle(Layer *layer, int normal, double *axle)
{
  int blend, blend0, blend1;
  double axle0[3], axle1[3];
  double uv[2], xyz[3];
  int faceId, node;

  if ( 0 < layerConstrained(layer,normal) ) {
    Grid *grid = layerGrid(layer);
    faceId = layerConstrained(layer,normal);
    node = layerNormalRoot(layer,normal);
    gridNodeUV(grid, node, faceId, uv);
    gridFaceNormalAtUV(grid,faceId,uv,xyz,axle);
    blend = adjItem(adjFirst(layerBlendAdj(layer),normal));
    if (node == layer->blend[blend].nodes[1]) {
      axle[0] = -axle[0]; axle[1] = -axle[1]; axle[2] = -axle[2];
    }
    return layer;
  }
  switch (layerBlendDegree(layer,normal)) {
  case 0: break;
  case 1:
    blend = adjItem(adjFirst(layerBlendAdj(layer),normal));
    layerBlendAxle(layer, blend, axle);
    break;
  case 2:
    blend0 = adjItem(adjFirst(layerBlendAdj(layer),normal));
    blend1 = adjItem(adjNext(adjFirst(layerBlendAdj(layer),normal)));
    layerBlendAxle(layer, blend0, axle0);
    layerBlendAxle(layer, blend1, axle1);
    if ( layer->blend[blend0].nodes[0] == layer->blend[blend1].nodes[1] ||
	 layer->blend[blend0].nodes[1] == layer->blend[blend1].nodes[0] ) {
      axle[0] = axle0[0] + axle1[0];
      axle[1] = axle0[1] + axle1[1];
      axle[2] = axle0[2] + axle1[2];
    }else{
      axle[0] = axle0[0] - axle1[0];
      axle[1] = axle0[1] - axle1[1];
      axle[2] = axle0[2] - axle1[2];
    }
    gridVectorNormalize(axle);
    break;
  default:
    printf( "ERROR: %s: %d: Cannot handle %d blends. use layerBlendAxle!\n",
	    __FILE__, __LINE__, layerBlendDegree(layer,normal));
    break;
  }
  return layer;
}
int layerBlendDegree(Layer *layer, int normal)
{
  if (NULL == layerBlendAdj(layer)) return 0;
  return adjDegree(layerBlendAdj(layer),normal);
}

int layerSubNormalDegree(Layer *layer, int normal)
{
  int degree;
  int blend, blendnormals[4];
  AdjIterator it;
  if (NULL == layerBlendAdj(layer)) return 0;
  degree =0;
  for ( it = adjFirst(layerBlendAdj(layer),normal); 
	adjValid(it); 
	it=adjNext(it) ){
    blend = adjItem(it);
    layerBlendNormals(layer, blend, blendnormals );
    degree++;
    if ( layerNormalRoot(layer, normal) == 
	 layerNormalRoot(layer, blendnormals[0]) ) {
      degree += layer->blend[blend].nSubNormal0;
    }else{
      degree += layer->blend[blend].nSubNormal1;      
    }
  }
  return degree;
}

Layer *layerSubBlendCount(Layer *layer, double maxNormalAngle)
{
  int normal;
  AdjIterator it;
  int blend;
  int blendnormals[4];
  double angle;
  int nSubNormal;

  if (layerNBlend(layer) <= 0) return layer;

  for (normal=0;normal<layer->originalnormal;normal++){
    switch (layerBlendDegree(layer,normal)) {
    case 0: break;
    case 1: case 2:
      blend = adjItem(adjFirst(layerBlendAdj(layer),normal));
      layerBlendNormals(layer, blend, blendnormals );
      if (normal == blendnormals[0] || normal == blendnormals[1] ) {
	angle = layerNormalAngle(layer,blendnormals[0], blendnormals[1]);
	nSubNormal = (int)(angle/maxNormalAngle)-1;
	nSubNormal = MIN(nSubNormal,MAXSUBNORMAL);
	layer->blend[blend].nSubNormal0 = nSubNormal;
      }else{
	angle = layerNormalAngle(layer,blendnormals[2], blendnormals[3]);
	nSubNormal = (int)(angle/maxNormalAngle)-1;
	nSubNormal = MIN(nSubNormal,MAXSUBNORMAL);
	layer->blend[blend].nSubNormal1 = nSubNormal;
      }
      break;
    case 3:
      for ( it = adjFirst(layerBlendAdj(layer),normal); 
	    adjValid(it); 
	    it=adjNext(it) ){
	blend = adjItem(it);
	layerBlendNormals(layer, blend, blendnormals );
	if ( layerNormalRoot(layer, normal) == 
	     layerNormalRoot(layer, blendnormals[0]) ) {
	  angle = layerNormalAngle(layer,blendnormals[0], blendnormals[1]);
	  nSubNormal = (int)(angle/maxNormalAngle)-1;
	  nSubNormal = MIN(nSubNormal,MAXSUBNORMAL);
	  layer->blend[blend].nSubNormal0 = nSubNormal;
	}else{
	  angle = layerNormalAngle(layer,blendnormals[2], blendnormals[3]);
	  nSubNormal = (int)(angle/maxNormalAngle)-1;
	  nSubNormal = MIN(nSubNormal,MAXSUBNORMAL);
	  layer->blend[blend].nSubNormal1 = nSubNormal;
	}
      }
      break;
    default:
      printf( "ERROR: %s: %d: Cannot handle %d blends. Write more code!\n",
	      __FILE__, __LINE__, layerBlendDegree(layer,normal));
      break;
    }
    
  }
  return layer;
}

Layer *layerSubBlendSmooth(Layer *layer)
{
  int normal;
  int blend;
  int blendnormals[4];
  int nSubNormal, origSubNormal;
  GridBool allSmooth;
  int slope = 2;

  if (layerNBlend(layer) <= 0) return layer;

  allSmooth = FALSE;
  while (!allSmooth) {
    allSmooth = TRUE;
    
    for (normal=0;normal<layer->originalnormal;normal++){
      if (2 == layerBlendDegree(layer,normal)) {
	blend = adjItem(adjFirst(layerBlendAdj(layer),normal));	
	layerBlendNormals(layer, blend, blendnormals );
	if (normal == blendnormals[0] || normal == blendnormals[1] ) {
	  origSubNormal = layer->blend[blend].nSubNormal0;
	  nSubNormal = layer->blend[blend].nSubNormal0;
	  nSubNormal = MIN(nSubNormal, layer->blend[blend].nSubNormal1+slope);
	}else{
	  origSubNormal = layer->blend[blend].nSubNormal1;
	  nSubNormal = layer->blend[blend].nSubNormal1;
	  nSubNormal = MIN(nSubNormal, layer->blend[blend].nSubNormal0+slope);
	}
	blend = adjItem(adjNext(adjFirst(layerBlendAdj(layer),normal)));	
	layerBlendNormals(layer, blend, blendnormals);
	if (normal == blendnormals[0] || normal == blendnormals[1] ) {
	  nSubNormal = MIN(nSubNormal, layer->blend[blend].nSubNormal1+slope);
	}else{
	  nSubNormal = MIN(nSubNormal, layer->blend[blend].nSubNormal0+slope);
	}
	if (origSubNormal != nSubNormal) {
	  printf("normal%11d subBlends reduced from%2d to%2d\n",
		 normal, origSubNormal, nSubNormal);
	  allSmooth = FALSE;
	  blend = adjItem(adjFirst(layerBlendAdj(layer),normal));	
	  layerBlendNormals(layer, blend, blendnormals );
	  if (normal == blendnormals[0] || normal == blendnormals[1] ) {
	    layer->blend[blend].nSubNormal0 = nSubNormal;
	  }else{
	    layer->blend[blend].nSubNormal1 = nSubNormal;
	  }
	  blend = adjItem(adjNext(adjFirst(layerBlendAdj(layer),normal)));	
	  layerBlendNormals(layer, blend, blendnormals);
	  if (normal == blendnormals[0] || normal == blendnormals[1] ) {
	    layer->blend[blend].nSubNormal0 = nSubNormal;
	  }else{
	    layer->blend[blend].nSubNormal1 = nSubNormal;
	  }
	}
      }
    }    
  }
  return layer;
}

Layer *layerSubBlendFill(Layer *layer)
{
  int normal;
  AdjIterator it;
  int blend;
  int blendnormals[4];
  double rotation;
  int nSubNormal, subNormal, subNormals[MAXSUBNORMAL];
  int startNormal;
  double axle[3];
  int i;

  if (layerNBlend(layer) <= 0) return layer;

  for (normal=0;normal<layer->originalnormal;normal++){
    switch (layerBlendDegree(layer,normal)) {
    case 0: break;
    case 1: case 2:
      blend = adjItem(adjFirst(layerBlendAdj(layer),normal));
      layerBlendNormals(layer, blend, blendnormals );
      if (normal == blendnormals[0] || normal == blendnormals[1] ) {
	startNormal = blendnormals[0];
	nSubNormal = layer->blend[blend].nSubNormal0;
	for(i=0;i<nSubNormal;i++){
	  subNormals[i] = layerDuplicateNormal(layer, normal );
	  rotation = (double)(i+1) / (double)(nSubNormal+1);
	  layerNormalBlendAxle(layer, normal, axle);
	  gridRotateDirection(layer->normal[blendnormals[0]].direction,
			      layer->normal[blendnormals[1]].direction,
			      axle, rotation,
			      layer->normal[subNormals[i]].direction);
	}
      }else{
	startNormal = blendnormals[2];
	nSubNormal = layer->blend[blend].nSubNormal1;
	for(i=0;i<nSubNormal;i++){
	  subNormals[i] = layerDuplicateNormal(layer, normal );
	  rotation = (double)(i+1) / (double)(nSubNormal+1);
	  layerNormalBlendAxle(layer, normal, axle);
	  gridRotateDirection(layer->normal[blendnormals[2]].direction,
			      layer->normal[blendnormals[3]].direction,
			      axle, rotation,
			      layer->normal[subNormals[i]].direction);
	}
      }
      for ( it = adjFirst(layerBlendAdj(layer),normal); 
	    adjValid(it); 
	    it=adjNext(it) ){
	blend = adjItem(it);
	layerBlendNormals(layer, blend, blendnormals );
	if (normal == blendnormals[0] || normal == blendnormals[1] ) {
	  layer->blend[blend].nSubNormal0 = nSubNormal;
	  if (startNormal == blendnormals[0]) {
	    for(i=0;i<nSubNormal;i++){
	      layer->blend[blend].subNormal0[i] = subNormals[i];
	    }
	  }else{
	    for(i=0;i<nSubNormal;i++){
	      layer->blend[blend].subNormal0[i] = subNormals[nSubNormal-i-1];
	    }
	  }
	}else{
	  layer->blend[blend].nSubNormal1 = nSubNormal;
	  if (startNormal == blendnormals[2]) {
	    for(i=0;i<nSubNormal;i++){
	      layer->blend[blend].subNormal1[i] = subNormals[i];
	    }
	  }else{
	    for(i=0;i<nSubNormal;i++){
	      layer->blend[blend].subNormal1[i] = subNormals[nSubNormal-i-1];
	    }
	  }
	}
      }
      break;
    case 3:
      for ( it = adjFirst(layerBlendAdj(layer),normal); 
	    adjValid(it); 
	    it=adjNext(it) ){
	blend = adjItem(it);
	layerBlendNormals(layer, blend, blendnormals );
	if ( layerNormalRoot(layer, normal) == 
	     layerNormalRoot(layer, blendnormals[0]) ) {
	  nSubNormal = layer->blend[blend].nSubNormal0;
	  for(i=0;i<nSubNormal;i++){
	    subNormal = layerDuplicateNormal(layer, normal );
	    layer->blend[blend].subNormal0[i] = subNormal;
	    rotation = (double)(i+1) / (double)(nSubNormal+1);
	    layerBlendAxle(layer, blend, axle);
	    gridRotateDirection(layer->normal[blendnormals[0]].direction,
				layer->normal[blendnormals[1]].direction,
				axle, rotation,
				layer->normal[subNormal].direction);
	  }
	}else{
	  nSubNormal = layer->blend[blend].nSubNormal1;
	  for(i=0;i<nSubNormal;i++){
	    subNormal = layerDuplicateNormal(layer, normal );
	    layer->blend[blend].subNormal1[i] = subNormal;
	    rotation = (double)(i+1) / (double)(nSubNormal+1);
	    layerBlendAxle(layer, blend, axle);
	    gridRotateDirection(layer->normal[blendnormals[2]].direction,
				layer->normal[blendnormals[3]].direction,
				axle, rotation,
				layer->normal[subNormal].direction);
	  }
	}
      }
      break;
    default:
      printf( "ERROR: %s: %d: Cannot handle %d blends. Write more code!\n",
	      __FILE__, __LINE__, layerBlendDegree(layer,normal));
      break;
    }
    
  }
  return layer;
}

Layer *layerSubBlend(Layer *layer, double maxNormalAngle)
{
  int blend;
  if (layerNBlend(layer) <= 0) return layer;
  if ( 0.0 >= maxNormalAngle ) maxNormalAngle = 30.0;
  if (layer != layerSubBlendCount(layer,maxNormalAngle)) return NULL;
  /* if (layer != layerSubBlendSmooth(layer)) return NULL; */
  for (blend=0; blend < layerNBlend(layer); blend++){
    if (EMPTY != layer->blend[blend].edgeId[0]) 
      layer->blend[blend].nSubNormal0 = layer->blend[blend].nSubNormal1;
    if (EMPTY != layer->blend[blend].edgeId[1]) 
      layer->blend[blend].nSubNormal1 = layer->blend[blend].nSubNormal0;	
  }
  if (layer != layerSubBlendFill(layer)) return NULL;
  return layer;
}   


int layerNSubBlend(Layer *layer, int blend )
{
  if (blend < 0 || blend >= layerNBlend(layer)) return EMPTY;
  return MAX( layer->blend[blend].nSubNormal0,
	      layer->blend[blend].nSubNormal1 ) + 1;
}

Layer *layerExtrudeBlend(Layer *layer, double dx, double dy, double dz )
{

  int blend, i, fix, fixblend, node, newnode, normal, newnormal;
  int normals[4];
  int edgeId, faceId;
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
    if (layer->blend[blend].normal[0]==EMPTY){
      node = layerNormalRoot(layer,layer->blend[blend].oldnormal[0]);
      gridNodeXYZ(grid,node,xyz);
      xyz[0] += dx; xyz[1] += dy; xyz[2] += dz;
      newnode = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);

      normal = layer->blend[blend].oldnormal[0];
      newnormal = layerDuplicateNormal(layer, normal);
      layer->normal[newnormal].root = newnode;
      for ( it = adjFirst(adj,normal); adjValid(it); it=adjNext(it) ){
	fixblend = adjItem(it);
	for (fix=0;fix<4;fix++)
	  if ( layer->blend[fixblend].oldnormal[fix] == normal )
	    layer->blend[fixblend].normal[fix] = newnormal;
      }

      normal = layer->blend[blend].oldnormal[1];
      newnormal = layerDuplicateNormal(layer, normal);
      layer->normal[newnormal].root = newnode;
      for ( it = adjFirst(adj,normal); adjValid(it); it=adjNext(it) ){
	fixblend = adjItem(it);
	for (fix=0;fix<4;fix++)
	  if ( layer->blend[fixblend].oldnormal[fix] == normal )
	    layer->blend[fixblend].normal[fix] = newnormal;
      }
    }

    if (layer->blend[blend].normal[2]==EMPTY){
      node = layerNormalRoot(layer,layer->blend[blend].oldnormal[2]);
      gridNodeXYZ(grid,node,xyz);
      xyz[0] += dx; xyz[1] += dy; xyz[2] += dz;
      newnode = gridAddNode(grid,xyz[0],xyz[1],xyz[2]);

      normal = layer->blend[blend].oldnormal[2];
      newnormal = layerDuplicateNormal(layer, normal);
      layer->normal[newnormal].root = newnode;
      for ( it = adjFirst(adj,normal); adjValid(it); it=adjNext(it) ){
	fixblend = adjItem(it);
	for (fix=0;fix<4;fix++)
	  if ( layer->blend[fixblend].oldnormal[fix] == normal )
	    layer->blend[fixblend].normal[fix] = newnormal;
      }

      normal = layer->blend[blend].oldnormal[3];
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

  for (blend=0; blend < layerNBlend(layer); blend++){
    layer->blend[blend].nodes[0] = 
      layerNormalRoot(layer,layer->blend[blend].normal[0]);
    layer->blend[blend].nodes[1] = 
      layerNormalRoot(layer,layer->blend[blend].normal[2]);
    normals[0] = layer->blend[blend].oldnormal[0];
    normals[1] = layer->blend[blend].oldnormal[2];
    normals[2] = layer->blend[blend].normal[0];
    normals[3] = layer->blend[blend].normal[2];

    layerForceTriangle(layer,normals[0],normals[2],normals[1]);
    layerForceTriangle(layer,normals[1],normals[2],normals[3]);
    normals[0] = layer->blend[blend].oldnormal[3];
    normals[1] = layer->blend[blend].oldnormal[1];
    normals[2] = layer->blend[blend].normal[3];
    normals[3] = layer->blend[blend].normal[1];

    layerForceTriangle(layer,normals[0],normals[2],normals[3]);
    layerForceTriangle(layer,normals[0],normals[3],normals[1]);
  }

  layerBuildNormalTriangleAdjacency(layer);

  for (blend=0; blend < layerNBlend(layer); blend++){
    edgeId = layer->blend[blend].edgeId[0];
    if ( edgeId > 0 ) {
      layerSetParentGeomEdge(layer,
			     layer->blend[blend].normal[0],
			     layer->blend[blend].oldnormal[0],
			     edgeId);
      layerSetParentGeomEdge(layer,
			     layer->blend[blend].normal[1],
			     layer->blend[blend].oldnormal[1],
			     edgeId);
    }
    edgeId = layer->blend[blend].edgeId[1];
    if ( edgeId > 0 ) {
      layerSetParentGeomEdge(layer,
			     layer->blend[blend].normal[2],
			     layer->blend[blend].oldnormal[2],
			     edgeId);
      layerSetParentGeomEdge(layer,
			     layer->blend[blend].normal[3],
			     layer->blend[blend].oldnormal[3],
			     edgeId);
    }
  }

  for (blend=0; blend < layerNBlend(layer); blend++){
    faceId = layerConstrained(layer, layer->blend[blend].normal[0]);
    if (faceId > 0 ) {
      layerConstrainTriangleSide(layer,
				 layer->blend[blend].normal[0],
				 layer->blend[blend].oldnormal[0],
				 faceId);
      layerConstrainTriangleSide(layer,
				 layer->blend[blend].normal[1],
				 layer->blend[blend].oldnormal[1],
				 faceId);
    }
    faceId = layerConstrained(layer, layer->blend[blend].normal[2]);
    if (faceId > 0 ) {
      layerConstrainTriangleSide(layer,
				 layer->blend[blend].normal[2],
				 layer->blend[blend].oldnormal[2],
				 faceId);
      layerConstrainTriangleSide(layer,
				 layer->blend[blend].normal[3],
				 layer->blend[blend].oldnormal[3],
				 faceId);
    }
  }

  return layer;
}

Layer *layerOrderedVertexBlends(Layer *layer, int normal, 
				 int *nVertexBlends, int *vertexBlends ){

  int node;
  AdjIterator it;
  int blend, nextNormal;
  int count;

  if (0 == layerBlendDegree(layer, normal)) return layer;

  node = layerNormalRoot(layer,normal);

  blend = adjItem(adjFirst(layerBlendAdj(layer),normal));
  vertexBlends[0] = blend;
  if (node == layer->blend[blend].nodes[0]) {
    nextNormal = layer->blend[blend].normal[0];
  }else{
    nextNormal = layer->blend[blend].normal[3];
  }

  for (count=1;count<*nVertexBlends;count++){
    for ( it = adjFirst(layerBlendAdj(layer),normal);
	  adjValid(it);
	  it=adjNext(it) ){
      blend = adjItem(it);
      if (node == layer->blend[blend].nodes[0]) {
	if (nextNormal == layer->blend[blend].normal[1]){
	  vertexBlends[count] = blend;
	  nextNormal = layer->blend[blend].normal[0];
	  break;
	}
      }else{
	if (nextNormal == layer->blend[blend].normal[2]){
	  vertexBlends[count] = blend;
	  nextNormal = layer->blend[blend].normal[3];
	  break;
	}
      }
    }
  }
  
  return layer;
}

Layer *layerOrderedVertexNormals(Layer *layer, int normal, 
				 int *nVertexNormals, int *vertexNormals ){
  int node;
  int i, blend, subNormal, nSubNormal;
  int count;
  int nVertexBlends, *vertexBlends;

  nVertexBlends = layerBlendDegree(layer, normal);

  if (0 == nVertexBlends) return layer;
  vertexBlends = malloc(nVertexBlends*sizeof(int));
  layerOrderedVertexBlends(layer, normal, &nVertexBlends, vertexBlends );
  
  node = layerNormalRoot(layer,normal);
  count = 0;
  for (i=0;i<nVertexBlends;i++){
    blend = vertexBlends[i];
    if (node == layer->blend[blend].nodes[0]) {
      vertexNormals[count] = layer->blend[blend].normal[1];
      count++;
      nSubNormal = layer->blend[blend].nSubNormal0;
      for (subNormal=0;subNormal<nSubNormal;subNormal++){
	vertexNormals[count] = 
	  layer->blend[blend].subNormal0[nSubNormal-subNormal-1];
	count++;
      }
    }else{
      vertexNormals[count] = layer->blend[blend].normal[2];
      count++;
      nSubNormal = layer->blend[blend].nSubNormal1;
      for (subNormal=0;subNormal<nSubNormal;subNormal++){
	vertexNormals[count] = 
	  layer->blend[blend].subNormal1[subNormal];
	count++;
      }
    }
  }
  
  free(vertexBlends);

  return layer;
}


Layer *layerPopulateNormalNearTree(Layer *layer)
{
  int normal;
  Grid *grid;
  double xyz[3], height, edgeLength, radius;

  grid = layerGrid(layer);
  if (NULL != layer->nearTree) free(layer->nearTree);

  layer->nearTree = malloc(layerNNormal(layer)*sizeof(Near));

  for(normal=0;normal<layerNNormal(layer);normal++){
    gridNodeXYZ(grid, layerNormalRoot(layer, normal ), xyz);
    layerGetNormalHeight(layer, normal, &height);
    layerNormalMaxEdgeLength(layer, normal, &edgeLength);
    radius = MAX(edgeLength, 2*height);
    nearInit(&layer->nearTree[normal], normal, xyz[0], xyz[1], xyz[2], radius);
    if (normal>0) nearInsert(layer->nearTree,&layer->nearTree[normal]);
  }

  return layer;
}

Layer *layerPopulateTriangleNearTree(Layer *layer)
{
  int i, triangle;
  Grid *grid;
  double node0[3], node1[3], node2[3], node3[3];
  double center[3], dist[3], radius;

  grid = layerGrid(layer);
  if (NULL != layer->nearTree) free(layer->nearTree);

  layer->nearTree = malloc(layerNTriangle(layer)*sizeof(Near));

  for(triangle=0;triangle<layerNTriangle(layer);triangle++){
    if (layer != layerTriangleInviscidTet(layer,triangle,
					  node0,node1,node2,node3)) return NULL;
    for (i=0;i<3;i++) center[i] = 0.25*(node0[i]+node1[i]+node2[i]+node3[i]);
    gridSubtractVector(node0,center,dist);
    radius = gridVectorLength(dist);
    gridSubtractVector(node1,center,dist);
    radius = MAX(radius,gridVectorLength(dist));
    gridSubtractVector(node2,center,dist);
    radius = MAX(radius,gridVectorLength(dist));
    gridSubtractVector(node3,center,dist);
    radius = MAX(radius,gridVectorLength(dist));
    nearInit(&layer->nearTree[triangle], 
	     triangle, center[0], center[1], center[2], radius);
    if (triangle>0) nearInsert(layer->nearTree,&layer->nearTree[triangle]);    
  }

  return layer;
}

Layer *layerTerminateCollidingNormals(Layer *layer)
{
  int normal;
  Near *target;
  int i, touched, maxTouched, *nearNormals;
  double dir1[3], dir2[3], dot;
  double xyz1[3], xyz2[3], view[3];
  double length, visible, reciprocal;

  if ( 0 < layerNBlend(layer) ) return NULL;

  //printf("layerPopulateNormalNearTree...\n");
  layerPopulateNormalNearTree(layer);

  maxTouched = layerNNormal(layer);
  nearNormals = malloc(maxTouched*sizeof(int));

  //printf("inspecting normal proximity...\n");
  for(normal=0;normal<layerNNormal(layer);normal++){
    target = &layer->nearTree[normal];
    touched = 0;
    nearTouched(layer->nearTree, target, &touched, maxTouched, nearNormals);
    layerNormalDirection(layer, normal, dir1);
    for(i=0;i<touched;i++){
      layerNormalDirection(layer, nearNormals[i], dir2);
      dot = gridDotProduct( dir1, dir2 );
      if ( dot < -0.2 ) {
	gridNodeXYZ(layerGrid(layer), layerNormalRoot(layer, normal ), xyz1);
	gridNodeXYZ(layerGrid(layer), layerNormalRoot(layer, nearNormals[i] ), 
		    xyz2);
	gridSubtractVector(xyz2, xyz1, view);
	length = sqrt ( gridDotProduct( view, view ) );
	view[0]= view[0]/length;
	view[1]= view[1]/length;
	view[2]= view[2]/length;
	visible = gridDotProduct( dir1, view );
	gridSubtractVector(xyz1, xyz2, view);
	view[0]= view[0]/length;
	view[1]= view[1]/length;
	view[2]= view[2]/length;
	reciprocal = gridDotProduct( dir2, view );
	if (visible > 0.2 && reciprocal > 0.2) {
	  layerTerminateNormal(layer,normal);
	  layerTerminateNormal(layer,nearNormals[i]);
	}
      }
    }
  }
  free(nearNormals);

  return layer;
}

Layer *layerTerminateFutureNegativeCellNormals(Layer *layer)
{
  Grid *grid = layer->grid;
  int triangle, normals[3], nodes[3], n[6], tet[4];
  int i, norm, normal;
  int  root, tip, node, nterminated;
  double xyz[3];
  int savedNodes[3];
  double volumeLimit = 1.0e-14;

  if (layerNNormal(layer) == 0 ) return NULL;

  for (norm=0;norm<3;norm++){
    savedNodes[norm] = gridAddNode(grid,0.0,0.0,0.0);
    if (EMPTY == savedNodes[norm]) return NULL;
  }
  
  for (triangle=0;triangle<layerNTriangle(layer);triangle++){
    layerTriangleNormals(layer, triangle, normals);
    layerTriangle(layer, triangle, nodes);
    
    for (norm=0;norm<3;norm++){
      normal = normals[norm];
      if (layerNormalTerminated(layer,normal)){
	layer->normal[normal].tip = layer->normal[normal].root;
      }else{
	root = layer->normal[normal].root;
	gridNodeXYZ(grid,root,xyz);
	for(i=0;i<3;i++)
	  xyz[i] = xyz[i] 
	    + layer->normal[normal].height
	    * layer->normal[normal].direction[i];
	tip = savedNodes[norm];
	gridSetNodeXYZ(grid,tip,xyz);
	layer->normal[normal].tip = tip;
      }
    }

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
    
    if (nodes[2]<nodes[1]){
      if (n[0]!=n[3]) {
	tet[0] = n[0]; tet[1] = n[4]; tet[2] = n[5]; tet[3] = n[3];
	if ( volumeLimit > gridVolume(grid, tet ) ) 
	  layerTerminateNormal( layer, normals[0] );
      }
      if (n[2]!=n[5]) {
	tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[4]; tet[3] = n[5];
	if ( volumeLimit > gridVolume(grid, tet ) ) 
	  layerTerminateNormal( layer, normals[1] );
      }
      if (n[1]!=n[4]) {
	tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[1]; tet[3] = n[4];
	if ( volumeLimit > gridVolume(grid, tet ) ) 
	  layerTerminateNormal( layer, normals[2] );
      }
    }else{
      if (n[0]!=n[3]) {
	tet[0] = n[0]; tet[1] = n[4]; tet[2] = n[5]; tet[3] = n[3];
	if ( volumeLimit > gridVolume(grid, tet ) ) 
	  layerTerminateNormal( layer, normals[0] );
      }
      if (n[1]!=n[4]) {
	tet[0] = n[0]; tet[1] = n[1]; tet[2] = n[5]; tet[3] = n[4];
	if ( volumeLimit > gridVolume(grid, tet ) ) 
	  layerTerminateNormal( layer, normals[1] );
      }
      if (n[2]!=n[5]) {
	tet[0] = n[2]; tet[1] = n[0]; tet[2] = n[1]; tet[3] = n[5];
	if ( volumeLimit > gridVolume(grid, tet ) ) 
	  layerTerminateNormal( layer, normals[2] );
      }
    }
    for (norm=0;norm<3;norm++){
      layer->normal[normals[norm]].tip = EMPTY;
    }
  }
  
  for (norm=0;norm<3;norm++) gridRemoveNode(grid,savedNodes[norm]);

  return layer;
}

GridBool layerTrianglesShareNormal(Layer *layer, int triangle1, int triangle2 )
{
  int normals1[3], normals2[3];
  int i, j;
  if ( NULL == layerTriangleNormals(layer, triangle1, normals1) ) return FALSE;
  if ( NULL == layerTriangleNormals(layer, triangle2, normals2) ) return FALSE;
  
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      if (normals1[i] == normals2[j]) return TRUE;

  return FALSE;
}


Layer *layerTerminateCollidingTriangles(Layer *layer)
{
  int triangle;
  Near *target;
  int i, touched, maxTouched, *nearTriangles;
  int nearTriangle;
  double a0[3], a1[3], a2[3], a3[3];
  double b0[3], b1[3], b2[3], b3[3];

  if ( 0 < layerNBlend(layer) ) return NULL;

  //printf("layerPopulateTriangleNearTree...\n");
  layerPopulateTriangleNearTree(layer);

  maxTouched = layerNTriangle(layer);
  nearTriangles = malloc(maxTouched*sizeof(int));

  //printf("inspecting triangle proximity...\n");
  for(triangle=0;triangle<layerNTriangle(layer);triangle++){
    target = &layer->nearTree[triangle];
    touched = 0;
    nearTouched(layer->nearTree, target, &touched, maxTouched, nearTriangles);
    //printf("triangle %d  touched %d\n",triangle,touched);
    layerTriangleInviscidTet(layer,triangle,a0,a1,a2,a3);
    for(i=0;i<touched;i++){
      nearTriangle = nearTriangles[i];
      if (!layerTrianglesShareNormal(layer, nearTriangle, triangle)) {
	layerTriangleInviscidTet(layer,nearTriangle,b0,b1,b2,b3);
	if (intersectTetTet(a0,a1,a2,a3,b0,b1,b2,b3)){
	  layerTerminateTriangleNormals(layer,triangle);
	  layerTerminateTriangleNormals(layer,nearTriangle);
	}
      }
    }
  }
  free(nearTriangles);

  return layer;
}

Layer *layerSmoothLayerWithHeight(Layer *layer)
{
  int normal, node;
  AdjIterator it;
  int triangle, nodes[3];
  double direction[3], xyz0[3], xyz1[3], edge[3];
  double avg, correction;
  int i, hits;
  Grid *grid;

  grid = layerGrid(layer);

  for ( normal = 0 ; normal < layerNNormal(layer) ; normal++ ) {
    node = layerNormalRoot(layer,normal);
    gridNodeXYZ(grid,node,xyz0);
    layerNormalDirection(layer,normal,direction);
    hits = 0;
    avg = 0;
    for ( it = adjFirst(layer->triangleAdj,normal);
	  adjValid(it);
	  it = adjNext(it) ){
      triangle = adjItem(it);
      layerTriangle(layer,triangle,nodes);
      for (i=0;i<3;i++){
	if (nodes[i]!=node) {
	  gridNodeXYZ(grid,nodes[i],xyz1);
	  gridSubtractVector(xyz1,xyz0,edge);
	  gridVectorNormalize(edge);
	  avg += gridDotProduct(edge,direction);
	  hits++;
	}
      }
    }
    avg = avg / (double)hits;
    //printf("%10.6f\n",avg);
    if (avg <= 0) {
      correction = 1 + avg;
    }else{
      correction = 1 + avg;
      correction = correction * correction;
      correction = MIN(correction,2.0);
    }
    layer->normal[normal].height =
      layer->normal[normal].height * correction;
  }

  return layer;
}

Layer *layerOffsetTriangleMR(Layer *layer, 
			     int triangle, int normal, double offset,
			     double *MR, double *dMRdX)
{
  double xyz1[3], xyz2[3], xyz3[3];
  int i, normals[3], temp;
  double direction[3];
  Grid *grid;

  grid = layerGrid(layer);

  if ( layer != layerTriangleNormals(layer, triangle, normals) ) return NULL;

  if (normals[0] != normal) {
    temp = normals[0];
    normals[0] = normals[1];
    normals[1] = normals[2];
    normals[2] = temp;
  }
  if (normals[0] != normal) {
    temp = normals[0];
    normals[0] = normals[1];
    normals[1] = normals[2];
    normals[2] = temp;
  }
  if (normals[0] != normal) return NULL;

  gridNodeXYZ(grid,layerNormalRoot(layer,normals[0]),xyz1);
  layerNormalDirection(layer,normals[0],direction);
  for(i=0;i<3;i++) xyz1[i] += offset* direction[i];

  gridNodeXYZ(grid,layerNormalRoot(layer,normals[1]),xyz2);
  layerNormalDirection(layer,normals[1],direction);
  for(i=0;i<3;i++) xyz2[i] += offset* direction[i];

  gridNodeXYZ(grid,layerNormalRoot(layer,normals[2]),xyz3);
  layerNormalDirection(layer,normals[2],direction);
  for(i=0;i<3;i++) xyz3[i] += offset* direction[i];

  FaceMRDerivative( xyz1[0], xyz1[1], xyz1[2],
		    xyz2[0], xyz2[1], xyz2[2],
		    xyz3[0], xyz3[1], xyz3[2],
		    MR, dMRdX );

  return layer;

}

Layer *layerOptimizeNormalDirection(Layer *layer, double offsetRatio)
{
  int i, normal, triangle;
  double worstMR, worstdMRdX[3], MR, dMRdX[3];
  double projection, relax, offset;
  AdjIterator it;
  
  relax  = 0.01;

  for (normal=0;normal<layerNNormal(layer);normal++){
    if ( 0 != layerConstrained(layer,normal) ){
      layerNormalMaxEdgeLength(layer,normal,&offset);
      offset *= offsetRatio;
      worstMR = 2.0;
      for ( it = adjFirst(layer->triangleAdj,normal); 
	    adjValid(it); 
	    it = adjNext(it) ){
	triangle = adjItem(it);
	layerOffsetTriangleMR(layer, triangle, normal, offset, &MR, dMRdX );
	if (MR<worstMR) {
	  worstMR = MR;
	  for(i=0;i<3;i++) worstdMRdX[i] = dMRdX[i];
	}
      }
      projection = gridDotProduct(layer->normal[normal].direction, worstdMRdX);
      for(i=0;i<3;i++) worstdMRdX[i] -= projection*layer->normal[normal].direction[i];
      gridVectorNormalize(worstdMRdX);
      for(i=0;i<3;i++) layer->normal[normal].direction[i] += relax*worstdMRdX[i];

      gridVectorNormalize(layer->normal[normal].direction);
    }
  }
  
  layerVisibleNormals(layer,0.25,1.0e-5);
  
  for (normal=0;normal<layerNNormal(layer);normal++){
    if ( 0 == layerConstrained(layer,normal) ){
      layerNormalMaxEdgeLength(layer,normal,&offset);
      offset *= offsetRatio;
      worstMR = 2.0;
      for ( it = adjFirst(layer->triangleAdj,normal); 
	    adjValid(it); 
	    it = adjNext(it) ){
	triangle = adjItem(it);
	layerOffsetTriangleMR(layer, triangle, normal, offset, &MR, dMRdX );
	if (MR<worstMR) {
	  worstMR = MR;
	  for(i=0;i<3;i++) worstdMRdX[i] = dMRdX[i];
	}
      }
      projection = gridDotProduct(layer->normal[normal].direction, worstdMRdX);
      for(i=0;i<3;i++) worstdMRdX[i] -= projection*layer->normal[normal].direction[i];
      gridVectorNormalize(worstdMRdX);
      for(i=0;i<3;i++) layer->normal[normal].direction[i] += relax*worstdMRdX[i];

      gridVectorNormalize(layer->normal[normal].direction);
    }
  }
  
  layerVisibleNormals(layer,0.25,1.0e-5);
  
  return layer;
}

Layer *layerWriteTecplotFrontWithData(Layer *layer, int nn )
{
  int i;
  double xyz[3];
  Grid *grid;

  double height, L, rate, Lmax, r;

  grid = layerGrid(layer);

  if ( NULL == layer->tecplotFile) {
    layer->tecplotFile = fopen("layer.plt","w");
    fprintf(layer->tecplotFile, "title=\"tecplot advancing layer\"\n");
    fprintf(layer->tecplotFile, "variables=\"X\",\"Y\",\"Z\",\"r\", \"id\", \"h\",\"L\",\"rate\",\"L_m_a_x\" \"term\" \"n\" \n");
  }

  fprintf(layer->tecplotFile, 
	  "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	  layerNNormal(layer), layerNTriangle(layer));

  for ( i=0; i<layerNNormal(layer) ; i++ ){
    int term=0;
    if( layerNormalTerminated(layer, i) ) term = 1;
    gridNodeXYZ(grid,layerNormalRoot(layer,i),xyz);
    r      = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
    height = layer->normal[i].height;
    L      = layer->normal[i].length;
    rate   = layer->normal[i].rate;
    Lmax   = layer->normal[i].maxlength;
    fprintf(layer->tecplotFile, "%23.15e %23.15e %23.15e %23.15e %d %23.15e %23.15e %23.15e %23.15e %d %d\n",xyz[0],xyz[1],xyz[2],r,i,height,L,rate,Lmax,term,nn);
  }

  fprintf(layer->tecplotFile, "\n");

  for ( i=0; i<layerNTriangle(layer) ; i++ ){
    fprintf(layer->tecplotFile, " %9d %9d %9d\n",
	    layer->triangle[i].normal[0]+1,
	    layer->triangle[i].normal[1]+1,
	    layer->triangle[i].normal[2]+1);
  }

  fflush(layer->tecplotFile);

  return layer;
}

Layer *layerWriteTecplotFrontGeometry(Layer *layer)
{
  int i;
  double xyz[3];
  Grid *grid;

  grid = layerGrid(layer);

  if ( NULL == layer->tecplotFile) {
    layer->tecplotFile = fopen("layer.t","w");
    fprintf(layer->tecplotFile, "title=\"tecplot advancing layer\"\n");
    fprintf(layer->tecplotFile, "variables=\"X\",\"Y\",\"Z\"\n");
  }

  fprintf(layer->tecplotFile, 
	  "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	  layerNNormal(layer), layerNTriangle(layer));

  for ( i=0; i<layerNNormal(layer) ; i++ ){
    gridNodeXYZ(grid,layerNormalRoot(layer,i),xyz);
    fprintf(layer->tecplotFile, "%23.15e%23.15e%23.15e\n",xyz[0],xyz[1],xyz[2]);
  }

  fprintf(layer->tecplotFile, "\n");

  for ( i=0; i<layerNTriangle(layer) ; i++ ){
    fprintf(layer->tecplotFile, " %9d %9d %9d\n",
	    layer->triangle[i].normal[0]+1,
	    layer->triangle[i].normal[1]+1,
	    layer->triangle[i].normal[2]+1);
  }

  fflush(layer->tecplotFile);

  return layer;
}

