
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "gridfiller.h"
#include "gridgeom.h"
#ifdef HAVE_SDK
#include "Goolache/CAPrI_IO.h"
#include "CADGeom/CADGeom.h"
#else
#include "FAKEGeom.h"
#endif

Grid *gridLoadPart( char *modeler, char *project, int maxnode )
{

  printf("calling CADGeom_Start ... \n");
  if ( ! CADGeom_Start( ) ){
    printf("ERROR: CADGeom_Start broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  

  printf("calling CADGeom_Load for project <%s> with %s kernel... \n",project,modeler);
#ifdef HAVE_CAPRI2
  if ( ! GeoMesh_LoadPart( modeler, project ) ){
#else
  if ( ! GeoMesh_LoadPart( project ) ){
#endif
    printf("ERROR: GeoMesh_LoadPart broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  return gridFillFromPart( 1, maxnode );

}
 
Grid *gridFillFromPart( int vol, int maxnode )
{

  Grid *grid;
  UGridPtr ugrid;
  CADCurvePtr edge;
  UGPatchPtr  localPatch, globalPatch;
  Iterator patchIterator;
  int gridDimensions[3];
  int patchDimensions[3];
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;
  int nedgenode;
  int nnode, nface, ncell, nedge;
  int maxface, maxcell, maxedge;
  int i, iedge, inode;
  int face, localNode, globalNode;
  double *xyz;
  int *c2n, *f2n, *faceId;
  double trange[2];
  int edgeEndPoint[2];

  if (NULL == (ugrid = CADGeom_VolumeGrid(vol)) ) {
    printf("ERROR: Can not find grid in restart. \n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  if( !UGrid_GetDims(ugrid,gridDimensions) ) {
    printf("ERROR: Could not get grid dimensions. \n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

  nnode = gridDimensions[0];
  nface = gridDimensions[1];
  ncell = gridDimensions[2];

  if( !CADGeom_GetVolume(vol,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) ) {
    printf("ERROR: CADGeom_GetVolume. \n%s\n",ErrMgr_GetErrStr());
  }

  nedge =0 ;
  for( iedge=1; iedge<=nGeomEdge; iedge++ ) {
    if( (edge=CADGeom_EdgeGrid(vol,iedge)) == NULL ) 
      printf("ERROR: CADGeom_EdgeGrid(%d).\n%s\n",iedge,ErrMgr_GetErrStr());
    nedge += (CADCURVE_NUMPTS(edge)-1);
  }

  printf("ugrid size: %d nodes %d faces %d cells %d edge elements.\n",
	 nnode,nface,ncell,nedge);

  maxface = maxnode;
  maxcell = maxnode * 6;
  maxedge = maxnode / 10;

  printf("max grid size: %d nodes %d faces %d cells %d edge elements.\n",
	 maxnode,maxface,maxcell,maxedge);

  if (  nnode>maxnode || nface>maxface || ncell>maxcell || nedge>maxedge){
    printf("ERROR: gridFillFromPart: max grid size is smaller then ugrid\n");
    return NULL;
  }

  c2n    = malloc( maxcell * 4 * sizeof(int) );
  f2n    = malloc( maxface * 3 * sizeof(int) );
  faceId = malloc( maxface *     sizeof(int) );
  xyz    = malloc( maxnode * 3 * sizeof(double) );

  for( i=0; i<nnode; i++ ) {
    xyz[0+3*i] = UGrid_PtValue(ugrid,i,0);
    xyz[1+3*i] = UGrid_PtValue(ugrid,i,1);
    xyz[2+3*i] = UGrid_PtValue(ugrid,i,2);
  }

  for( i=0; i<nface; i++ ) {
    f2n[0+3*i] = UGrid_VertValue(ugrid,i,0);
    f2n[1+3*i] = UGrid_VertValue(ugrid,i,1);
    f2n[2+3*i] = UGrid_VertValue(ugrid,i,2);
    faceId[i]  = UGrid_FlagValue(ugrid,i);
  }

  for( i=0; i<ncell; i++ ) {
    c2n[0+4*i] = UGrid_TetValue(ugrid,i,0);
    c2n[1+4*i] = UGrid_TetValue(ugrid,i,1);
    c2n[2+4*i] = UGrid_TetValue(ugrid,i,2);
    c2n[3+4*i] = UGrid_TetValue(ugrid,i,3);
  }

  grid = gridImport( maxnode, nnode, maxface, nface, 
		     maxcell, ncell, maxedge,
		     xyz, f2n, faceId, c2n );

  gridSetNGeomNode( grid, nGeomNode );
  gridSetNGeomEdge( grid, nGeomEdge );
  gridSetNGeomFace( grid, nGeomFace );

  inode = nGeomNode;

  for( iedge=1; iedge<=nGeomEdge; iedge++ ) {
    if( (edge=CADGeom_EdgeGrid(vol,iedge)) == NULL ) 
      printf("ERROR: CADGeom_EdgeGrid(%d).\n%s\n",iedge,ErrMgr_GetErrStr());
 
    nedgenode = CADCURVE_NUMPTS(edge);

    CADGeom_GetEdge( vol, iedge, trange, edgeEndPoint );

    edgeEndPoint[0]--; /* convert from fortran to c numbers */
    edgeEndPoint[1]--;

    gridAddGeomEdge( grid, iedge, edgeEndPoint[0], edgeEndPoint[1]);

    if (nedgenode == 2) {
      gridAddEdge(grid, edgeEndPoint[0], edgeEndPoint[1], 
		  iedge, trange[0], trange[1]);
    }else{
      gridAddEdge(grid, edgeEndPoint[0], inode, iedge,
		  edge->param[0], edge->param[1]);
      for( i=1 ; i < (nedgenode-2) ; i++ ) { /* skip end segments */
	gridAddEdge(grid, inode, inode+1, iedge,
		  edge->param[i], edge->param[i+1]);
	inode++;
      }
      gridAddEdge(grid, inode, edgeEndPoint[1], iedge,
		  edge->param[nedgenode-2], 
		  edge->param[nedgenode-1]);
      inode++;
    }
  }

  if ( nedge != gridNEdge(grid) )
    printf("ERROR: gridLoadPart: %s: %d: nedge != gridNEdge(grid)\n",
	   __FILE__,__LINE__);

  /* get uv vals for surface(s) */
  /* we use globalPatch to track with the localPatch so that we can get global
   * node numbering relative the volume grid and NOT the face grid as would
   * be the case of global index of upp
   */

  globalPatch = DList_SetIteratorToHead(UGrid_PatchList(ugrid),&patchIterator);

  for( face=1; face<=nGeomFace; face++ ) {
    char filename[265];

    localPatch = CADGeom_FaceGrid(vol,face);
    sprintf(filename,"localUGPatch%03d.t",face);
    sprintf(filename,"globalUGPatch%03d.t",face);
    UGPatch_GetDims(localPatch,patchDimensions);
    for( localNode=0; localNode<patchDimensions[0]; localNode++ ) {	     
      globalNode = UGPatch_GlobalIndex(globalPatch,localNode);
      gridSetNodeUV( grid, globalNode, face,
		     UGPatch_Parameter(localPatch,localNode,0), 
		     UGPatch_Parameter(localPatch,localNode,1));
    }

    globalPatch = DList_GetNextItem(&patchIterator);
  }

  return grid;
}

int gridSavePart( Grid *grid, char *project )
{
  int vol=1;
  UGridPtr ugrid;  
  int nnode, nface, ncell;
  double *xyz;
  int *f2n, *faceId, *c2n;
  int i, ixyz, iface;
  int iedge, curveEndPoint[2], nCurveNode, *curve;
  double trange[2];
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;

  double *temp_xyz, *temp_tuv;

  UGPatchPtr  localPatch, globalPatch;
  Iterator patchIterator;
  int patchDimensions[3];
  int face, localNode, globalNode;



  if( !CADGeom_GetVolume(vol,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) )
    printf("ERROR: CADGeom_GetVolume, line %d of %s\n.",__LINE__, __FILE__);

  gridSortNodeGridEx( grid );

  gridExport( grid, &nnode, &nface, &ncell,
	      &xyz, &f2n, &faceId, &c2n );
  
  /* edge stuff */
  for (iedge=1; iedge<=nGeomEdge; iedge++){

    CADGeom_GetEdge( vol, iedge, trange, curveEndPoint);
    curveEndPoint[0]--; curveEndPoint[1]--;/* fortran to c numbering */
    
    nCurveNode = gridGeomCurveSize( grid, iedge, curveEndPoint[0]);
    curve =    malloc( nCurveNode *     sizeof(int) );
    temp_xyz = malloc( nCurveNode * 3 * sizeof(double) );
    temp_tuv = malloc( nCurveNode *     sizeof(double) );

    gridGeomCurve( grid, iedge, curveEndPoint[0], curve );
    gridGeomCurveT( grid, iedge, curveEndPoint[0], temp_tuv );

    for ( i=0; i<nCurveNode; i++) /* include end points */
      for ( ixyz=0; ixyz<3 ; ixyz++)
	temp_xyz[ixyz+3*i] = xyz[ixyz+3*curve[i]];

    CADGeom_UpdateEdgeGrid( vol, iedge, nCurveNode, temp_xyz, temp_tuv );

    free(curve);
  }

  if ( !UGrid_FromArrays( &ugrid, nnode, xyz, nface, f2n, ncell, c2n  )) {
    printf(" Could not make UGridPtr, line %d of %s\n", __LINE__, __FILE__);
    return(-1);
  }
  
  for (iface = 0 ; iface < nface ; iface++ ) {
    UGrid_FlagValue(ugrid,iface) = faceId[iface];
  }

  printf("Rebuilding Element Connectivity...");
  UGrid_BuildConnectivity(ugrid);		/* Build Connectivity */
  printf("Complete\n");

  if( !UGPatch_InitSurfacePatches(ugrid) ) {
    printf(" Could not make surface patches for new UGridPtr, line %d of %s\n",
	   __LINE__, __FILE__);    
    return(-1);
  }

  UGrid_TIMESTAMP(ugrid) = time( NULL );	/* Updated time */
  UGrid_ALGORITHM(ugrid) = UGrid_ALGORITHM(CADGeom_VolumeGrid(vol));

  if( !CADGeom_SetVolumeGrid( vol, ugrid ) ) {
    printf(" Could not replace CADGeom volume grid, line %d of %s\n",
	   __LINE__, __FILE__);    
    return(-1);
  }

  if( !CADTopo_FlushPatches(vol,ugrid) ) {
    printf(" Could not flush patches CADTopo, line %d of %s\n",
	   __LINE__, __FILE__);    
    return(-1);
  }

  /* get uv vals for surface(s) */
  /* we use globalPatch to track with the localPatch so that we can get global
   * node numbering relative the volume grid and NOT the face grid as would
   * be the case of global index of upp
   */

  globalPatch = DList_SetIteratorToHead(UGrid_PatchList(ugrid),&patchIterator);

  for( face=1; face<=nGeomFace; face++ ) {
    double uv[2];
    char filename[265];

    localPatch = CADGeom_FaceGrid(vol,face);
    sprintf(filename,"localUGPatch%03d.t",face);
    sprintf(filename,"globalUGPatch%03d.t",face);
    UGPatch_GetDims(localPatch,patchDimensions);
    for( localNode=0; localNode<patchDimensions[0]; localNode++ ) {	     
      globalNode = UGPatch_GlobalIndex(globalPatch,localNode);
      gridNodeUV( grid, globalNode, face, uv );
      UGPatch_Parameter(localPatch,localNode,0) = uv[0];
      UGPatch_Parameter(localPatch,localNode,1) = uv[1];
   }

    globalPatch = DList_GetNextItem(&patchIterator);
  }

  if ( NULL != project ) {
    GeoMesh_UseDefaultIOCallbacks();
    if( !CADGeom_SavePart(vol,project) ) {
      printf("Yo! Could NOT save \"%s\".\n",project);
    }
  }

  return 0;
}

 
int gridSavePartExplicitly( Grid *grid, char *project )
{
  int vol=1;
  UGridPtr ugrid;  
  int nnode, nface, ncell;
  double *xyz;
  int *f2n, *faceId, *c2n;
  int i, ixyz, iface;
  int iedge, curveEndPoint[2], nCurveNode, *curve;
  double trange[2];
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;

  double *temp_xyz, *temp_tuv;

  UGPatchPtr  localPatch, globalPatch;
  Iterator patchIterator;
  int patchDimensions[3];
  int face, localNode, globalNode;

  int local_nnode, local_nface;
  int *local_f2n;
  double *local_xyz, *local_uv;
  int *faceg2l, *facel2g;
  int face_edge_count;
  int node;
  int min_face_interior_node, max_face_interior_node;

  if( !CADGeom_GetVolume(vol,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) )
    printf("ERROR: CADGeom_GetVolume, line %d of %s\n.",__LINE__, __FILE__);

  gridSortNodeGridEx( grid );

  gridExport( grid, &nnode, &nface, &ncell,
	      &xyz, &f2n, &faceId, &c2n );
  
  /* edge stuff */
  for (iedge=1; iedge<=nGeomEdge; iedge++){

    CADGeom_GetEdge( vol, iedge, trange, curveEndPoint);
    curveEndPoint[0]--; curveEndPoint[1]--;/* fortran to c numbering */
    
    nCurveNode = gridGeomCurveSize( grid, iedge, curveEndPoint[0]);
    curve =    malloc( nCurveNode *     sizeof(int) );
    temp_xyz = malloc( nCurveNode * 3 * sizeof(double) );
    temp_tuv = malloc( nCurveNode *     sizeof(double) );

    gridGeomCurve( grid, iedge, curveEndPoint[0], curve );
    gridGeomCurveT( grid, iedge, curveEndPoint[0], temp_tuv );

    for ( i=0; i<nCurveNode; i++) /* include end points */
      for ( ixyz=0; ixyz<3 ; ixyz++)
	temp_xyz[ixyz+3*i] = xyz[ixyz+3*curve[i]];

    CADGeom_UpdateEdgeGrid( vol, iedge, nCurveNode, temp_xyz, temp_tuv );

    free(curve);
  }

  /* face stuff */
  faceg2l = malloc( nnode * sizeof(int) );
  for( face=1; face<=nGeomFace; face++ ) {
    face_edge_count = gridFaceEdgeCount( face );
    facel2g = malloc( face_edge_count * sizeof(int) );
    for (node=0;node<face_edge_count;node++) facel2g[node] = EMPTY;
    gridFaceEdgeLocal2Global( grid, face, face_edge_count, facel2g );
    for (node=0;node<nnode;node++) faceg2l[node] = EMPTY;
    for (node=0;node<face_edge_count;node++) faceg2l[facel2g[node]] = node;

    min_face_interior_node = nnode;
    max_face_interior_node = -1;
    local_nface = 0;
    for (iface=0;iface<nface;iface++){
      if ( face == faceId[iface]) {
	local_nface++;
	for (node=0;node<3;node++){
	  if ( faceg2l[f2n[node+3*iface]] < 0 ) {
	    min_face_interior_node=MIN( min_face_interior_node,
					f2n[node+3*iface] );
	    max_face_interior_node=MAX( max_face_interior_node,
					f2n[node+3*iface] );
	  }
	}
      }
    }
    local_nnode = face_edge_count;
    for ( node  = min_face_interior_node ;
	  node <= max_face_interior_node ; node++ ) {
      faceg2l[node] = local_nnode;
      local_nnode++;
    }

    local_xyz = malloc( local_nnode * 3 * sizeof(double) );
    local_uv  = malloc( local_nnode * 2 * sizeof(double) );
    local_f2n = malloc( local_nface * 3 * sizeof(int) );

    for ( node = 0 ; node < nnode ; node++ ) {
      if ( EMPTY != faceg2l[node] ) {
	gridNodeXYZ( grid, node, &(local_xyz[3*faceg2l[node]]) );
	gridNodeUV( grid, node, face, &(local_uv[2*faceg2l[node]]) );
      }
    }
    local_nface = 0;
    for (iface=0;iface<nface;iface++){
      if ( face == faceId[iface]) {
	for (node=0;node<3;node++){
	  local_f2n[node+3*local_nface] = f2n[node+3*iface];
	}
	local_nface++;
      }
    }  

    gridUpdateGeometryFace( grid, face, 
			    local_nnode, local_xyz, local_uv, 
			    local_nface, local_f2n );

    free(local_xyz);
    free(local_uv);
    free(local_f2n);

    free( facel2g );
  }
  free( faceg2l );

  if ( !UGrid_FromArrays( &ugrid, nnode, xyz, nface, f2n, ncell, c2n  )) {
    printf(" Could not make UGridPtr, line %d of %s\n", __LINE__, __FILE__);
    return(-1);
  }
  
  for (iface = 0 ; iface < nface ; iface++ ) {
    UGrid_FlagValue(ugrid,iface) = faceId[iface];
  }

  printf("Rebuilding Element Connectivity...");
  UGrid_BuildConnectivity(ugrid);		/* Build Connectivity */
  printf("Complete\n");

  UGrid_TIMESTAMP(ugrid) = time( NULL );	/* Updated time */
  UGrid_ALGORITHM(ugrid) = UGrid_ALGORITHM(CADGeom_VolumeGrid(vol));

  if( !CADGeom_SetVolumeGrid( vol, ugrid ) ) {
    printf(" Could not replace CADGeom volume grid, line %d of %s\n",
	   __LINE__, __FILE__);    
    return(-1);
  }

  /* get uv vals for surface(s) */
  /* we use globalPatch to track with the localPatch so that we can get global
   * node numbering relative the volume grid and NOT the face grid as would
   * be the case of global index of upp
   */

  globalPatch = DList_SetIteratorToHead(UGrid_PatchList(ugrid),&patchIterator);

  for( face=1; face<=nGeomFace; face++ ) {
    double uv[2];
    char filename[265];

    localPatch = CADGeom_FaceGrid(vol,face);
    sprintf(filename,"localUGPatch%03d.t",face);
    sprintf(filename,"globalUGPatch%03d.t",face);
    UGPatch_GetDims(localPatch,patchDimensions);
    for( localNode=0; localNode<patchDimensions[0]; localNode++ ) {	     
      globalNode = UGPatch_GlobalIndex(globalPatch,localNode);
      gridNodeUV( grid, globalNode, face, uv );
      UGPatch_Parameter(localPatch,localNode,0) = uv[0];
      UGPatch_Parameter(localPatch,localNode,1) = uv[1];
    }

    globalPatch = DList_GetNextItem(&patchIterator);
  }

  if ( NULL != project ) {
    GeoMesh_UseDefaultIOCallbacks();
    if( !CADGeom_SavePart(vol,project) ) {
      printf("Yo! Could NOT save \"%s\".\n",project);
    }
  }

  return 0;
}

 
