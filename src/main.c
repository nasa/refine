
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <values.h>
#include "grid.h"
#include <CADGeom/CADGeom.h>

Grid *gridLoadPart( char *project );
int gridSavePart( Grid *grid, char *project );

int main( int argc, char *argv[] )
{
  Grid *grid;
  char *filename;
  char *project;
  char *output;


  filename = "../test/om6_inv08.fgrid";
  printf("running default filename %s\n",filename);
  project = "../test/om6";
  printf("running default project %s\n",project);
  
  grid = gridImportFAST( filename );
  
  printf("orig size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));
  printf("minimum Aspect Ratio %12f\n",gridMinAR(grid));
  printf("minimum Volume %12.8e\n",gridMinVolume(grid));

  printf("edge swapping grid...\n");
  gridSwap(grid);

  printf("new size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));
  printf("minimum Aspect Ratio %12f\n",gridMinAR(grid));
  printf("minimum Volume %12.8e\n",gridMinVolume(grid));

  gridFree(grid);

  grid = gridLoadPart( project );

  printf("restart grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  printf("minimum Aspect Ratio %12f\n",gridMinAR(grid));
  printf("minimum Volume %12.8e\n",gridMinVolume(grid));

  printf("thrashing grid...\n");
  gridThrash(grid);
  printf("new size: %d nodes %d faces %d cells %d edge elements.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid),gridNEdge(grid));
  printf("minimum Aspect Ratio %12f\n",gridMinAR(grid));
  printf("minimum Volume %12.8e\n",gridMinVolume(grid));
  printf("edge swapping grid...\n");
  gridSwap(grid);

  printf("new size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));
  printf("minimum Aspect Ratio %12f\n",gridMinAR(grid));
  printf("minimum Volume %12.8e\n",gridMinVolume(grid));


  output = "../test/om6_out";
  printf("writing output project %s\n",output);

  gridSavePart( grid, output );

  printf("Done.\n");

  return;
}

Grid *gridLoadPart( char *project )
{
  Grid *grid;
  int vol=1;
  UGridPtr ugrid;
  CADCurvePtr edge;
  int gridDimensions[3];
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;
  int nedgenode;
  int nnode, nface, ncell, nedge;
  int maxnode, maxface, maxcell, maxedge;
  int i, iedge, inode;
  double *xyz;
  int *c2n, *f2n, *faceId;
  double trange[2];
  int edgeEndPoint[2];

  printf("calling CADGeom_Start ... \n");
  if ( ! CADGeom_Start( ) ){
    printf("ERROR: CADGeom_Start broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  

  printf("calling CADGeom_Load for project <%s> ... \n",project);
  if ( ! CADGeom_LoadPart( project ) ){
    printf("ERROR: CADGeom_LoadPart broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }

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

  maxnode = nnode * 100;
  maxface = nface * 100;
  maxcell = ncell * 100;
  maxedge = nedge * 100;

  printf("max grid size: %d nodes %d faces %d cells %d edge elements.\n",
	 maxnode,maxface,maxcell,maxedge);

  c2n    = malloc( maxcell * 4 * sizeof(int) );
  f2n    = malloc( maxface * 3 * sizeof(int) );
  faceId = malloc( maxface *     sizeof(int) );
  xyz    = malloc( maxnode * 3 * sizeof(double) );

  for( i=0; i<nnode; i++ ) {
    xyz[0+3*i] = UGrid_PtValue(ugrid,i,X);
    xyz[1+3*i] = UGrid_PtValue(ugrid,i,Y);
    xyz[2+3*i] = UGrid_PtValue(ugrid,i,Z);
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

  inode = nGeomNode;

  for( iedge=1; iedge<=nGeomEdge; iedge++ ) {
    if( (edge=CADGeom_EdgeGrid(vol,iedge)) == NULL ) 
      printf("ERROR: CADGeom_EdgeGrid(%d).\n%s\n",iedge,ErrMgr_GetErrStr());
 
    nedgenode = CADCURVE_NUMPTS(edge);

    CADGeom_GetEdge( vol, iedge, trange, edgeEndPoint );

    edgeEndPoint[0]--; /* convert from fortran to c numbers */
    edgeEndPoint[1]--;

    if (nedgenode == 2) {
      gridAddEdge(grid, edgeEndPoint[0], edgeEndPoint[1], iedge);
    }else{
      gridAddEdge(grid, edgeEndPoint[0], inode, iedge);
      for( i=0 ; i < (nedgenode-3) ; i++ ) {  
	gridAddEdge(grid, inode, inode+1, iedge);
	inode++;
      }
      gridAddEdge(grid, inode, edgeEndPoint[1], iedge);
    }
  }

  if ( nedge != gridNEdge(grid) )
    printf("ERROR: nedge != gridNEdge(grid)\n");

  return grid;
}
 
int gridSavePart( Grid *grid, char *project )
{
  int vol=1;
  UGridPtr ugrid;  
  int nnode, nface, ncell;
  double *xyz;
  int *f2n, *faceId, *c2n;
  int *o2n;
  int i, ixyz, iface, newnode;
  int iedge, curveEndPoint[2], nCurveNode, *curve;
  double trange[2];
  int patchDimensions[4]; // check on 4
  int nGeomNode, nGeomEdge, nGeomFace, nGeomGroups;

  double *temp_xyz, *temp_tuv;
  int *temp_face;

  Iterator    it;        /* DList Iterator */
  UGPatchPtr  patch;     /* UGPatch of Face relative to Volume */

  if( !CADGeom_GetVolume(vol,&nGeomNode,&nGeomEdge,&nGeomFace,&nGeomGroups) )
    printf("ERROR: CADGeom_GetVolume, line %d of %s\n.",__LINE__, __FILE__);

  gridExport( grid, &nnode, &nface, &ncell,
	      &xyz, &f2n, &faceId, &c2n );

  o2n = malloc( nnode * sizeof(int) );
  for (i=0;i<nnode;i++) o2n[i] = EMPTY;

  // geom nodes
  for (i=0;i<nGeomNode;i++) o2n[i] = i;
  newnode = nGeomNode;
  
  // edge stuff
  for (iedge=1; iedge<=nGeomEdge; iedge++){

    CADGeom_GetEdge( vol, iedge, trange, curveEndPoint);
    curveEndPoint[0]--; curveEndPoint[1]--;// fortran to c numbering
    
    nCurveNode = gridGeomCurveSize( grid, iedge, curveEndPoint[0]);
    curve =    malloc( nCurveNode *     sizeof(int) );
    temp_xyz = malloc( nCurveNode * 3 * sizeof(double) );
    temp_tuv = malloc( nCurveNode *     sizeof(double) );

    gridGeomCurve( grid, iedge, curveEndPoint[0], curve );

    for ( i=1; i<(nCurveNode-1); i++){ // skip end points
      o2n[curve[i]] = newnode;
      newnode++;
    }

    for ( i=0; i<nCurveNode; i++){ // include end points
      for ( ixyz=0; ixyz<3 ; ixyz++)
	temp_xyz[ixyz+3*i] = xyz[ixyz+3*curve[i]];
      temp_tuv[i] = DBL_MAX;
    }

    CADGeom_UpdateEdgeGrid( vol, iedge, nCurveNode, temp_xyz, temp_tuv );

    free(curve);
  }

  // face stuff goes here


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

  /* Face Stuff */

  patch = DList_SetIteratorToHead(UGrid_PatchList(ugrid),&it);

  for( iface=1; iface<=nGeomFace; iface++ ) {
    UGPatch_GetDims(patch,patchDimensions);
    temp_xyz = malloc( patchDimensions[0] * 3 * sizeof( double ) );
    temp_tuv = malloc( patchDimensions[0] * 2 * sizeof( double ) );
    temp_face = malloc( patchDimensions[1] * 3 * sizeof( int ) );
    for( i=0; i<patchDimensions[0]; i++ ) {
      temp_xyz[3*i+X] = UGPatch_PtValue(patch,i,X);
      temp_xyz[3*i+Y] = UGPatch_PtValue(patch,i,Y);
      temp_xyz[3*i+Z] = UGPatch_PtValue(patch,i,Z);
    }

    for( i=0; i<patchDimensions[1]; i++ ) {
      temp_face[3*i+0] =  UGPatch_VertValue(patch,i,0);
      temp_face[3*i+1] =  UGPatch_VertValue(patch,i,1);
      temp_face[3*i+2] =  UGPatch_VertValue(patch,i,2);
    }
    
    if( !CADGeom_UpdateFaceGrid(vol, iface, patchDimensions[0], temp_xyz, temp_tuv, patchDimensions[1], temp_face) ) {
      fprintf(stderr, "%s\nWarning face %d update in %s\n",ErrMgr_GetErrStr(),iface,__FILE__);
    }

    patch = (UGPatchPtr)DList_GetNextItem(&it);
  }

  CADGeom_UseDefaultIOCallbacks();

  if( !CADGeom_SavePart(vol,project) ) {
    printf("Yo! Could NOT save \"%s\".\n",project);
  }

  printf("calling CADGeom_Stop ... \n");
  if ( ! CADGeom_Stop( ) ){
    printf("ERROR: CADGeom_Stop broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  
  
}

 
