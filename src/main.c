
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
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
  printf("new size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));
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
  int gridDimensions[3];
  int nnode, nface, ncell;
  int maxnode, maxface, maxcell;
  int i;
  double *xyz;
  int *c2n, *f2n, *faceId;

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

  printf("ugrid size: %d nodes %d faces %d cells.\n",nnode,nface,ncell);

  maxnode = nnode * 100;
  maxface = nface * 100;
  maxcell = ncell * 100;

  printf("max grid size: %d nodes %d faces %d cells.\n",
	 maxnode,maxface,maxcell);

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

  return gridImport( maxnode,nnode, maxface, nface, maxcell, ncell,
		     xyz, f2n, faceId, c2n );
}
 
int gridSavePart( Grid *grid, char *project )
{
  int vol=1;
  UGridPtr ugrid;  
  int nnode, nface, ncell;
  double *xyz;
  int *f2n, *faceId, *c2n;
  int iface; 

  gridExport( grid, &nnode, &nface, &ncell,
	      &xyz, &f2n, &faceId, &c2n );

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

 
