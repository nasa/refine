
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
#include <values.h>
#include "grid.h"
#include "gridmetric.h"
#include "gridswap.h"
#include "gridcad.h"
#include "layer.h"
#include <CADGeom/CADGeom.h>

Grid *gridLoadPart( char *project, int maxnode );
int gridSavePart( Grid *grid, char *project );

#define PRINT_STATUS printf("minimum Thawed Aspect Ratio %8.6f Mean Ratio %8.6f Volume %10.6e\n", gridMinThawedAR(grid),gridMinThawedFaceMR(grid), gridMinVolume(grid));

#define DUMP_TEC if (!boundaryLayerGrid) {iview++;printf("Frame %d\n",iview);gridWriteTecplotSurfaceZone(grid);}

#define STATUS DUMP_TEC PRINT_STATUS

int main( int argc, char *argv[] )
{
  Grid *grid;
  Layer *layer;
  int bcs[2], jmax;
  double height;
  char project[256];
  char adaptfile[256], outputProject[256], outputFAST[256];
  int i, j, oldSize, newSize;
  int wiggleSteps, wiggle;
  double ratio=0.3;
  double ratioRefine, ratioCollapse;
  bool projected;
  bool boundaryLayerGrid = FALSE;
  bool debugInsert = FALSE;
  int iview = 0;

  sprintf( project,       "" );
  sprintf( outputProject, "" );
  sprintf( adaptfile,     "" );    
  sprintf( outputFAST,    "" );

  i = 1;
  while( i < argc ) {
    if( strcmp(argv[i],"-p") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-p argument %d: %s\n",i, project);
    } else if( strcmp(argv[i],"-o") == 0 ) {
      i++; sprintf( outputProject, "%s", argv[i] );
      printf("-o argument %d: %s\n",i, outputProject);
    } else if( strcmp(argv[i],"-a") == 0 ) {
      i++; sprintf( adaptfile,"%s",argv[i]  );
      printf("-a argument %d: %s\n",i, adaptfile);
    } else if( strcmp(argv[i],"-l") == 0 ) {
      boundaryLayerGrid = TRUE;
      printf("-l argument %d, ignoring -a \n",i);
    } else if( strcmp(argv[i],"-r") == 0 ) {
      i++; ratio = atof(argv[i]);
      printf("-r argument %d: %f\n",i, ratio);
    } else if( strcmp(argv[i],"-i") == 0 ) {
      debugInsert = TRUE;
      printf("-i argument %d\n",i);
    } else if( strcmp(argv[i],"-h") == 0 ) {
      printf("Usage: flag value pairs:\n");
      printf(" -p input project name\n");
      printf(" -o output project name\n");
      printf(" -a party project_adapt_hess file name\n");
      printf(" -l make a boundary layer grid -a ignored\n");
      printf(" -r initial edge length ratio for adapt\n");
      printf(" -i debug general insert nodes\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }
  
  if(debugInsert)                 sprintf(project,"../test/box1" );
  if(strcmp(project,"")==0)       sprintf(project,"../test/om6" );
  if(strcmp(outputProject,"")==0) sprintf(outputProject,"%s_out", project );
  if(strcmp(adaptfile,"")==0)     sprintf(adaptfile,"%s_adapt_hess",project);
  if(strcmp(outputFAST,"")==0)    sprintf(outputFAST,"%s.fgrid",outputProject);

  if(boundaryLayerGrid || debugInsert ) sprintf(adaptfile,"none");

  printf("running project %s\n",project);
  grid = gridLoadPart( project, 500000 );

  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: loaded part does not have right handed boundaries\n");

  printf("restart grid size: %d nodes %d faces %d cells.\n",
	 gridNNode(grid),gridNFace(grid),gridNCell(grid));

  if(strcmp(adaptfile,"none")==0) {
    printf("adapt parameter >none< selected.\n");
    gridResetSpacing(grid);
    if (boundaryLayerGrid) {
      printf("freezing distant volume nodes.\n");
      gridFreezeAll(grid);
      printf("thaw bc 1.\n");
      gridThawNearBC(grid,0.5,1);
      printf("thaw bc 2.\n");
      gridThawNearBC(grid,0.5,2);
      gridFreezeBCFace(grid,1);
      gridFreezeBCFace(grid,2);
      printf("make advancing layer object.\n");
      layer = layerCreate(grid);
      bcs[0]=1;
      bcs[1]=2;
      printf("make advancing layer front.\n");
      layerMakeFront(layer,2,bcs);
      printf("make advancing layer front normals.\n");
      layerMakeNormal(layer);
      layerConstrainNormal(layer,5);
      printf("make advancing layer front normals visible to front.\n");
      layerVisibleNormals(layer);
    }else{
      if (!debugInsert)
	gridScaleSpacingSphere(grid, 0.0, 0.0, 0.0, 1.0, 0.7 );
    }
  }else{
    printf("reading adapt parameter from file %s ...\n",adaptfile);
    gridImportAdapt(grid, adaptfile); // Do not sort nodes before this call.
  }
  STATUS;

  for (i=0;i<3;i++){
    projected = ( grid == gridRobustProject(grid));
    if (projected) {
      printf("edge swapping grid...\n");gridSwap(grid);
      printf("node smoothing grid...\n");gridSmooth(grid);
    }else{
      printf("node smoothing volume grid...\n");gridSmoothVolume(grid);
    }
  }
  STATUS;

  oldSize = 1;
  newSize = gridNNode(grid);
  jmax = 40;
  if (debugInsert) jmax = -1;
  for ( j=0; (j<jmax) && (
	(ratio < 0.99) || 
	  (((double)ABS(newSize-oldSize)/(double)oldSize)>0.001) ||
	  !projected );
	j++){

    if (boundaryLayerGrid) {
      height = 0.00001*pow(1.2,j);
      layerTerminateNormalWithSpacing(layer,height*5.);
      if (layerNActiveNormal(layer) == 0 ) jmax=0;
      printf("insert layer height = %f\n",height);
      wiggleSteps = MIN(5,(int)(height/0.001)+1);
      height = height / (double)wiggleSteps;
      layerVisibleNormals(layer);
      layerAdvance(layer,height);
      for (i=1;i<wiggleSteps;i++) {
	printf("edge swapping grid...\n");gridSwap(grid);
	printf("node smoothing grid...\n");gridSmooth(grid);
	printf("edge swapping grid...\n");gridSwap(grid);
	printf("node smoothing grid...\n");gridSmooth(grid);
	gridAdapt(grid,0.4,1.5);
	printf("edge swapping grid...\n");gridSwap(grid);
	printf("node smoothing grid...\n");gridSmooth(grid);
	printf("edge swapping grid...\n");gridSwap(grid);
	printf("node smoothing grid...\n");gridSmooth(grid);
	printf("wiggle step %d of %d, minAR %8.5f\n",i+1,wiggleSteps,gridMinThawedAR(grid));
	layerWiggle(layer,height);
	//printf("minimum Volume %12.8e\n", gridMinVolume(grid));
      }
      printf("edge swapping grid...\n");gridSwap(grid);
      printf("node smoothing grid...\n");gridSmooth(grid);
      printf("edge swapping grid...\n");gridSwap(grid);
      printf("node smoothing grid...\n");gridSmooth(grid);
      STATUS;
    }
    if (ratio<0.01) ratio = 0.01;
    if (ratio>1.0) ratio = 1.0;
    ratioCollapse = 0.4*ratio;
    ratioRefine   = 1.5/ratio;
    if (boundaryLayerGrid) {
      printf("adapt, ratio %4.2f, collapse limit %8.5f, refine limit %10.5f\n",
	     ratio, 0.4, 1.5 );
      gridAdapt(grid,0.4,1.5);
    }else{
      printf("adapt, ratio %4.2f, collapse limit %8.5f, refine limit %10.5f\n",
             ratio, ratioCollapse, ratioRefine );
      gridAdapt(grid,ratioCollapse,ratioRefine);
    }
    oldSize = newSize;
    newSize = gridNNode(grid) ;
    printf("%02d new size: %d nodes %d faces %d cells %d edge elements.\n",
	   j, gridNNode(grid),gridNFace(grid),gridNCell(grid),gridNEdge(grid));
    STATUS;
        
    for (i=0;i<2;i++){
      projected = ( grid == gridRobustProject(grid));
      if (projected) {
	printf("edge swapping grid...\n");gridSwap(grid);
	printf("node smoothing grid...\n");gridSmooth(grid);
	if (((double)ABS(newSize-oldSize)/(double)oldSize)<0.3)
	  ratio = ratio + 0.025;
      }else{
	printf("node smoothing volume grid...\n");gridSmoothVolume(grid);
	ratio = ratio - 0.05;
      }
    }
    STATUS;
    if (boundaryLayerGrid) {
    }else{
      gridFreezeGoodNodes(grid,0.6,0.4,1.5);
      printf("nodes frozen %d\n",gridNFrozen(grid));
    }
  }

  if (debugInsert) {
    gridInsertInToGeomEdge(grid, 0.532, 0.0, 0.0);
    gridInsertInToGeomEdge(grid, 1.0, 0.567, 0.0);
    gridInsertInToGeomEdge(grid, 1.0, 0.0, 0.33);
  }

  if (!gridRightHandedBoundary(grid)) 
    printf("ERROR: modifed grid does not have right handed boundaries\n");

  printf("writing output project %s\n",outputProject);
  gridSavePart( grid, outputProject );

  printf("writing output FAST file %s\n",outputFAST);
  gridExportFAST( grid, outputFAST );

  printf("Done.\n");

  return;
}

Grid *gridLoadPart( char *project, int maxnode )
{
  Grid *grid;
  int vol=1;
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

  maxface = maxnode;
  maxcell = maxnode * 6;
  maxedge = maxnode / 10;

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

  gridSetNGeomNode( grid, nGeomNode );
  gridSetNGeomEdge( grid, nGeomEdge );

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
      for( i=1 ; i < (nedgenode-2) ; i++ ) { // skip end segments  
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
    localPatch = CADGeom_FaceGrid(vol,face);
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
  int i, ixyz, iface, icell, inode, newnode, node;
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

  gridSortNodeGridEx( grid );

  gridExport( grid, &nnode, &nface, &ncell,
	      &xyz, &f2n, &faceId, &c2n );
  
  // edge stuff
  for (iedge=1; iedge<=nGeomEdge; iedge++){

    CADGeom_GetEdge( vol, iedge, trange, curveEndPoint);
    curveEndPoint[0]--; curveEndPoint[1]--;// fortran to c numbering
    
    nCurveNode = gridGeomCurveSize( grid, iedge, curveEndPoint[0]);
    curve =    malloc( nCurveNode *     sizeof(int) );
    temp_xyz = malloc( nCurveNode * 3 * sizeof(double) );
    temp_tuv = malloc( nCurveNode *     sizeof(double) );

    gridGeomCurve( grid, iedge, curveEndPoint[0], curve );
    gridGeomCurveT( grid, iedge, curveEndPoint[0], temp_tuv );

    for ( i=0; i<nCurveNode; i++) // include end points
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

  CADGeom_UseDefaultIOCallbacks();

  if( !CADGeom_SavePart(vol,project) ) {
    printf("Yo! Could NOT save \"%s\".\n",project);
  }

  printf("WARNING: skipping CADGeom_Stop ... \n");
  /*
  if ( ! CADGeom_Stop( ) ){
    printf("ERROR: CADGeom_Stop broke.\n%s\n",ErrMgr_GetErrStr());
    return NULL;
  }  
  */
}

 
