
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include "mesherx.c"
#include "Goolache/CAPrI_IO.h"
#include "CADGeom/CADGeom.h"
#include "UG_API/UGMgr.h"

#ifdef PROE_MAIN
int GridEx_Main( int argc, char *argv[] )
#else
int main( int argc, char *argv[] )
#endif
{

  char   project[256];
  int    vol=1;
  int    i;
  int    npo;
  int    nel;
  int    *iel;
  double *xyz;
  double scale;
  int maxnode;
  GridBool mixedElement, blendElement, qualityImprovement, copyGridY;
  GridBool bil;
  char     modeler[81];

  sprintf( modeler,       "" );
  sprintf( project,       "" );
  scale = 1.0;
  maxnode = 50000;
  mixedElement = FALSE;
  blendElement = FALSE;
  qualityImprovement = FALSE;
  copyGridY = FALSE;
  bil = FALSE;

  i = 1;
  while( i < argc ) {
    if( strcmp(argv[i],"-p") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-p argument %d: %s\n",i, project);
    } else if( strcmp(argv[i],"-felisa") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-felisa argument %d: %s\n",i, project);
      sprintf( modeler, "FELISA" );
    } else if( strcmp(argv[i],"-open") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-open argument %d: %s\n",i, project);
      sprintf( modeler, "OpenCASCADE" );
    } else if( strcmp(argv[i],"-proe") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-proe argument %d: %s\n",i, project);
      sprintf( modeler, "Pro/ENGINEER" );
    } else if( strcmp(argv[i],"-parasolid") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-parasolid argument %d: %s\n",i, project);
      sprintf( modeler, "Parasolid" );
    } else if( strcmp(argv[i],"-catia") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-catia argument %d: %s\n",i, project);
      sprintf( modeler, "CatiaV5" );
    } else if( strcmp(argv[i],"-ug") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-ug argument %d: %s\n",i, project);
      sprintf( modeler, "UniGraphics" );
    } else if( strcmp(argv[i],"-sw") == 0 ) {
      i++; sprintf( project, "%s", argv[i] );
      printf("-sw argument %d: %s\n",i, project);
      sprintf( modeler, "SolidWorks" );
    } else if( strcmp(argv[i],"-s") == 0 ) {
      i++; scale = atof(argv[i]);
      printf("-s argument %d: %f\n",i, scale);
    } else if( strcmp(argv[i],"-n") == 0 ) {
      i++; maxnode = atoi(argv[i]);
      printf("-n argument %d: %d\n",i, maxnode);
    } else if( strcmp(argv[i],"-m") == 0 ) {
      mixedElement = TRUE;
      printf("-m argument %d: activated mixed element layers\n",i);
    } else if( strcmp(argv[i],"-q") == 0 ) {
      qualityImprovement = TRUE;
      printf("-q argument %d: activated grid quality improvement\n",i);
    } else if( strcmp(argv[i],"-b") == 0 ) {
      blendElement = TRUE;
      printf("-b argument %d: activated blend elements\n",i);
    } else if( strcmp(argv[i],"-cy") == 0 ) {
      copyGridY = TRUE;
      printf("-cy argument %d: activated grid copy about y=0 elements\n",i);
    } else if( strcmp(argv[i],"-bil") == 0 ) {
      bil = TRUE;
      printf("-bil argument %d: activated Bil Kleb's case\n",i);
    } else if( strcmp(argv[i],"-vgbg") == 0 ) {
      UGMgr_SetSizerFromIdentity( UG_VGRID );
      printf("-vgbg argument %d: activated VGRID background spacing\n",i);
    } else if( strcmp(argv[i],"-h") == 0 ) {
      printf("Usage: flag value pairs:\n");
#ifdef HAVE_CAPRI2
      printf(" -felisa input FELISA project name\n");
      printf(" -open input OpenCASCADE project name\n");
      printf(" -proe input Pro/E project name\n");
      printf(" -parasolid input Parasolid project name\n");
      printf(" -catia input Catia V5 project name\n");
      printf(" -ug input Unigraphics project name\n");
      printf(" -sw input SolidWorks project name\n");
#else
      printf(" -p input project name\n");
#endif
      printf(" -s scale background grid\n");
      printf(" -n maximum number of nodes\n");
      printf(" -m mixed element layers\n");
      printf(" -q use edge swapping to improve grid quality\n");
      printf(" -b use blend elements on first layer\n");
      printf(" -cy copy grid about the y=0 plane\n");
      printf(" -vgbg set background spacing to VGRID\n");
      printf(" -bil Bil Kleb's case\n");
      return(0);
    } else {
      fprintf(stderr,"Argument \"%s %s\" Ignored\n",argv[i],argv[i+1]);
      i++;
    }
    i++;
  }

  if(strcmp(modeler,"")==0)       sprintf(modeler,"FELISA" );
  if(strcmp(project,"")==0)       sprintf(project,"../test/box1" );

  printf("calling MeshMgr_Initialize ... \n");
  if ( ! UGMgr_LoadLibs( ) ){
    printf("ERROR: UGMgr_LoadLibs broke.\n%s\n",ErrMgr_GetErrStr());
    return 1;
  }  

  printf("calling CADGeom_Start ... \n");
  if ( ! CADGeom_Start( ) ){
    printf("ERROR: CADGeom_Start broke.\n%s\n",ErrMgr_GetErrStr());
    return 1;
  }  

  printf("calling GeoMesh_Load for project <%s> ... \n",project);
#ifdef HAVE_CAPRI2
  if ( ! GeoMesh_LoadPart( modeler, project ) ){
#else
  if ( ! GeoMesh_LoadPart( project ) ){
#endif
    printf("ERROR: GeoMesh_LoadPart broke.\n%s\n",ErrMgr_GetErrStr());
    return 1;
  }

  if ( scale != 1.0 ) {
    UGMgr_SetSizerScale( scale );
    if ( !CAPrIMesh_CreateTShell( vol )) {
      printf("ERROR: could not create shell\n");
      return 0;
    }
  }else{
    if (NULL == CADGeom_VolumeGrid(vol) ) {
      if ( !CAPrIMesh_CreateTShell( vol )) {
	printf("ERROR: could not create shell\n");
	return 0;
      }   
    }
  }


  MesherX_DiscretizeVolume( maxnode, scale, project, 
			    mixedElement, blendElement, qualityImprovement,
			    copyGridY, bil );

  return(0);
}
    
