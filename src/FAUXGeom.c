
/* Karen Bibb
 * Aerothermodynamics Branch
 * NASA Langley Research Center
 * Phone:(757)864-8005
 * Email:k.l.bibb@larc.nasa.gov 
 */

/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* a more general FAKEGeom
 * requires file faux_input
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include "FAKEGeom.h"

#ifdef __APPLE__       /* Not needed on Mac OS X */
#else
#include <malloc.h>
#endif

#include "FAKEGeom.h"

/*  original definitions, using general face definition.
    
static double x0 = -135.02;
static double x1 =  405.61;
static double y0 =    0.00;
static double y1 =  231.00;
static double z0 = -231.00;
static double z1 =  231.00;

#define facex0 (3)
#define facex1 (4)
#define facey0 (2)
#define facey1 (1)
#define facez0 (6)
#define facez1 (5)
*/
int nfaux = 0;

struct face {
  int faceid;
  char faceType[10];
  double normal(3);
  double offset;
};

face *faux_faces = NULL;

int *faceId2Index = NULL;

static void initialize_faux(void)
{

  FILE *f;
  f= fopen("faux_input","r");
  if (NULL==f){
    printf("could not open faux_input file\n");
    return;
  }
  fscanf("%d",&nfaux);
  faux_faces = (face *) malloc( nfaux * sizeof(face) );
  for (i=0;i<nfaux;i++){
    faux_faces[i].normal[0] = 0.0;
    faux_faces[i].normal[1] = 0.0;
    faux_faces[i].normal[2] = 0.0;
    fscanf("%d %s %ld",&(faux_faces[i].faceid),flavor,&(faux_faces[i].offset));
    if( strcmp(flavor,"xmin") == 0 ) { faux_faces[i].normal[0] = 1.0;
      faux_faces[i].faceType = "xplane"};
    if( strcmp(flavor,"xmax") == 0 ) { faux_faces[i].normal[0] = -1.0;
      faux_faces[i].faceType = "xplane"};
    if( strcmp(flavor,"ymin") == 0 ) { faux_faces[i].normal[1] = 1.0;
      faux_faces[i].faceType = "yplane"};
    if( strcmp(flavor,"ymax") == 0 ) { faux_faces[i].normal[1] = -1.0;
      faux_faces[i].faceType = "yplane"};
    if( strcmp(flavor,"zmin") == 0 ) { faux_faces[i].normal[2] = 1.0;
      faux_faces[i].faceType = "zplane"};
    if( strcmp(flavor,"zmax") == 0 ) { faux_faces[i].normal[2] = -1.0;
      faux_faces[i].faceType = "zplane"};
  }
}

GridBool CADGeom_NearestOnEdge(int vol, int edgeId, 
			   double *xyz, double *t, double *xyznew)
{
  *t = 0;
  xyznew[0] = xyz[0];
  xyznew[1] = xyz[1];
  xyznew[2] = xyz[2];

  return TRUE;
}

GridBool CADGeom_NearestOnFace(int vol, int faceId, 
			   double *xyz, double *uv, double *xyznew)
{

  if (0 == nfaux) initialize_faux( );

  switch (faux_faces[i].faceType) {
  case "xplane":
    uv[0] = xyz[1];
    uv[1] = xyz[2];
    xyznew[0] = faux_faces[i].offset;
    xyznew[1] = xyz[1];
    xyznew[2] = xyz[2];
    break;
  case "yplane":
    uv[0] = xyz[0];
    uv[1] = xyz[2];
    xyznew[0] = xyz[0];
    xyznew[1] = faux_faces[i].offset;
    xyznew[2] = xyz[2];
    break;
  case "zplane":
    uv[0] = xyz[0];
    uv[1] = xyz[1];
    xyznew[0] = xyz[0];
    xyznew[1] = xyz[1];
    xyznew[2] = faux_faces[i].offset;
    break;
  default:
    printf("ERROR: %s: %d: face %d unknown.\n",__FILE__,__LINE__,faceId);
    return FALSE;
  }


  return TRUE;
}

GridBool CADGeom_PointOnEdge(int vol, int edgeId, 
			 double t, double *xyz, 
			 int derivativeFlag, double *dt, double *dtdt )
{
  xyz[0] = t;
  xyz[1] = 0.0;
  xyz[2] = 0.0;

  if (derivativeFlag > 0){
    dt[0] = 1.0;
    dt[1] = 0.0;
    dt[2] = 0.0;
  }

  if (derivativeFlag > 1){
    dtdt[0] = 0.0;
    dtdt[1] = 0.0;
    dtdt[2] = 0.0;
  }

  return TRUE;
}

GridBool CADGeom_PointOnFace(int vol, int faceId, 
			 double *uv, double *xyz, 
			 int derivativeFlag, double *du, double *dv,
			 double *dudu, double *dudv, double *dvdv )
{
  switch (faux_faces[i].faceType) {
  case "xplane":
    xyz[0] = faux_faces[i].offset;
    xyz[1] = uv[0];
    xyz[2] = uv[1];
    break;
  case "yplane":
    xyz[0] = uv[0];
    xyz[1] = faux_faces[i].offset;
    xyz[2] = uv[1];
    break;
  case "zplane":
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = faux_faces[i].offset;
    break;
    break;
  default:
    printf("ERROR: %s: %d: face %d unknown.\n",__FILE__,__LINE__,faceId);
    return FALSE;
  }

  if (derivativeFlag > 0){
    switch (faux_faces[i].faceType) {
    case "xplane":
      du[0] = 0.0; du[1] = 1.0; du[2] = 0.0;
      dv[0] = 0.0; dv[1] = 0.0; dv[2] = 1.0;
      break;
    case "yplane":
      du[0] = 1.0; du[1] = 0.0; du[2] = 0.0;
      dv[0] = 0.0; dv[1] = 0.0; dv[2] = 1.0;
      break;
    case "zplane":
      du[0] = 1.0; du[1] = 0.0; du[2] = 0.0;
      dv[0] = 0.0; dv[1] = 1.0; dv[2] = 0.0;
      break;
    default:
      printf("ERROR: %s: %d: face %d unknown.\n",__FILE__,__LINE__,faceId);
      return FALSE;
    }
  }

  if (derivativeFlag > 1){
    dudu[0] = 0.0;
    dudu[1] = 0.0;
    dudu[2] = 0.0;
    dudv[0] = 0.0;
    dudv[1] = 0.0;
    dudv[2] = 0.0;
    dvdv[0] = 0.0;
    dvdv[1] = 0.0;
    dvdv[2] = 0.0;
  }

  return TRUE;
}


GridBool CADGeom_NormalToFace(int vol, int faceId, 
			 double *uv, double *xyz, double *normal )
{
  /* NOTE:  negative on normal, used to match existing FAKEGeom to normals assigned at beginning.
   *        Initialization seems to be correct, but this is the direction that has been working...
   */
  normal[0] = -faux_faces[i].normal[0];
  normal[1] = -faux_faces[i].normal[2];
  normal[2] = -faux_faces[i].normal[2];
  
  switch (faux_faces[i].faceType) {
  case "xplane":
    xyz[0] = faux_faces[i].offset;
    xyz[1] = uv[0];
    xyz[2] = uv[1];
  case "yplane":
    xyz[0] = uv[0];
    xyz[1] = faux_faces[i].offset;
    xyz[2] = uv[1];
  case "zplane":
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = faux_faces[i].offset;
    break;
  default:
    printf("ERROR: %s: %d: face %d unknown.\n",__FILE__,__LINE__,faceId);
    return FALSE;
  }
}
