
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#include <stdlib.h>
#include <stdio.h>

#ifdef __APPLE__       /* Not needed on Mac OS X */
#else
#include <malloc.h>
#endif

#include "FAKEGeom.h"

int nfaux = 0;

struct face {
  int faceid;
  double normal(3);
  double offset;
};

face *faux_faces = NULL;

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
    if( strcmp(flavor,"xmin") == 0 ) faux_faces[i].normal[0] = 1.0;
    if( strcmp(flavor,"xmax") == 0 ) faux_faces[i].normal[0] = -1.0;
    if( strcmp(flavor,"ymin") == 0 ) faux_faces[i].normal[1] = 1.0;
    if( strcmp(flavor,"ymax") == 0 ) faux_faces[i].normal[1] = -1.0;
    if( strcmp(flavor,"zmin") == 0 ) faux_faces[i].normal[2] = 1.0;
    if( strcmp(flavor,"zmax") == 0 ) faux_faces[i].normal[2] = -1.0;
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
  
  switch (faceId) {
  case 3:
    uv[0] = xyz[1];
    uv[1] = xyz[2];
    xyznew[0] = x0;
    xyznew[1] = xyz[1];
    xyznew[2] = xyz[2];
    break;
  case 4:
    uv[0] = xyz[1];
    uv[1] = xyz[2];
    xyznew[0] = x1;
    xyznew[1] = xyz[1];
    xyznew[2] = xyz[2];
    break;
  case 28:
    uv[0] = xyz[0];
    uv[1] = xyz[2];
    xyznew[0] = xyz[0];
    xyznew[1] = y0;
    xyznew[2] = xyz[2];
    break;
  case 2:
    uv[0] = xyz[0];
    uv[1] = xyz[2];
    xyznew[0] = xyz[0];
    xyznew[1] = y1;
    xyznew[2] = xyz[2];
    break;
  case 29: 
    uv[0] = xyz[0];
    uv[1] = xyz[1];
    xyznew[0] = xyz[0];
    xyznew[1] = xyz[1];
    xyznew[2] = z0;
    break;
  case 1: 
    uv[0] = xyz[0];
    uv[1] = xyz[1];
    xyznew[0] = xyz[0];
    xyznew[1] = xyz[1];
    xyznew[2] = z1;
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
  switch (faceId) {
  case 3:
    xyz[0] = x0;
    xyz[1] = uv[0];
    xyz[2] = uv[1];
    break;
  case 4:
    xyz[0] = x1;
    xyz[1] = uv[0];
    xyz[2] = uv[1];
    break;
  case 28:
    xyz[0] = uv[0];
    xyz[1] = y0;
    xyz[2] = uv[1];
    break;
  case 2:
    xyz[0] = uv[0];
    xyz[1] = y1;
    xyz[2] = uv[1];
    break;
  case 29: 
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = z0;
    break;
  case 1: 
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = z1;
    break;
  default:
    printf("ERROR: %s: %d: face %d unknown.\n",__FILE__,__LINE__,faceId);
    return FALSE;
  }

  if (derivativeFlag > 0){
    switch (faceId) {
    case 3: case 4:
      du[0] = 0.0; du[1] = 1.0; du[2] = 0.0;
      dv[0] = 0.0; dv[1] = 0.0; dv[2] = 1.0;
      break;
    case 2: case 28:
      du[0] = 1.0; du[1] = 0.0; du[2] = 0.0;
      dv[0] = 0.0; dv[1] = 0.0; dv[2] = 1.0;
      break;
    case 1: case 29:
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
  switch (faceId) {
  case 3:
    xyz[0] = x0;
    xyz[1] = uv[0];
    xyz[2] = uv[1];
    normal[0] = -1.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    break;
  case 4:
    xyz[0] = x1;
    xyz[1] = uv[0];
    xyz[2] = uv[1];
    normal[0] = 1.0;
    normal[1] = 0.0;
    normal[2] = 0.0;
    break;
  case 28:
    xyz[0] = uv[0];
    xyz[1] = y0;
    xyz[2] = uv[1];
    normal[0] = 0.0;
    normal[1] = -1.0;
    normal[2] = 0.0;
    break;
  case 2:
    xyz[0] = uv[0];
    xyz[1] = y1;
    xyz[2] = uv[1];
    normal[0] = 0.0;
    normal[1] = 1.0;
    normal[2] = 0.0;
    break;
  case 29: 
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = z0;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = -1.0;
    break;
  case 1: 
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = z1;
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 1.0;
    break;
  default:
    printf("ERROR: %s: %d: face %d unknown.\n",__FILE__,__LINE__,faceId);
    return FALSE;
  }
}
