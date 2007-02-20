
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

int nfaux = 0;

#define xplane (0)
#define yplane (1)
#define zplane (2)

typedef struct face face;
struct face {
  int faceid;
  int faceType;
  double normal[3];
  double offset;
};

static struct face *faux_faces = NULL;

static GridBool initialize_faux(void)
{
  char flavor[1025];
  int i;
 
  FILE *f;
  f = fopen("faux_input","r");
  if (NULL==f){
    printf("could not open faux_input file\n");
    return FALSE;
  }
  if ( 1 != fscanf(f,"%d",&nfaux) ) 
    {
      printf("error parsing line 1 of faux_input file\n");
      nfaux = 0;
      return FALSE;      
    }
  faux_faces = (face *) malloc( nfaux * sizeof(face) );
  for (i=0;i<nfaux;i++){
    faux_faces[i].normal[0] = 0.0;
    faux_faces[i].normal[1] = 0.0;
    faux_faces[i].normal[2] = 0.0;
    if ( 3 != fscanf(f,"%d %s %lf",
		     &(faux_faces[i].faceid),
		     flavor,
		     &(faux_faces[i].offset) ) )
      {
	printf("error parsing line %d of faux_input file\n",i+2);
	nfaux = 0;
	free( faux_faces );
	faux_faces = NULL;
	return FALSE;
      }

    if(      strcmp(flavor,"xplane") == 0 ) 
      { 
	faux_faces[i].normal[0] = 1.0;
	faux_faces[i].faceType = xplane; 
      }
    else if( strcmp(flavor,"yplane") == 0 ) 
      { 
	faux_faces[i].normal[1] = 1.0;
	faux_faces[i].faceType = yplane; 
      }
    else if( strcmp(flavor,"zplane") == 0 ) 
      { 
	faux_faces[i].normal[2] = 1.0;
	faux_faces[i].faceType = zplane; 
      }
    else
      {
	printf("error parsing flavor %s on line %d of faux_input file\n",
	       flavor,i+2);
	nfaux = 0;
	free( faux_faces );
	faux_faces = NULL;
	return FALSE;	   
      }

    printf("%4d: %4d of %4d type %3d offset %30.15f\n",
	   i, faux_faces[i].faceid, nfaux, 
	   faux_faces[i].faceType, faux_faces[i].offset);
   
  }

  fclose(f);
  return TRUE;
}

static int faux_faceId( int faceId )
{
  int i;

  for (i=0;i<nfaux;i++){
    if( faux_faces[i].faceid == faceId ){
      return i;
    }
  }

  printf("ERROR: %s: %d: face %d unknown in FAUXGeom.\n",
	 __FILE__,__LINE__,faceId);
  return -1;

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
  int id;

  if (0 == nfaux) { 
    if ( !initialize_faux( ) )
      {
	return FALSE;	
      }
  }

  id = faux_faceId(faceId);
  if ( id < 0  ) return FALSE;	
    
  switch (faux_faces[id].faceType) {
  case xplane:
    uv[0] = xyz[1];
    uv[1] = xyz[2];
    xyznew[0] = faux_faces[id].offset;
    xyznew[1] = xyz[1];
    xyznew[2] = xyz[2];
    break;
  case yplane:
    uv[0] = xyz[0];
    uv[1] = xyz[2];
    xyznew[0] = xyz[0];
    xyznew[1] = faux_faces[id].offset;
    xyznew[2] = xyz[2];
    break;
  case zplane:
    uv[0] = xyz[0];
    uv[1] = xyz[1];
    xyznew[0] = xyz[0];
    xyznew[1] = xyz[1];
    xyznew[2] = faux_faces[id].offset;
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
  int id;

  if (0 == nfaux) { 
    if ( !initialize_faux( ) )
      {
	return FALSE;	
      }
  }

  id = faux_faceId(faceId);
  if ( id < 0 ) return FALSE;	

  switch (faux_faces[id].faceType) {
  case xplane:
    xyz[0] = faux_faces[id].offset;
    xyz[1] = uv[0];
    xyz[2] = uv[1];
    break;
  case yplane:
    xyz[0] = uv[0];
    xyz[1] = faux_faces[id].offset;
    xyz[2] = uv[1];
    break;
  case zplane:
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = faux_faces[id].offset;
    break;
    break;
  default:
    printf("ERROR: %s: %d: face %d unknown.\n",__FILE__,__LINE__,faceId);
    return FALSE;
  }

  if (derivativeFlag > 0){
    switch (faux_faces[id].faceType) {
    case xplane:
      du[0] = 0.0; du[1] = 1.0; du[2] = 0.0;
      dv[0] = 0.0; dv[1] = 0.0; dv[2] = 1.0;
      break;
    case yplane:
      du[0] = 1.0; du[1] = 0.0; du[2] = 0.0;
      dv[0] = 0.0; dv[1] = 0.0; dv[2] = 1.0;
      break;
    case zplane:
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
  int id;

  if (0 == nfaux) { 
    if ( !initialize_faux( ) )
      {
	return FALSE;	
      }
  }

  id = faux_faceId(faceId);
  if ( id < 0 ) return FALSE;	

  normal[0] = faux_faces[id].normal[0];
  normal[1] = faux_faces[id].normal[2];
  normal[2] = faux_faces[id].normal[2];
  
  switch (faux_faces[id].faceType) {
  case xplane:
    xyz[0] = faux_faces[id].offset;
    xyz[1] = uv[0];
    xyz[2] = uv[1];
  case yplane:
    xyz[0] = uv[0];
    xyz[1] = faux_faces[id].offset;
    xyz[2] = uv[1];
  case zplane:
    xyz[0] = uv[0];
    xyz[1] = uv[1];
    xyz[2] = faux_faces[id].offset;
    break;
  default:
    printf("ERROR: %s: %d: face %d unknown.\n",__FILE__,__LINE__,faceId);
    return FALSE;
  }
}
