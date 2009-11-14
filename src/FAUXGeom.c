
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
#include <string.h>
#include <math.h>
#include "FAKEGeom.h"
#include "gridmath.h"

#ifdef __APPLE__       /* Not needed on Mac OS X */
#else
#include <malloc.h>
#endif

#include "FAKEGeom.h"

typedef struct face face;
struct face {
  int faceid;
  int faceType;
  double normal[3];
  double center[3];
  double offset;
  double u_dir[3];
  double v_dir[3];
};

#define nfaux_is_uninitialized (-1)

static int nfaux = nfaux_is_uninitialized;
static struct face *faux_faces = NULL;

#define xplane (0)
#define yplane (1)
#define zplane (2)
#define general_plane  (3)
#define cylinder  (4)

static GridBool initialize_faux(void)
{
  char flavor[1025];
  int i;
  int smallest_dir;
 
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

  if ( nfaux < 0 )
    {
      printf("line 1 of faux_input file negative %d, ignoring\n",nfaux);
      nfaux = 0;
      return FALSE;      
    }

  if ( 0 == nfaux ) return TRUE;

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
	printf("error parsing item %d of faux_input file\n",i+1);
	nfaux = 0;
	free( faux_faces );
	faux_faces = NULL;
	return FALSE;
      }

    if(      strncmp(flavor,"xplane",6) == 0 ) 
      { 
	faux_faces[i].normal[0] = 1.0;
	faux_faces[i].faceType = xplane; 
      }
    else if( strncmp(flavor,"yplane",6) == 0 ) 
      { 
	faux_faces[i].normal[1] = 1.0;
	faux_faces[i].faceType = yplane; 
      }
    else if( strncmp(flavor,"zplane",6) == 0 ) 
      { 
	faux_faces[i].normal[2] = 1.0;
	faux_faces[i].faceType = zplane; 
      }
    else if( strncmp(flavor,"general_plane",13) == 0 ) 
      { 
	faux_faces[i].faceType = general_plane; 
	if ( 3 != fscanf(f,"%lf %lf %lf",
			 &(faux_faces[i].normal[0]),
			 &(faux_faces[i].normal[1]),
			 &(faux_faces[i].normal[2]) ) )
	  {
	    printf("error parsing item %d of faux_input file\n",i+1);
	    printf("missing general_plane normal\n");
	    nfaux = 0;
	    free( faux_faces );
	    faux_faces = NULL;
	    return FALSE;
	  }
      }
    else if( strncmp(flavor,"cylinder",8) == 0 ) 
      { 
	faux_faces[i].faceType = general_plane; 
	if ( 3 != fscanf(f,"%lf %lf %lf",
			 &(faux_faces[i].center[0]),
			 &(faux_faces[i].center[1]),
			 &(faux_faces[i].center[2]) ) )
	  {
	    printf("error parsing line %d of faux_input file\n",i+1);
	    printf("missing first point of cylinder\n");
	    nfaux = 0;
	    free( faux_faces );
	    faux_faces = NULL;
	    return FALSE;
	  }
	if ( 3 != fscanf(f,"%lf %lf %lf",
			 &(faux_faces[i].normal[0]),
			 &(faux_faces[i].normal[1]),
			 &(faux_faces[i].normal[2]) ) )
	  {
	    printf("error parsing line %d of faux_input file\n",i+1);
	    printf("missing second point of cylinder\n");
	    nfaux = 0;
	    free( faux_faces );
	    faux_faces = NULL;
	    return FALSE;
	  }
	faux_faces[i].normal[0] -= faux_faces[i].center[0];
	faux_faces[i].normal[1] -= faux_faces[i].center[1];
	faux_faces[i].normal[2] -= faux_faces[i].center[2];
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

    /* make sure normal is normalized */
    gridVectorNormalize(faux_faces[i].normal);

    /* find the smallest normal component */
    smallest_dir = 0;
    if (ABS(faux_faces[i].normal[1]) < 
	ABS(faux_faces[i].normal[smallest_dir])) smallest_dir = 1;
    if (ABS(faux_faces[i].normal[2]) < 
	ABS(faux_faces[i].normal[smallest_dir])) smallest_dir = 2;
    /* choose u direction in the smallest normal component */
    faux_faces[i].u_dir[0] = 0.0;
    faux_faces[i].u_dir[1] = 0.0;
    faux_faces[i].u_dir[2] = 0.0;
    faux_faces[i].u_dir[smallest_dir] = 1.0;
    /* orthogonalize u_dir to normal */
    gridVectorOrthogonalize(faux_faces[i].u_dir, faux_faces[i].normal);
    /* v_dir is orthogonal to u_dir and normal */
    gridCrossProduct(faux_faces[i].normal,faux_faces[i].u_dir,
		     faux_faces[i].v_dir);
    /* make sure v_dir is normal */
    gridVectorNormalize(faux_faces[i].v_dir);

    /*
    printf("%4d: %4d of %4d type %3d offset %15.8f\n",
	   i, faux_faces[i].faceid, nfaux, 
	   faux_faces[i].faceType, faux_faces[i].offset);
    printf("%4d: normal %15.8f %15.8f %15.8f\n", i, 
	   faux_faces[i].normal[0],
	   faux_faces[i].normal[1],
	   faux_faces[i].normal[2]);
    printf("%4d: u_dir  %15.8f %15.8f %15.8f\n", i, 
	   faux_faces[i].u_dir[0],
	   faux_faces[i].u_dir[1],
	   faux_faces[i].u_dir[2]);
    printf("%4d: v_dir  %15.8f %15.8f %15.8f\n", i, 
	   faux_faces[i].v_dir[0],
	   faux_faces[i].v_dir[1],
	   faux_faces[i].v_dir[2]);
    */

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

  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(edgeId);

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
  double dxyz[3];
  double x, y;

  SUPRESS_UNUSED_COMPILER_WARNING(vol);

  if (nfaux_is_uninitialized == nfaux) { 
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
  case general_plane:
    uv[0] = faux_faces[id].u_dir[0]*xyz[0]
          + faux_faces[id].u_dir[1]*xyz[1]
          + faux_faces[id].u_dir[2]*xyz[2];
    uv[1] = faux_faces[id].v_dir[0]*xyz[0]
          + faux_faces[id].v_dir[1]*xyz[1]
          + faux_faces[id].v_dir[2]*xyz[2];
    xyznew[0] = faux_faces[id].offset * faux_faces[id].normal[0];
    xyznew[1] = faux_faces[id].offset * faux_faces[id].normal[1];
    xyznew[2] = faux_faces[id].offset * faux_faces[id].normal[2];
    xyznew[0] += uv[0] * faux_faces[id].u_dir[0];
    xyznew[1] += uv[0] * faux_faces[id].u_dir[1];
    xyznew[2] += uv[0] * faux_faces[id].u_dir[2];
    xyznew[0] += uv[1] * faux_faces[id].v_dir[0];
    xyznew[1] += uv[1] * faux_faces[id].v_dir[1];
    xyznew[2] += uv[1] * faux_faces[id].v_dir[2];
    break;
  case cylinder:
    /* vector to the surface point from center */
    dxyz[0] = xyz[0] - faux_faces[id].center[0];
    dxyz[1] = xyz[1] - faux_faces[id].center[1];
    dxyz[2] = xyz[2] - faux_faces[id].center[2];
    /* "z" direction along the center axis */
    uv[0] = dxyz[0]*faux_faces[id].normal[0]
          + dxyz[1]*faux_faces[id].normal[1]
          + dxyz[2]*faux_faces[id].normal[2];
    /* put dxyz in the plane normal to center axis */
    dxyz[0] -= uv[0]*faux_faces[id].normal[0];
    dxyz[1] -= uv[0]*faux_faces[id].normal[1];
    dxyz[2] -= uv[0]*faux_faces[id].normal[2];
    x = dxyz[0] * faux_faces[id].u_dir[0] 
      + dxyz[1] * faux_faces[id].u_dir[1] 
      + dxyz[2] * faux_faces[id].u_dir[2];
    y = dxyz[0] * faux_faces[id].v_dir[0] 
      + dxyz[1] * faux_faces[id].v_dir[1] 
      + dxyz[2] * faux_faces[id].v_dir[2];
    uv[1] = atan2( y, x );
    xyznew[0] = faux_faces[id].center[0] + uv[0]*faux_faces[id].normal[0];
    xyznew[1] = faux_faces[id].center[1] + uv[0]*faux_faces[id].normal[1];
    xyznew[2] = faux_faces[id].center[2] + uv[0]*faux_faces[id].normal[2];
    x = faux_faces[id].offset * cos( uv[1] );
    xyznew[0] += x * faux_faces[id].u_dir[0];
    xyznew[1] += x * faux_faces[id].u_dir[1];
    xyznew[2] += x * faux_faces[id].u_dir[2];
    y = faux_faces[id].offset * sin( uv[1] );
    xyznew[0] += y * faux_faces[id].v_dir[0];
    xyznew[1] += y * faux_faces[id].v_dir[1];
    xyznew[2] += y * faux_faces[id].v_dir[2];
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

  SUPRESS_UNUSED_COMPILER_WARNING(vol);
  SUPRESS_UNUSED_COMPILER_WARNING(edgeId);

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
  double x, y;

  SUPRESS_UNUSED_COMPILER_WARNING(vol);

  if (nfaux_is_uninitialized == nfaux) { 
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
  case general_plane:
    xyz[0] = faux_faces[id].offset * faux_faces[id].normal[0];
    xyz[1] = faux_faces[id].offset * faux_faces[id].normal[1];
    xyz[2] = faux_faces[id].offset * faux_faces[id].normal[2];
    xyz[0] += uv[0] * faux_faces[id].u_dir[0];
    xyz[1] += uv[0] * faux_faces[id].u_dir[1];
    xyz[2] += uv[0] * faux_faces[id].u_dir[2];
    xyz[0] += uv[1] * faux_faces[id].v_dir[0];
    xyz[1] += uv[1] * faux_faces[id].v_dir[1];
    xyz[2] += uv[1] * faux_faces[id].v_dir[2];    
    break;
  case cylinder:
    xyz[0] = faux_faces[id].center[0] + uv[0]*faux_faces[id].normal[0];
    xyz[1] = faux_faces[id].center[1] + uv[0]*faux_faces[id].normal[1];
    xyz[2] = faux_faces[id].center[2] + uv[0]*faux_faces[id].normal[2];
    x = faux_faces[id].offset * cos( uv[1] );
    xyz[0] += x * faux_faces[id].u_dir[0];
    xyz[1] += x * faux_faces[id].u_dir[1];
    xyz[2] += x * faux_faces[id].u_dir[2];
    y = faux_faces[id].offset * sin( uv[1] );
    xyz[0] += y * faux_faces[id].v_dir[0];
    xyz[1] += y * faux_faces[id].v_dir[1];
    xyz[2] += y * faux_faces[id].v_dir[2];
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

  SUPRESS_UNUSED_COMPILER_WARNING(vol);

  if (nfaux_is_uninitialized == nfaux) { 
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

  return TRUE;
}
