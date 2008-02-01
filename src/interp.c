
/* Interp, a list items ranked in priorty
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

#include <stdlib.h>
#include <stdio.h>
#ifndef __APPLE__       /* Not needed on Mac OS X */
#include <malloc.h>
#endif
#ifdef __APPLE__       /* Not needed on Mac OS X */
#include <float.h>
#else
#include <values.h>
#endif
#include "interp.h"
#include "gridmath.h"

  /* Rule 4 Solin, Segeth and Dolezel (SSD): order 6 */
static int nq = 24;
static double xq[] = { -0.570794257481696, -0.287617227554912, -0.570794257481696,
		  -0.570794257481696, -0.918652082930777, 0.755956248792332,
		  -0.918652082930777, -0.918652082930777, -0.355324219715449,
		  -0.934027340853653, -0.355324219715449, -0.355324219715449,
		  -0.872677996249965, -0.872677996249965, -0.872677996249965,
		  -0.872677996249965, -0.872677996249965, -0.872677996249965,
		  -0.460655337083368, -0.460655337083368, -0.460655337083368,
		  0.206011329583298, 0.206011329583298, 0.206011329583298 };
static double yq[] = { -0.570794257481696, -0.570794257481696, -0.287617227554912,
		  -0.570794257481696, -0.918652082930777, -0.918652082930777,
		  0.755956248792332, -0.918652082930777, -0.355324219715449,
		  -0.355324219715449, -0.934027340853653, -0.355324219715449,
		  -0.872677996249965, -0.460655337083368, -0.872677996249965,
		  0.206011329583298, -0.460655337083368, 0.206011329583298,
		  -0.872677996249965, -0.872677996249965, 0.206011329583298,
		  -0.872677996249965, -0.872677996249965, -0.460655337083368 };
static double zq[] = { -0.570794257481696, -0.570794257481696, -0.570794257481696,
		  -0.287617227554912, -0.918652082930777, -0.918652082930777,
		  -0.918652082930777, 0.755956248792332, -0.355324219715449,
		  -0.355324219715449, -0.355324219715449, -0.934027340853653,
		  -0.460655337083368, -0.872677996249965, 0.206011329583298,
		  -0.872677996249965, 0.206011329583298, -0.460655337083368,
		  -0.872677996249965, 0.206011329583298, -0.872677996249965,
		  -0.460655337083368, -0.872677996249965, -0.872677996249965 };
static double wq[] = { 0.053230333677557, 0.053230333677557, 0.053230333677557,
		  0.053230333677557, 0.013436281407094, 0.013436281407094,
		  0.013436281407094, 0.013436281407094, 0.073809575391540,
		  0.073809575391540, 0.073809575391540, 0.073809575391540,
		  0.064285714285714, 0.064285714285714, 0.064285714285714,
		  0.064285714285714, 0.064285714285714, 0.064285714285714,
		  0.064285714285714, 0.064285714285714, 0.064285714285714,
		  0.064285714285714, 0.064285714285714, 0.064285714285714 };

static double tet_volume6( double *a, double *b, double *c, double *d )
{
  double m11, m12, m13;
  double det;

  m11 = (a[0]-d[0])*((b[1]-d[1])*(c[2]-d[2])-(c[1]-d[1])*(b[2]-d[2]));
  m12 = (a[1]-d[1])*((b[0]-d[0])*(c[2]-d[2])-(c[0]-d[0])*(b[2]-d[2]));
  m13 = (a[2]-d[2])*((b[0]-d[0])*(c[1]-d[1])-(c[0]-d[0])*(b[1]-d[1]));
  det = ( m11 - m12 + m13 );

  return(-det);
}

Interp* interpCreate( Grid *grid, int function_id, int order, int error_order )
{
  Interp *interp;
  int cell;
  int nodes[4];
  int node;
  double xyz[3];
  double weight;

  interp = (Interp *)malloc( sizeof(Interp) );

  interp->grid = gridDup(grid);
  interp->function_id = function_id;
  interp->error_order = error_order;	
  interp->dim = 1;
  interp->w = NULL;
  if ( order < 0 ) {
    interp->f = NULL;
  }else{
    if ( 1 != order ) 
      {
	printf("interpCreate %d order not implemented\n",order);
	return NULL;
      }
    interp->order = EMPTY; /* for temporary exact interp */	
    interp->f = (double *)malloc( 4*gridNCell(interp->grid)*sizeof(double) );
    for(cell=0;cell<gridNCell(interpGrid(interp));cell++)
      {
	gridCell(interpGrid(interp),cell,nodes);
	for(node=0;node<4;node++)
	  {
	    gridNodeXYZ(grid, nodes[node], xyz );
	    interpFunction( interp, xyz, &(interp->f[node+4*cell]), &weight );
	  }
	if (TRUE) {
	  double a[20];
	  int row, col;
	  int iq,j;
	  double xyz0[3];
	  double xyz1[3];
	  double xyz2[3];
	  double xyz3[3];
	  double bary[4];
	  double func;
	  gridNodeXYZ(grid, nodes[0], xyz0 );
	  gridNodeXYZ(grid, nodes[1], xyz1 );
	  gridNodeXYZ(grid, nodes[2], xyz2 );
	  gridNodeXYZ(grid, nodes[3], xyz3 );
	  for(row=0;row<20;row++)a[row]=0.0;
	  for(iq=0;iq<nq;iq++)
	    {
	      bary[1] = 0.5*(1.0+xq[iq]); 
	      bary[2] = 0.5*(1.0+yq[iq]); 
	      bary[3] = 0.5*(1.0+zq[iq]); 
	      bary[0] = 1.0-bary[1]-bary[2]-bary[3];
	      for(j=0;j<3;j++)
		xyz[j] = bary[0]*xyz0[j] + bary[1]*xyz1[j] + 
		  bary[2]*xyz2[j] + bary[3]*xyz3[j];
	      interpFunction( interp, xyz, &func, &weight );
	      for(row=0;row<4;row++)
		{
		  for(col=0;col<4;col++)
		    a[row+col*4] += bary[row]*bary[col];
		  a[row+4*4] += bary[row]*func;
		}
	    }
	  gridGaussianElimination( 4, 5, a );
	  gridGaussianBacksolve( 4, 5, a );
	  for(row=0;row<4;row++)interp->f[row+4*cell] = a[row+4*4];
	  
	}
      }
  }
  interp->order = order;
  return interp;
}

Interp* interpDirect( Grid *grid )
{
  Interp *interp;
  int i;
  int cell;
  int nodes[4];
  int node, line;
  int nrow, nb;
  int loc2row[10];
  int node0, node1, conn;
  double func[5];
  double *f;
  FILE *file;

  interp = (Interp *)malloc( sizeof(Interp) );

  interp->grid = gridDup(grid);
  interp->function_id = EMPTY;
  interp->error_order = 1;	
  interp->order = 2;	
  interp->dim = 10;

  gridCreateConn(interpGrid(interp));

  nb = interpNB(interp);
  nrow = gridNNode(interpGrid(interp)) + gridNConn(interpGrid(interp));
  f = (double *)malloc( interpDim(interp)*nrow*sizeof(double) );

  interp->f = (double *)malloc( interpDim(interp)*nb*gridNCell(interp->grid)*
				sizeof(double) );

  file = fopen("direct_flow","r");
  for(node=0;node<gridNNode(interpGrid(interp));node++)
    fscanf(file,"%lf %lf %lf %lf %lf",
	   &(f[0+interpDim(interp)*node]),
	   &(f[1+interpDim(interp)*node]),
	   &(f[2+interpDim(interp)*node]),
	   &(f[3+interpDim(interp)*node]),
	   &(f[4+interpDim(interp)*node]) );
  for(line=0;line<gridNConn(interpGrid(interp));line++)
    {
      fscanf(file,"%d %d %lf %lf %lf %lf %lf",&node0,&node1,
	     &(func[0]),
	     &(func[1]),
	     &(func[2]),
	     &(func[3]),
	     &(func[4]) );
      node0--;
      node1--;
      conn = gridFindConn(interpGrid(interp), node0, node1 );
      if ( EMPTY == conn )
	{
	  printf("%s: %d: EMPTY conn %d %d.\n",__FILE__,__LINE__,node0,node1);
	  return NULL;
	}
      for(i=0;i<5;i++)
	f[i+interpDim(interp)*(gridNNode(interpGrid(interp))+conn)] = func[i];
    }
  fclose(file);

  for(cell=0;cell<gridNCell(interpGrid(interp));cell++)
    {
      gridCell(interpGrid(interp),cell,nodes);
      if ( !interpLoc2Row( interp, cell, loc2row )) 
	{
	  printf("%s: %d: interpLoc2Row false.\n",__FILE__,__LINE__);
	  return NULL;
	}
      for(node=0;node<nb;node++)
	for(i=0;i<5;i++)
	  interp->f[i+interpDim(interp)*(node+nb*cell)] = 
	    f[i+interpDim(interp)*loc2row[node]];
    }

  file = fopen("direct_sadj","r");
  for(node=0;node<gridNNode(interpGrid(interp));node++)
    fscanf(file,"%lf %lf %lf %lf %lf",
	   &(f[0+interpDim(interp)*node]),
	   &(f[1+interpDim(interp)*node]),
	   &(f[2+interpDim(interp)*node]),
	   &(f[3+interpDim(interp)*node]),
	   &(f[4+interpDim(interp)*node]) );
  for(line=0;line<gridNConn(interpGrid(interp));line++)
    {
      fscanf(file,"%d %d %lf %lf %lf %lf %lf",&node0,&node1,
	     &(func[0]),
	     &(func[1]),
	     &(func[2]),
	     &(func[3]),
	     &(func[4]) );
      node0--;
      node1--;
      conn = gridFindConn(interpGrid(interp), node0, node1 );
      if ( EMPTY == conn )
	{
	  printf("%s: %d: EMPTY conn %d %d.\n",__FILE__,__LINE__,node0,node1);
	  return NULL;
	}
      for(i=0;i<5;i++)
	f[i+interpDim(interp)*(gridNNode(interpGrid(interp))+conn)] = func[i];
    }
  fclose(file);

  for(cell=0;cell<gridNCell(interpGrid(interp));cell++)
    {
      gridCell(interpGrid(interp),cell,nodes);
      if ( !interpLoc2Row( interp, cell, loc2row )) 
	{
	  printf("%s: %d: interpLoc2Row false.\n",__FILE__,__LINE__);
	  return NULL;
	}
      for(node=0;node<nb;node++)
	for(i=0;i<5;i++)
	  interp->f[5+i+interpDim(interp)*(node+nb*cell)] = 
	    f[i+interpDim(interp)*loc2row[node]];
    }

  interp->w = (double *)malloc( interpDim(interp)*nb*gridNCell(interp->grid)*
				sizeof(double) );

  file = fopen("direct_sres","r");
  for(node=0;node<gridNNode(interpGrid(interp));node++)
    fscanf(file,"%lf %lf %lf %lf %lf",
	   &(f[0+interpDim(interp)*node]),
	   &(f[1+interpDim(interp)*node]),
	   &(f[2+interpDim(interp)*node]),
	   &(f[3+interpDim(interp)*node]),
	   &(f[4+interpDim(interp)*node]) );
  for(line=0;line<gridNConn(interpGrid(interp));line++)
    {
      fscanf(file,"%d %d %lf %lf %lf %lf %lf",&node0,&node1,
	     &(func[0]),
	     &(func[1]),
	     &(func[2]),
	     &(func[3]),
	     &(func[4]) );
      node0--;
      node1--;
      conn = gridFindConn(interpGrid(interp), node0, node1 );
      if ( EMPTY == conn )
	{
	  printf("%s: %d: EMPTY conn %d %d.\n",__FILE__,__LINE__,node0,node1);
	  return NULL;
	}
      for(i=0;i<5;i++)
	f[i+interpDim(interp)*(gridNNode(interpGrid(interp))+conn)] = func[i];
    }
  fclose(file);

  for(cell=0;cell<gridNCell(interpGrid(interp));cell++)
    {
      gridCell(interpGrid(interp),cell,nodes);
      if ( !interpLoc2Row( interp, cell, loc2row )) 
	{
	  printf("%s: %d: interpLoc2Row false.\n",__FILE__,__LINE__);
	  return NULL;
	}
      for(node=0;node<nb;node++)
	for(i=0;i<5;i++)
	  interp->w[i+interpDim(interp)*(node+nb*cell)] = 
	    f[i+interpDim(interp)*loc2row[node]];
    }

  file = fopen("direct_fres","r");
  for(node=0;node<gridNNode(interpGrid(interp));node++)
    fscanf(file,"%lf %lf %lf %lf %lf",
	   &(f[0+interpDim(interp)*node]),
	   &(f[1+interpDim(interp)*node]),
	   &(f[2+interpDim(interp)*node]),
	   &(f[3+interpDim(interp)*node]),
	   &(f[4+interpDim(interp)*node]) );
  for(line=0;line<gridNConn(interpGrid(interp));line++)
    {
      fscanf(file,"%d %d %lf %lf %lf %lf %lf",&node0,&node1,
	     &(func[0]),
	     &(func[1]),
	     &(func[2]),
	     &(func[3]),
	     &(func[4]) );
      node0--;
      node1--;
      conn = gridFindConn(interpGrid(interp), node0, node1 );
      if ( EMPTY == conn )
	{
	  printf("%s: %d: EMPTY conn %d %d.\n",__FILE__,__LINE__,node0,node1);
	  return NULL;
	}
      for(i=0;i<5;i++)
	f[i+interpDim(interp)*(gridNNode(interpGrid(interp))+conn)] = func[i];
    }
  fclose(file);

  for(cell=0;cell<gridNCell(interpGrid(interp));cell++)
    {
      gridCell(interpGrid(interp),cell,nodes);
      if ( !interpLoc2Row( interp, cell, loc2row )) 
	{
	  printf("%s: %d: interpLoc2Row false.\n",__FILE__,__LINE__);
	  return NULL;
	}
      for(node=0;node<nb;node++)
	for(i=0;i<5;i++)
	  interp->w[5+i+interpDim(interp)*(node+nb*cell)] = 
	    f[i+interpDim(interp)*loc2row[node]];
    }

  free(f);

  return interp;
}

void interpFree( Interp *interp )
{
  if ( NULL != interp->w ) free(interp->w);
  if ( NULL != interp->f ) free(interp->f);
  if ( NULL != interp->grid ) gridFree(interp->grid);
  free( interp );
}

int interpNB(Interp *interp)
{
  switch ( interpOrder(interp) ) 
    {
    case 1:
      return 4;
      break;
    case 2:
      return 10;
      break;
    default:
      printf("interpNB %d order not translated to nb\n",interpOrder(interp));
      return EMPTY;
    }
}

GridBool interpPhi( Interp *interp, double *bary, double *phi )
{
  int j;
  double x, y, z;
  x = bary[1];
  y = bary[2];
  z = bary[3];
  switch ( interpOrder(interp) ) 
    {
    case 1:
      for(j=0;j<4;j++) phi[j] = bary[j];
      break;
    case 2:
      phi[0] = 1.0-3.0*z-3.0*y-3.0*x+2.0*z*z+4.0*y*z+4.0*x*z+2.0*y*y+4.0*x*y+2.0*x*x;
      phi[1] = -x+2.0*x*x;
      phi[2] = -y+2.0*y*y;
      phi[3] = -z+2.0*z*z;

      phi[4] = 4.0*x-4.0*x*z-4.0*x*y-4.0*x*x;
      phi[5] = 4.0*y-4.0*y*z-4.0*y*y-4.0*x*y;
      phi[6] = 4.0*z-4.0*z*z-4.0*y*z-4.0*x*z;

      phi[7] = 4.0*x*y;
      phi[8] = 4.0*x*z;
      phi[9] = 4.0*y*z;
      break;
    default:
      printf("interpPhi %d order has no phi\n",interpOrder(interp));
      return FALSE;
    }
  return TRUE;
}

GridBool interpLoc2Row( Interp *interp, int cell, int *loc2row )
{
  int edge;
  Grid *grid = interpGrid(interp);
  switch ( interpOrder(interp) ) 
    {
    case 1:
      if ( grid != gridCell(grid,cell,loc2row)) return FALSE;
      break;
    case 2:
      if ( grid != gridCell(grid,cell,loc2row)) return FALSE;
      for(edge=0;edge<6;edge++)
	loc2row[4+edge] = gridNNode(grid) + gridCell2Conn(grid, cell, edge );
      break;
    default:
      printf("interpLoc2Row %d order has no phi\n",interpOrder(interp));
      return FALSE;
    }
  return TRUE;
}

double interpVectProduct( int nrow, double *v1, double *v2 )
{
  double norm;
  int i;

  norm = 0.0;
  for (i=0;i<nrow;i++) norm += v1[i]*v2[i];

  return norm;
}

GridBool interpReconstructB( Interp *orig, Interp *interp, int nrow, double *b )
{
  int cell, nodes[4];
  int iq,j,row;
  double xyz0[3];
  double xyz1[3];
  double xyz2[3];
  double xyz3[3];
  double xyz[3];
  double volume6;
  double bary[4];
  double func;
  int nb;
  int loc2row[10];
  double phi[10];
  double weight;

  Grid *grid = interpGrid(interp);

  for(row=0;row<nrow;row++) b[row] = 0.0;

  for(cell=0;cell<gridNCell(grid);cell++)
    {
      gridCell(grid,cell,nodes);
      gridNodeXYZ(grid, nodes[0], xyz0 );
      gridNodeXYZ(grid, nodes[1], xyz1 );
      gridNodeXYZ(grid, nodes[2], xyz2 );
      gridNodeXYZ(grid, nodes[3], xyz3 );
      volume6 = tet_volume6(xyz0,xyz1,xyz2,xyz3);
      for(iq=0;iq<nq;iq++)
	{
	  bary[1] = 0.5*(1.0+xq[iq]); 
	  bary[2] = 0.5*(1.0+yq[iq]); 
	  bary[3] = 0.5*(1.0+zq[iq]); 
	  bary[0] = 1.0-bary[1]-bary[2]-bary[3];
	  for(j=0;j<3;j++)
	    xyz[j] = bary[0]*xyz0[j] + bary[1]*xyz1[j] + 
	             bary[2]*xyz2[j] + bary[3]*xyz3[j];
	  interpFunction( orig, xyz, &func, &weight );
	  nb = interpNB(interp);
	  interpPhi( interp, bary, phi );
	  /*	  for(j=0;j<nb;j++) phi[j] *= 0.125 * volume6 * wq[iq]; */
	  if ( !interpLoc2Row( interp, cell, loc2row )) 
	    return FALSE;
	  for(row=0;row<nb;row++)
	    b[loc2row[row]] += func * phi[row];
	}
    }

  return TRUE;  
}

GridBool interpReconstructAx( Interp *interp, int nrow, double *x, double *ax )
{
  int cell, nodes[4];
  int iq,j,row, col;
  double xyz0[3];
  double xyz1[3];
  double xyz2[3];
  double xyz3[3];
  double volume6;
  double bary[4];
  int nb;
  int loc2row[10];
  double phi[10];

  Grid *grid = interpGrid(interp);

  for(row=0;row<nrow;row++) ax[row] = 0.0;

  for(cell=0;cell<gridNCell(grid);cell++)
    {
      gridCell(grid,cell,nodes);
      gridNodeXYZ(grid, nodes[0], xyz0 );
      gridNodeXYZ(grid, nodes[1], xyz1 );
      gridNodeXYZ(grid, nodes[2], xyz2 );
      gridNodeXYZ(grid, nodes[3], xyz3 );
      volume6 = tet_volume6(xyz0,xyz1,xyz2,xyz3);
      for(iq=0;iq<nq;iq++)
	{
	  bary[1] = 0.5*(1.0+xq[iq]); 
	  bary[2] = 0.5*(1.0+yq[iq]); 
	  bary[3] = 0.5*(1.0+zq[iq]); 
	  bary[0] = 1.0-bary[1]-bary[2]-bary[3];
	  nb = interpNB(interp);
	  interpPhi( interp, bary, phi );
	  /*	  for(j=0;j<nb;j++) phi[j] *= 0.125 * volume6 * wq[iq]; */
	  if ( !interpLoc2Row( interp, cell, loc2row ))
	    return FALSE;
	  for(row=0;row<nb;row++)
	    for(col=0;col<nb;col++)
	      ax[loc2row[row]] += phi[row] * phi[col] * x[loc2row[col]];
	}
    }

  return TRUE;  
}

Interp *interpContinuousReconstruction( Interp *orig, 
					int order, int error_order )
{
  Interp *interp;
  int cell;
  int nodes[4];
  int node;
  double xyz[3];
  double *x, *p, *r, *b, *ap;
  int nrow, nb;
  int loc2row[10];
  int iteration;
  double resid0,resid,last_resid,pap,alpha;

  interp = (Interp *)malloc( sizeof(Interp) );

  interp->grid = gridDup(interpGrid(orig));
  interp->function_id = orig->function_id;
  interp->error_order = error_order;
  interp->order = order;

  nb = interpNB(interp);
  switch ( interpOrder(interp) ) 
    {
    case 1:
      nrow = gridNNode(interpGrid(interp));
      break;
    case 2:
      gridCreateConn(interpGrid(interp));
      gridCreateConn(interpGrid(orig));
      nrow = gridNNode(interpGrid(interp)) + gridNConn(interpGrid(interp));
      break;
    default:
      printf( "interpContinuousReconstruction %d order not implemented\n",
	      interpOrder(interp));
      return NULL;
    }
  x  = (double *)malloc( nrow*sizeof(double) );
  p  = (double *)malloc( nrow*sizeof(double) );
  r  = (double *)malloc( nrow*sizeof(double) );
  b  = (double *)malloc( nrow*sizeof(double) );
  ap = (double *)malloc( nrow*sizeof(double) );
  
  interpReconstructB( orig, interp, nrow, b );
  for(node=0;node<nrow;node++) x[node] = 0.0;
  for(node=0;node<nrow;node++) r[node] = b[node];
  for(node=0;node<nrow;node++) p[node] = r[node];

  resid = interpVectProduct( nrow, r, r );
  printf("resid0 %e\n",resid);
  resid0 = resid;
  for (iteration =0; ((iteration<100)&&((resid/resid0)>1.0e-15)); iteration++)
    {
      interpReconstructAx( interp, nrow, p, ap );
      pap = interpVectProduct( nrow, p, ap );
      alpha = resid/pap;
      for(node=0;node<nrow;node++) x[node] += alpha*p[node];
      for(node=0;node<nrow;node++) r[node] -= alpha*ap[node];
      last_resid = resid;
      resid = interpVectProduct( nrow, r, r );
      printf("%3d %e\n",iteration+1,resid/resid0);
      for(node=0;node<nrow;node++) 
	p[node] = r[node] + resid/last_resid*p[node];
    }

  free(p);
  free(r);
  free(b);
  free(ap);

  interp->f = (double *)malloc( nb*gridNCell(interp->grid)*sizeof(double) );
  for(cell=0;cell<gridNCell(interpGrid(interp));cell++)
    {
      gridCell(interpGrid(interp),cell,nodes);
      if ( !interpLoc2Row( interp, cell, loc2row )) return FALSE;
      for(node=0;node<nb;node++)
	  interp->f[node+nb*cell] = x[loc2row[node]];
    }
  free(x);
  return interp;
}

GridBool interpFunction( Interp *interp, double *xyz, 
			 double *func, double *weight )
{
  double a, c, tanhaz;
  double bary[4];
  int cell;
  int i;
  static int last_cell = 0;

  if ( interpOrder(interp) > 0 ) {
    if ( !gridCellValid( interpGrid(interp), last_cell ) )
      for (last_cell = 0; 
	   last_cell < gridMaxCell(interpGrid(interp)); 
	   last_cell++ )
	if ( gridCellValid( interpGrid(interp), last_cell ) ) break;
    cell = gridFindEnclosingCell(interpGrid(interp), last_cell, xyz, bary);
    if ( EMPTY == cell )
      {
	for(i=0;i<interpDim(interp);i++)
	  { 
	    func[i] = 0.0; 
	    if ( NULL != interp->w ) weight[i] = 0.0;
	}
	return FALSE;
      }
    last_cell = cell;
    return interpFunctionInCell( interp, cell, bary, func, weight );
  } else {
    (*weight) = 1.0;
    switch (interpFunctionId(interp)) {
    case 0:
      (*func) = 8.0*xyz[0]*xyz[0] + 32.0*xyz[1]*xyz[1] + 128.0*xyz[2]*xyz[2];
      break;
    case 1:
      a = 10.0;
      c = 1000.0;
      (*func) = 8.0*xyz[0]*xyz[0] + 32.0*xyz[1]*xyz[1] +
	c*tanh(a*(xyz[2]-0.5));
      break;
    }
  }

  return TRUE;
}

GridBool interpFunctionInCell( Interp *interp, 
			       int cell, double *bary, 
			       double *func, double *weight )
{
  int nb, node, i;
  double phi[10];

  for(i=0;i<interpDim(interp);i++)
    func[i] = 0.0;
  for(i=0;i<interpDim(interp);i++)
    weight[i] = 1.0;

  if ( !gridCellValid( interpGrid(interp), cell ) )
    {
      printf("%s: %d: interpFunctionInCell: %s %d\n",
	     __FILE__,__LINE__,"invalid cell",cell);
      return FALSE;      
    }

  nb = interpNB(interp);
  if (EMPTY == nb) {
    printf(" nb %d in interpFunctionInCell\n",nb);
    return FALSE;
  }

  if (!interpPhi( interp, bary, phi ) ) return FALSE;

  for(node=0;node<nb;node++)
    for(i=0;i<interpDim(interp);i++)
      func[i] += interp->f[i+interpDim(interp)*(node+nb*cell)]*phi[node];

  if ( NULL != interp->w )
    {
      for(i=0;i<interpDim(interp);i++)
	weight[i] = 0.0;
      for(node=0;node<nb;node++)
	for(i=0;i<interpDim(interp);i++)
	  weight[i] += interp->w[i+interpDim(interp)*(node+nb*cell)]*phi[node];
    }

  return TRUE;
}

GridBool interpMetric( Interp *interp, double *xyz, double *m )
{
  double a, c, tanhaz;
  switch (interpFunctionId(interp)) {
  case 0:
    m[0] =  16.0;
    m[1] =   0.0;
    m[2] =   0.0;
    m[3] =  64.0;
    m[4] =   0.0;
    m[5] = 256.0;
    break;
  case 1:
    a = 10.0;
    c = 1000.0;
    m[0] =  16.0;
    m[1] =   0.0;
    m[2] =   0.0;
    m[3] =  64.0;
    m[4] =   0.0;
    tanhaz = tanh(a*(xyz[2]-0.5));
    m[5] = -c*a*2.0*tanhaz*(1.0-tanhaz*tanhaz)*a;
    m[5] = MAX(ABS(m[5]),1.0);
    break;
  }
    
  return TRUE;
}

#define TRYIF(worked) {if (!(worked)){ (*error) = DBL_MAX; return TRUE;} } 

GridBool interpError( Interp *interp,
		      double *xyz0, double *xyz1, double *xyz2, double *xyz3, 
		      double *error )
{
  int i,j,k;
  double xyz[3];
  double b0,b1,b2,b3;
  double pinterp[10], func[10], weight[10];
  double phi;
  double n0[10],n1[10],n2[10],n3[10];
  double e01[10],e02[10],e03[10],e12[10],e13[10],e23[10];
  double diff,volume6;

  volume6 = tet_volume6(xyz0,xyz1,xyz2,xyz3);
  if ( volume6 <= 6.0e-12 ) {
    if (ABS(volume6)< 1.0e-15) volume6 = 1.0e-15;
    if ( ABS(volume6) > 0.0 )
      {
	(*error) = 1.0/ABS(volume6);
      }else{
	(*error) = 1.0e40;
      }
    return TRUE; 
  }

  (*error) = 0.0;
  TRYIF( interpFunction( interp, xyz0, n0, weight ) );
  TRYIF( interpFunction( interp, xyz1, n1, weight ) );
  TRYIF( interpFunction( interp, xyz2, n2, weight ) );
  TRYIF( interpFunction( interp, xyz3, n3, weight ) );
  if ( 2 == interpErrorOrder(interp) )
    {
      gridAverageVector(xyz0,xyz1,xyz); 
      TRYIF( interpFunction( interp, xyz, e01, weight ) );
      gridAverageVector(xyz0,xyz2,xyz); 
      TRYIF( interpFunction( interp, xyz, e02, weight ) );
      gridAverageVector(xyz0,xyz3,xyz); 
      TRYIF( interpFunction( interp, xyz, e03, weight ) );
      gridAverageVector(xyz1,xyz2,xyz); 
      TRYIF( interpFunction( interp, xyz, e12, weight ) );
      gridAverageVector(xyz1,xyz3,xyz); 
      TRYIF( interpFunction( interp, xyz, e13, weight ) );
      gridAverageVector(xyz2,xyz3,xyz); 
      TRYIF( interpFunction( interp, xyz, e23, weight ) );
    }

  for(i=0;i<nq;i++)
    for(k=0;k<interpDim(interp);k++)
    {
      b1 = 0.5*(1.0+xq[i]); 
      b2 = 0.5*(1.0+yq[i]); 
      b3 = 0.5*(1.0+zq[i]); 
      b0 = 1.0-b1-b2-b3;
      for(j=0;j<3;j++)
	xyz[j] = b0*xyz0[j] + b1*xyz1[j] + b2*xyz2[j] + b3*xyz3[j];
      TRYIF( interpFunction( interp, xyz, func, weight ) );
      switch ( interpErrorOrder(interp) ) 
	{
	case 1:
	  pinterp[k] = b0*n0[k] + b1*n1[k] + b2*n2[k] + b3*n3[k]; break;
	case 2:
	  pinterp[k] = 0.0;
	  phi = 1.0-3.0*b3-3.0*b2-3.0*b1+2.0*b3*b3+4.0*b2*b3+4.0*b1*b3+2.0*b2*b2+4.0*b1*b2+2.0*b1*b1;
	  pinterp[k] += phi*n0[k];
	  phi = 4.0*b1-4.0*b1*b3-4.0*b1*b2-4.0*b1*b1;
	  pinterp[k] += phi*e01[k];
	  phi = -b1+2.0*b1*b1;
	  pinterp[k] += phi*n1[k];
	  phi = 4.0*b2-4.0*b2*b3-4.0*b2*b2-4.0*b1*b2;
	  pinterp[k] += phi*e02[k];
	  phi = 4.0*b1*b2;
	  pinterp[k] += phi*e12[k];
	  phi = -b2+2.0*b2*b2;
	  pinterp[k] += phi*n2[k];
	  phi = 4.0*b3-4.0*b3*b3-4.0*b2*b3-4.0*b1*b3;
	  pinterp[k] += phi*e03[k];
	  phi = 4.0*b1*b3;
	  pinterp[k] += phi*e13[k];
	  phi = 4.0*b2*b3;
	  pinterp[k] += phi*e23[k];
	  phi = -b3+2.0*b3*b3;
	  pinterp[k] += phi*n3[k];
	  break;
	default:
	  printf("%s: %d: interpError: error_order %d not implemented\n",
		 __FILE__,__LINE__,interpErrorOrder(interp));
	  return FALSE;
	}
      diff = weight[k]*(pinterp[k]-func[k]);
      (*error) += 0.125 * volume6 * wq[i] * diff * diff;
    }

  (*error) = sqrt(*error);

  return TRUE;
}

static int nq1 = 6;
static double tq1[] = {
  0.0337652428984240,
  0.1693953067668678,
  0.3806904069584016,
  0.6193095930415985,
  0.8306046932331322,
  0.9662347571015760 };
static double wq1[] = {
  0.0856622461895852,
  0.1803807865240693,
  0.2339569672863455,
  0.2339569672863455,
  0.1803807865240693,
  0.0856622461895852 };

GridBool interpError1D( Interp *interp,
			double *xyz0, double *xyz1, 
			double *error )
{
  double n0, n1;
  double xyz[3];
  double length;
  int iq, j;
  double t;
  double pinterp, func;
  double diff;
  double weight;

  if ( 1 != interpErrorOrder(interp) )
    {
      printf("%s: %d: interpErroriD: error_order %d not implemented\n",
	     __FILE__,__LINE__,interpErrorOrder(interp));
    }

  (*error) = 0.0;
  interpFunction( interp, xyz0, &n0, &weight );
  interpFunction( interp, xyz1, &n1, &weight );
  gridSubtractVector(xyz1,xyz0,xyz);
  length = gridVectorLength(xyz);
  for(iq=0;iq<nq1;iq++)
    {
      t = tq1[iq];
      for(j=0;j<3;j++)
	xyz[j] = t*xyz1[j] + (1.0-t)*xyz0[j];
      pinterp = t*n1 + (1.0-t)*n0;
      interpFunction( interp, xyz, &func, &weight );
      diff = pinterp-func;
      (*error) += length * wq1[iq] * diff * diff;
    }

  (*error) = sqrt(*error);

  return TRUE;
}

GridBool interpSplitImprovement1D( Interp *interp, 
				   double *xyz0, double *xyz1,
				   double *error_before, double *error_after )
{
  double mid[3];
  double error;

  if ( !interpError1D( interp,xyz0,xyz1,error_before ) ) return FALSE;
  gridAverageVector(xyz0,xyz1,mid);

  if ( !interpError1D( interp,xyz0,mid,&error ) ) return FALSE;
  (*error_after) = error*error;
  if ( !interpError1D( interp,mid,xyz1,&error ) ) return FALSE;
  (*error_after) += error*error;

  (*error_after) = sqrt(*error_after);

  return TRUE;
}

GridBool interpSplitImprovement( Interp *interp, 
				 double *xyz0, double *xyz1,
				 double *xyz2, double *xyz3,
				 double *error_before, double *error_after )
{
  double mid[3];
  double error;

  if ( !interpError( interp,xyz0,xyz1,xyz2,xyz3,error_before ) ) return FALSE;
  gridAverageVector(xyz0,xyz1,mid);

  if ( !interpError( interp,mid,xyz1,xyz2,xyz3,&error ) ) return FALSE;
  (*error_after) = error*error;
  if ( !interpError( interp,xyz0,mid,xyz2,xyz3,&error ) ) return FALSE;
  (*error_after) += error*error;

  (*error_after) = sqrt(*error_after);

  return TRUE;
}

void interpTecLine( Interp *interp,  FILE *f, 
		    double *xyz, double *val, double *weight )
{
  int i;
  fprintf( f, "%25.17e %25.17e %25.17e", xyz[0], xyz[1], xyz[2] );
  for(i=0;i<interpDim(interp);i++)
    fprintf( f, " %25.17e", val[i] );
  for(i=0;i<interpDim(interp);i++)
    fprintf( f, " %25.17e", weight[i] );
  fprintf( f, "\n" );
}

GridBool interpTecplot( Interp *interp, char *filename )
{

  FILE *f;
  int nodes[3], face, faceId;
  int i;
  int cell, cellnodes[4];
  double xyz0[3],xyz1[3],xyz2[3],xyz[3];
  double cell0[3],cell1[3],cell2[3],cell3[3];
  double bary[4];
  double val[10];
  double weight[10];
  Grid *grid = interpGrid(interp);

  if (NULL == filename)
    {
      f = fopen( "interp.t", "w" );
    } else {
      f = fopen( filename, "w" );
    }

  if ( NULL == f ) return FALSE;

  fprintf( f, "title=\"tecplot interp file\"\n" );
  fprintf( f, "variables=\"x\",\"y\",\"z\"" );
  for(i=0;i<interpDim(interp);i++)
    fprintf( f, ",\"s%d\"", i );
  for(i=0;i<interpDim(interp);i++)
    fprintf( f, ",\"w%d\"", i );
  fprintf( f, "\n" );

  fprintf( f, "zone t=surf, i=%d, j=%d, f=fepoint, et=triangle\n",
	   6*gridNFace(grid), 4*gridNFace(grid));

  for (face=0;face<gridMaxFace(grid);face++)
    if ( gridFace(grid, face, nodes, &faceId) )
      if ( interpOrder(interp) > 0 ) {
	cell = gridFindCellWithFace(grid, face );
	if ( EMPTY == cell ) return FALSE;
	gridCell(grid,cell,cellnodes);
	gridNodeXYZ(grid, cellnodes[0], cell0 );
	gridNodeXYZ(grid, cellnodes[1], cell1 );
	gridNodeXYZ(grid, cellnodes[2], cell2 );
	gridNodeXYZ(grid, cellnodes[3], cell3 );

	gridNodeXYZ(grid, nodes[0], xyz0 );
	gridBarycentricCoordinate(cell0, cell1, cell2, cell3, xyz0, bary );
	interpFunctionInCell( interp, cell, bary, val, weight );
	interpTecLine( interp, f, xyz0, val, weight );

	gridNodeXYZ(grid, nodes[1], xyz1 );
	gridBarycentricCoordinate(cell0, cell1, cell2, cell3, xyz1, bary );
	interpFunctionInCell( interp, cell, bary, val, weight );
	interpTecLine( interp, f, xyz1, val, weight );

	gridNodeXYZ(grid, nodes[2], xyz2 );
	gridBarycentricCoordinate(cell0, cell1, cell2, cell3, xyz2, bary );
	interpFunctionInCell( interp, cell, bary, val, weight );
	interpTecLine( interp, f, xyz2, val, weight );

	gridAverageVector(xyz1,xyz2,xyz)
	gridBarycentricCoordinate(cell0, cell1, cell2, cell3, xyz, bary );
	interpFunctionInCell( interp, cell, bary, val, weight );
	interpTecLine( interp, f, xyz, val, weight );

	gridAverageVector(xyz0,xyz2,xyz)
	gridBarycentricCoordinate(cell0, cell1, cell2, cell3, xyz, bary );
	interpFunctionInCell( interp, cell, bary, val, weight );
	interpTecLine( interp, f, xyz, val, weight );

	gridAverageVector(xyz0,xyz1,xyz)
	gridBarycentricCoordinate(cell0, cell1, cell2, cell3, xyz, bary );
	interpFunctionInCell( interp, cell, bary, val, weight );
	interpTecLine( interp, f, xyz, val, weight );

      }else{

	gridNodeXYZ(grid, nodes[0], xyz0 );
	interpFunction( interp, xyz0, val, weight );
	interpTecLine( interp, f, xyz0, val, weight );

	gridNodeXYZ(grid, nodes[1], xyz1 );
	interpFunction( interp, xyz1, val, weight );
	interpTecLine( interp, f, xyz1, val, weight );

	gridNodeXYZ(grid, nodes[2], xyz2 );
	interpFunction( interp, xyz2, val, weight );
	interpTecLine( interp, f, xyz2, val, weight );

	gridAverageVector(xyz1,xyz2,xyz);
	interpFunction( interp, xyz, val, weight );
	interpTecLine( interp, f, xyz, val, weight );

	gridAverageVector(xyz0,xyz2,xyz);
	interpFunction( interp, xyz, val, weight );
	interpTecLine( interp, f, xyz, val, weight );

	gridAverageVector(xyz0,xyz1,xyz);
	interpFunction( interp, xyz, val, weight );
	interpTecLine( interp, f, xyz, val, weight );

      }

  for (face=0;face<gridNFace(grid);face++)
    {
      fprintf( f, "%d %d %d\n",6*face+1,6*face+6,6*face+5);
      fprintf( f, "%d %d %d\n",6*face+2,6*face+4,6*face+6);
      fprintf( f, "%d %d %d\n",6*face+3,6*face+5,6*face+4);
      fprintf( f, "%d %d %d\n",6*face+4,6*face+5,6*face+6);
    }
    
  fclose(f);

  return TRUE;  
}

double interpTotalError(Grid *grid) 
{
  int cell, nodes[4];
  double xyz0[3], xyz1[3], xyz2[3], xyz3[3];
  double total_error, cell_error;
  total_error = 0.0;
  Interp *interp= gridInterp(grid);
  for (cell=0;cell<gridMaxCell(grid);cell++){ 
    if (grid==gridCell(grid, cell, nodes)) { 
      gridNodeXYZ(grid,nodes[0],xyz0);
      gridNodeXYZ(grid,nodes[1],xyz1);
      gridNodeXYZ(grid,nodes[2],xyz2);
      gridNodeXYZ(grid,nodes[3],xyz3);
      interpError( interp,
		   xyz0,xyz1,xyz2,xyz3, 
		   &cell_error );
      total_error += cell_error*cell_error;
    }
  }
  return sqrt(total_error);
}
