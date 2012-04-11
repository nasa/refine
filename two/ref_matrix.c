
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_matrix.h"
#include "ref_math.h"
#include "ref_malloc.h"

REF_STATUS ref_matrix_diag_m( REF_DBL *m, 
			      REF_DBL *d )
{
  REF_DBL L,u,v,s;
  REF_DBL e[3];

  REF_DBL f,tst1,tst2;
  REF_INT i,j,l,mm;
  REF_INT l1, l2;
  REF_DBL h, g, p, r, dl1;
  REF_DBL c, c2, c3, el1;
  REF_INT mml, ii, k;
  REF_DBL s2;

  /* one rotation to make tridiagonal ( zero out m[2] ) */
  /* http://www.geometrictools.com/Documentation/EigenSymmetricNxN.pdf */
  /* Eigen System Solvers for Symmetric Matricies, David Eberly */

  L = sqrt( m[1]*m[1] + m[2]*m[2] );

  if ( ref_math_divisible(m[1],L) && ref_math_divisible(m[2],L) )
    {
      u = m[1] / L;
      v = m[2] / L;
      s = 2.0 * u * m[4] + v * ( m[5] - m[3] );

      ref_matrix_eig( d, 0 ) = m[0];
      ref_matrix_eig( d, 1 ) = m[3] + v*s;
      ref_matrix_eig( d, 2 ) = m[5] - v*s;

      ref_matrix_vec( d, 0, 0 ) = 1.0;
      ref_matrix_vec( d, 1, 0 ) = 0.0;
      ref_matrix_vec( d, 2, 0 ) = 0.0;

      ref_matrix_vec( d, 0, 1 ) = 0.0;
      ref_matrix_vec( d, 1, 1 ) = u;
      ref_matrix_vec( d, 2, 1 ) = v;

      ref_matrix_vec( d, 0, 2 ) = 0.0;
      ref_matrix_vec( d, 1, 2 ) = v;
      ref_matrix_vec( d, 2, 2 ) = -u;

      e[0] = L;
      e[1] = m[4]- u*s;
      e[2] = 0.0;
    }
  else
    {
      ref_matrix_eig( d, 0 ) = m[0];
      ref_matrix_eig( d, 1 ) = m[3];
      ref_matrix_eig( d, 2 ) = m[5];

      ref_matrix_vec( d, 0, 0 ) = 1.0;
      ref_matrix_vec( d, 1, 0 ) = 0.0;
      ref_matrix_vec( d, 2, 0 ) = 0.0;

      ref_matrix_vec( d, 0, 1 ) = 0.0;
      ref_matrix_vec( d, 1, 1 ) = 1.0;
      ref_matrix_vec( d, 2, 1 ) = 0.0;

      ref_matrix_vec( d, 0, 2 ) = 0.0;
      ref_matrix_vec( d, 1, 2 ) = 0.0;
      ref_matrix_vec( d, 2, 2 ) = 1.0;

      e[0] = m[1];
      e[1] = m[4];
      e[2] = 0.0;
    }

  c3 = 0; s2 = 0; /* quiet -Wall used without set compiler warning */

  f = 0.0;
  tst1 = 0.0;
  e[2] = 0.0;

#define gridSign(a,b) (b>=0?ABS(a):-ABS(a))

  for( l = 0; l < 3; l++){ /* row_loop  */
    j = 0;
    h = ABS(d[l]) + ABS(e[l]);
    if (tst1 < h) tst1 = h;
    /* look for small sub-diagonal element */
    for( mm = l; mm<3;mm++) { /*test_for_zero_e */
      tst2 = tst1 + ABS(e[mm]);
      if (ABS(tst2 - tst1) < 1.0e-14 ) break;
      /* e[2] is always zero, so there is no exit through the bottom of loop*/
    }
    if (mm != l) { /* l_not_equal_mm */
      do {
	j = j + 1;
	/* set error -- no convergence to an eigenvalue after 30 iterations */
	if (j > 30 ) {
	  RSS( REF_FAILURE, "not converged" );
        }
	/* form shift */
	l1 = l + 1;
	l2 = l1 + 1;
	g = d[l];
	p = (d[l1] - g) / (2.0 * e[l]);
	r = sqrt(p*p+1.0);
       	d[l] = e[l] / (p + gridSign(r,p));
	d[l1] = e[l] * (p + gridSign(r,p));
	dl1 = d[l1];
	h = g - d[l];
	if (l2<=2) {
	  for (i = l2; i<3; i++) {
	    d[i] = d[i] - h;
	  }
	}
	f = f + h;
	/* ql transformation */
	p = d[mm];
	c = 1.0;
	c2 = c;
	el1 = e[l1];
	s = 0.0;
	mml = mm - l;
	for (ii = 0;ii< mml; ii++ ) {
	  c3 = c2;
	  c2 = c;
	  s2 = s;
	  i = mm - ii - 1;
	  g = c * e[i];
	  h = c * p;
	  r = sqrt(p*p+e[i]*e[i]);
	  e[i+1] = s * r;
	  s = e[i] / r;
	  c = p / r;
	  p = c * d[i] - s * g;
	  d[i+1] = h + s * (c * g + s * d[i]);
	  /* form vector */
	  for (k = 0;k< 3;k++){
	    h = ref_matrix_vec( d, k, i+1);
	    ref_matrix_vec( d, k, i+1) =  
	      s * ref_matrix_vec( d, k, i) + c * h;
	    ref_matrix_vec( d, k, i) =  
	      c * ref_matrix_vec( d, k, i) - s * h;
	  }
	}
	p = -s * s2 * c3 * el1 * e[l] / dl1;
	e[l] = s * p;
	d[l] = c * p;
	tst2 = tst1 + ABS(e[l]);
	if (ABS(tst2 - tst1) < 1.0e-14) break;
      } while (REF_TRUE); /* iterate */
    } /* l_not_equal_mm */
    d[l] = d[l] + f;
  } /* row_loop */

  if (REF_FALSE)
    {
      REF_DBL m2[6];
      REF_DBL tol = -1.0;
      RSS( ref_matrix_form_m( d, m2 ), "reform m" );
      RWDS( m2[0], m[0], tol, "m[0]");
      RWDS( m2[1], m[1], tol, "m[1]");
      RWDS( m2[2], m[2], tol, "m[2]");
      RWDS( m2[3], m[3], tol, "m[3]");
      RWDS( m2[4], m[4], tol, "m[4]");
      RWDS( m2[5], m[5], tol, "m[5]");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_ascending_eig( REF_DBL * d )
{
  REF_DBL temp;
  REF_INT i;

    if ( ref_matrix_eig( d, 1 ) > ref_matrix_eig( d, 0 ) )
      {
	temp = ref_matrix_eig( d, 0 ); 
	ref_matrix_eig( d, 0 ) = ref_matrix_eig( d, 1 ); 
	ref_matrix_eig( d, 1 ) = temp;
	for (i=0;i<3;i++)
	  {
	    temp = ref_matrix_vec( d, i, 0 ); 
	    ref_matrix_vec( d, i, 0 ) = ref_matrix_vec( d, i, 1 ); 
	    ref_matrix_vec( d, i, 1 ) = temp;
	  }
      }

    if ( ref_matrix_eig( d, 2 ) > ref_matrix_eig( d, 0 ) )
      {
	temp = ref_matrix_eig( d, 0 ); 
	ref_matrix_eig( d, 0 ) = ref_matrix_eig( d, 2 ); 
	ref_matrix_eig( d, 2 ) = temp;
	for (i=0;i<3;i++)
	  {
	    temp = ref_matrix_vec( d, i, 0 ); 
	    ref_matrix_vec( d, i, 0 ) = ref_matrix_vec( d, i, 2 ); 
	    ref_matrix_vec( d, i, 2 ) = temp;
	  }
      }

    if ( ref_matrix_eig( d, 2 ) > ref_matrix_eig( d, 1 ) )
      {
	temp = ref_matrix_eig( d, 1 ); 
	ref_matrix_eig( d, 1 ) = ref_matrix_eig( d, 2 ); 
	ref_matrix_eig( d, 2 ) = temp;
	for (i=0;i<3;i++)
	  {
	    temp = ref_matrix_vec( d, i, 1 ); 
	    ref_matrix_vec( d, i, 1 ) = ref_matrix_vec( d, i, 2 ); 
	    ref_matrix_vec( d, i, 2 ) = temp;
	  }
      }

   return REF_SUCCESS;
}

REF_STATUS ref_matrix_form_m( REF_DBL *d,
			      REF_DBL *m )
{
  /* m = d * e * d' */

  /*  d * e
     d[ 3] d[ 6] d[ 9]   d[0]   0    0  
     d[ 4] d[ 7] d[10] *   0   d[1]  0
     d[ 5] d[ 8] d[11]     0    0   d[2]
  */

  /* (d*e) * d'
     d[ 3]*d[0]  d[ 6]*d[1]  d[ 9]*d[2]     d[ 3] d[ 4] d[ 5]
     d[ 4]*d[0]  d[ 7]*d[1]  d[10]*d[2]  *  d[ 6] d[ 7] d[ 8]
     d[ 5]*d[0]  d[ 8]*d[1]  d[11]*d[2]     d[ 9] d[10] d[11]
  */

  m[0] =  d[ 3]*d[0]*d[ 3] + d[ 6]*d[1]*d[ 6] +  d[ 9]*d[2]*d[ 9];
  m[1] =  d[ 3]*d[0]*d[ 4] + d[ 6]*d[1]*d[ 7] +  d[ 9]*d[2]*d[10];
  m[2] =  d[ 3]*d[0]*d[ 5] + d[ 6]*d[1]*d[ 8] +  d[ 9]*d[2]*d[11];
  m[3] =  d[ 4]*d[0]*d[ 4] + d[ 7]*d[1]*d[ 7] +  d[10]*d[2]*d[10];
  m[4] =  d[ 4]*d[0]*d[ 5] + d[ 7]*d[1]*d[ 8] +  d[10]*d[2]*d[11];
  m[5] =  d[ 5]*d[0]*d[ 5] + d[ 8]*d[1]*d[ 8] +  d[11]*d[2]*d[11];

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_inv_m( REF_DBL *m,
			     REF_DBL *inv_m_upper_tri)
{
  REF_DBL det;

  det = ref_matrix_det_m(m);

  inv_m_upper_tri[0] = (m[3]*m[5]-m[4]*m[4]);
  inv_m_upper_tri[1] = (m[2]*m[4]-m[1]*m[5]);
  inv_m_upper_tri[2] = (m[1]*m[4]-m[2]*m[3]);
  inv_m_upper_tri[3] = (m[0]*m[5]-m[2]*m[2]);
  inv_m_upper_tri[4] = (m[2]*m[1]-m[0]*m[4]);
  inv_m_upper_tri[5] = (m[0]*m[3]-m[1]*m[1]);

  if ( !ref_math_divisible( inv_m_upper_tri[0], det ) ||
       !ref_math_divisible( inv_m_upper_tri[1], det ) ||
       !ref_math_divisible( inv_m_upper_tri[2], det ) ||
       !ref_math_divisible( inv_m_upper_tri[3], det ) ||
       !ref_math_divisible( inv_m_upper_tri[4], det ) ||
       !ref_math_divisible( inv_m_upper_tri[5], det ) ) return REF_DIV_ZERO;
       
  inv_m_upper_tri[0] /= det;
  inv_m_upper_tri[1] /= det;
  inv_m_upper_tri[2] /= det;
  inv_m_upper_tri[3] /= det;
  inv_m_upper_tri[4] /= det;
  inv_m_upper_tri[5] /= det;

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_log_m( REF_DBL *m_upper_tri,
			     REF_DBL *log_m_upper_tri)
{
  REF_DBL d[12];

  RSS( ref_matrix_diag_m( m_upper_tri, d ), "diag");

  d[0] = log(d[0]);
  d[1] = log(d[1]);
  d[2] = log(d[2]);

  RSS( ref_matrix_form_m( d, log_m_upper_tri ), "form m");

  return REF_SUCCESS;
}
REF_STATUS ref_matrix_exp_m( REF_DBL *m_upper_tri,
			     REF_DBL *exp_m_upper_tri)
{
  REF_DBL d[12];

  RSS( ref_matrix_diag_m( m_upper_tri, d ), "diag");

  d[0] = exp(d[0]);
  d[1] = exp(d[1]);
  d[2] = exp(d[2]);

  RSS( ref_matrix_form_m( d, exp_m_upper_tri ), "form m");


  return REF_SUCCESS;
}
REF_STATUS ref_matrix_average_m( REF_DBL *m0_upper_tri,
				 REF_DBL *m1_upper_tri,
				 REF_DBL *avg_m_upper_tri)
{
  REF_INT i;

  for( i=0;i<6;i++)
    avg_m_upper_tri[i] = 0.5*(m0_upper_tri[i]+m1_upper_tri[i]);

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_mult_m( REF_DBL *m1, REF_DBL *m2,
			      REF_DBL *product )
{

  /* first col */
  product[0] = m1[0]*m2[0] + m1[1]*m2[1] + m1[2]*m2[2];
  product[1] = m1[1]*m2[0] + m1[3]*m2[1] + m1[4]*m2[2];
  product[2] = m1[2]*m2[0] + m1[4]*m2[1] + m1[5]*m2[2];

  /* mid col */
  product[3] = m1[0]*m2[1] + m1[1]*m2[3] + m1[2]*m2[4];
  product[4] = m1[1]*m2[1] + m1[3]*m2[3] + m1[4]*m2[4];
  product[5] = m1[2]*m2[1] + m1[4]*m2[3] + m1[5]*m2[4];

  /* last col */
  product[6] = m1[0]*m2[2] + m1[1]*m2[4] + m1[2]*m2[5];
  product[7] = m1[1]*m2[2] + m1[3]*m2[4] + m1[4]*m2[5];
  product[8] = m1[2]*m2[2] + m1[4]*m2[4] + m1[5]*m2[5];
  
  return REF_SUCCESS;
}

REF_STATUS ref_matrix_show_ab( REF_INT rows, REF_INT cols, REF_DBL *ab )
{
  REF_INT row, col;
  char format[] = "%10.5f" ;
  for ( row = 0; row < rows; row++)
    {
      for (col = 0; col < cols; col++ )
	{
	  printf(format,ab[row+rows*col]);
	  if ( col < cols-1 ) printf(" ");
	  if ( col == rows-1 ) printf("| ");
	}
      printf("\n");
    }
  return REF_SUCCESS;
}

REF_STATUS ref_matrix_solve_ab( REF_INT rows, REF_INT cols, REF_DBL *ab )
{
  REF_INT row, col;
  REF_INT i, j, k;
  REF_INT pivot_row;
  REF_DBL largest_pivot, pivot;
  REF_DBL temp;
  REF_DBL factor;
  REF_DBL rhs;

  for (col=0; col<rows; col++) 
    {  
      /* find largest pivot */
      pivot_row = col;
      largest_pivot = ABS(ab[pivot_row+rows*col]);    
      for (i=col+1;i<rows;i++) 
	{
	  pivot = ABS(ab[i+rows*col]);
	  if ( pivot > largest_pivot ) 
	    {
	      largest_pivot = pivot;
	      pivot_row = i;
	    }
	}

      /* exchange rows to get the best pivot on the diagonal, 
	 unless it is already there */
      if ( pivot_row != col ) 
	for (j=col;j<cols;j++) 
	  {
	    temp = ab[pivot_row+j*rows];
	    ab[pivot_row+j*rows] = ab[col+j*rows];
	    ab[col+j*rows] = temp;
	  }

      /* normalize pivot row */
      pivot = ab[col+rows*col];
      for (j=col;j<cols;j++) 
	{
	  if ( !ref_math_divisible( ab[col+j*rows], pivot ))return REF_DIV_ZERO;
	  ab[col+j*rows] /= pivot;
	}

      /* elimate sub diagonal terms */
      for (i=col+1;i<rows;i++) 
	{
	  factor = ab[i+col*rows];
	  for (j=col;j<cols;j++) 
	    ab[i+j*rows] -= ab[col+j*rows]*factor;
	  
	} 
    }

  for (col=rows; col<cols; col++)
    for (row=rows-1;row>-1;row--) 
      {
	rhs = ab[row+col*rows];
	for (k = row+1; k<rows; k++)
	  rhs -= ab[row+k*rows]*ab[k+col*rows];
	if ( !ref_math_divisible( rhs, ab[row+row*rows] ) ) return REF_DIV_ZERO;
	ab[row+col*rows] = rhs/ab[row+row*rows];
      }

  return REF_SUCCESS;
}

#define fill_ab(row,n1,n0)						\
  ab[(row)+0*6] =     ((n1)[0]-(n0)[0])*((n1)[0]-(n0)[0]);		\
  ab[(row)+1*6] = 2.0*((n1)[0]-(n0)[0])*((n1)[1]-(n0)[1]);		\
  ab[(row)+2*6] = 2.0*((n1)[0]-(n0)[0])*((n1)[2]-(n0)[2]);		\
  ab[(row)+3*6] =     ((n1)[1]-(n0)[1])*((n1)[1]-(n0)[1]);		\
  ab[(row)+4*6] = 2.0*((n1)[1]-(n0)[1])*((n1)[2]-(n0)[2]);		\
  ab[(row)+5*6] =     ((n1)[2]-(n0)[2])*((n1)[2]-(n0)[2]);

REF_STATUS ref_matrix_imply_m( REF_DBL *m, 
			       REF_DBL *xyz0, REF_DBL *xyz1, 
			       REF_DBL *xyz2, REF_DBL *xyz3 )
{
  REF_DBL ab[42];
  REF_INT i;

  fill_ab(0,xyz1,xyz0);
  fill_ab(1,xyz2,xyz0);
  fill_ab(2,xyz3,xyz0);
  fill_ab(3,xyz2,xyz1);
  fill_ab(4,xyz3,xyz1);
  fill_ab(5,xyz3,xyz2);

  for ( i = 0; i<6; i++ )
    ab[i+6*6] = 1.0;

  RSS( ref_matrix_solve_ab( 6, 7, ab ), "matrix singular" );

  for ( i = 0; i<6; i++ )
    m[i] = ab[i+6*6];

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_show_aqr( REF_INT n, REF_DBL *a, REF_DBL *q, REF_DBL *r)
{
  REF_INT row, col;
  char format[] = "%10.5f" ;
  
  for (row = 0; row<n ; row++ )
    {
      for (col=0;col<n;col++)
	{
	  printf(format,a[row+n*col]);
	  printf(" ");
	  if ( col == n-1 ) printf("= ");
	}
      for (col=0;col<n;col++)
	{
	  printf(format,q[row+n*col]);
	  printf(" ");
	  if ( col == n-1 ) printf("x ");
	}
      for (col=0;col<n;col++)
	{
	  printf(format,r[row+n*col]);
	  if ( col < n-1 ) printf(" ");
	}
      printf("\n");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_qr( REF_INT n, REF_DBL *a, REF_DBL *q, REF_DBL *r )
{
  REF_INT i, j, k;

  for (j = 0; j<n ; j++ )
    for (i=0;i<n;i++)
      q[i+n*j] = a[i+n*j];

  for (j = 0; j<n ; j++ )
    for (i=0;i<n;i++)
      r[i+n*j] = 0.0;

  for (k = 0; k<n ; k++ )
    {
      for (i=0;i<n;i++)
	r[k+n*k] += q[i+n*k]*q[i+n*k];
      r[k+n*k] = sqrt( r[k+n*k] );
      for (i=0;i<n;i++)
	{
	  if ( !ref_math_divisible(  q[i+n*k], r[k+n*k] )) return REF_DIV_ZERO;
	  q[i+n*k] /= r[k+n*k];
	}
      for (j=k+1;j<n;j++)
	{
	  for (i=0;i<n;i++)
	    r[k+n*j] += a[i+n*j]*q[i+n*k];
	  for (i=0;i<n;i++)
	    q[i+n*j] -= r[k+n*j]*q[i+n*k];
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_diag_gen( REF_INT n, REF_DBL *a, 
				REF_DBL *values, REF_DBL *vectors )
{
  REF_DBL *q, *r, *rq, *qq, *ab;
  REF_INT i,j,k,iter;
  REF_DBL max_lower, trace, conv, convm;
  REF_DBL len;
  
  ref_malloc( ab, n*(n+1), REF_DBL );
  ref_malloc( qq, n*n, REF_DBL );
  ref_malloc( rq,  n*n, REF_DBL );
  ref_malloc( q,  n*n, REF_DBL );
  ref_malloc( r,  n*n, REF_DBL );

  for (j = 0; j<n ; j++ )
    for (i=0;i<n;i++)
      rq[i+j*n]=a[i+j*n];

  for (j = 0; j<n ; j++ )
    for (i=0;i<n;i++)
      vectors[i+j*n]=0.0;
  for (i=0;i<n;i++)
    vectors[i+i*n]=1.0;

  iter = 0;
  conv = 1.0;
  while (conv > 1.0e-13)
    {
      iter++;

      RSS( ref_matrix_qr( n, rq, q, r ), "qr");
      ref_matrix_mult_gen( n, r, q, rq );

      for (j = 0; j<n ; j++ )
	for (i=0;i<n;i++)
	  qq[i+j*n]=vectors[i+j*n];
      ref_matrix_mult_gen( n, qq, q, vectors );

      max_lower=0.0;
      for (j = 0; j<n ; j++ )
	for (i=j+1;i<n;i++)
	  max_lower = MAX( max_lower, ABS( rq[i+j*n] ) );
      trace = 0.0;
      for (i=0;i<n;i++)trace+= ABS(rq[i+i*n]);
      conv = max_lower/trace;

      if ( iter > 10000 ) {
	printf("value conv %e used %d\n",conv,iter);
	return REF_FAILURE;
      }

    }

  for (i=0;i<n;i++)
    values[i]=rq[i+i*n];

  
  for (k=0;k<n;k++)
    {
      iter = 0;
      conv = 1.0;
      while (conv > 1.0e-13)
	{
	  iter++;
	  for (j=0;j<n;j++)
	    for (i=0;i<n;i++)
	      ab[i+j*n]=a[i+j*n];
	  for (i=0;i<n;i++)
	    ab[i+i*n] -= 1.00001*values[k];

	  for (i=0;i<n;i++)
	    ab[i+n*n] = vectors[i+k*n];

	  RSS( ref_matrix_solve_ab( n, n+1, ab ), "solve" );

	  len = 0.0;
	  for (i=0;i<n;i++)
	    len += ab[i+n*n]*ab[i+n*n];
	  len = sqrt(len);
	  for (i=0;i<n;i++)
	    {
	      if ( !ref_math_divisible( ab[i+n*n],len )) return REF_DIV_ZERO;
	      ab[i+n*n] = ab[i+n*n]/len;
	    }

	  conv = 0.0;
	  for (i=0;i<n;i++) 
	    conv += (vectors[i+k*n]-ab[i+n*n])*(vectors[i+k*n]-ab[i+n*n]);
	  convm = 0.0;
	  for (i=0;i<n;i++) 
	    convm += (vectors[i+k*n]+ab[i+n*n])*(vectors[i+k*n]+ab[i+n*n]);

	  conv = MIN(conv,convm);

	  for (i=0;i<n;i++)
	    vectors[i+k*n] = ab[i+n*n];
	  
	  if ( iter > 10000 ) {
	    printf("vectr %d conv %e used %d\n",k,conv,iter);	  
	    return REF_FAILURE;
	  }

	}
    }

  ref_free( r );
  ref_free( q );
  ref_free( rq );
  ref_free( qq );
  ref_free( ab );

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_mult_gen( REF_INT n, REF_DBL *a, REF_DBL *b, REF_DBL *r )
{
  REF_INT i, j, k;

  for (j = 0; j<n ; j++ )
    for (i=0;i<n;i++)
      {
	r[i+j*n] = 0.0;
	for (k=0;k<n;k++)
	  r[i+j*n] += a[i+k*n]*b[k+j*n];
      }
  return REF_SUCCESS;
}

REF_STATUS ref_matrix_inv_gen( REF_INT n, REF_DBL *orig, REF_DBL *inv )
{
  REF_INT i, j, k;
  REF_DBL *a;
  REF_DBL pivot, scale;

  ref_malloc( a, n*n, REF_DBL );

  for (j = 0; j<n ; j++ )
    for (i=0;i<n;i++)
      a[i+n*j] = orig[i+n*j];

  for (j = 0; j<n ; j++ )
    for (i=0;i<n;i++)
      inv[i+n*j] = 0.0;
  for (i=0;i<n;i++)
    inv[i+n*i] = 1.0;

  for (j = 0; j<n ; j++ )
    {
      /* scale row so a[j+n*j] is 1.0 */
      pivot = a[j+n*j];
      for (k=0;k<n;k++)
	{
	  if ( !ref_math_divisible( a[j+k*n], pivot )) return REF_DIV_ZERO;
	  a[j+k*n] /= pivot;
	  if ( !ref_math_divisible( inv[j+k*n], pivot )) return REF_DIV_ZERO;
	  inv[j+k*n] /= pivot;
	}
      /* eliminate lower triangle */
      for (i=j+1;i<n;i++)
	{
	  if ( !ref_math_divisible( a[i+j*n],a[j+j*n]  )) return REF_DIV_ZERO;
	  scale = a[i+j*n] / a[j+j*n];
	  for (k=0;k<n;k++)
	    a[i+k*n] -= scale * a[j+k*n];
	  for (k=0;k<n;k++)
	    inv[i+k*n] -= scale * inv[j+k*n];
	}

      /* eliminate upper triangle */
      for (i=0;i<j;i++)
	{
	  if ( !ref_math_divisible( a[i+j*n],a[j+j*n]  )) return REF_DIV_ZERO;
	  scale = a[i+j*n] / a[j+j*n];
	  for (k=0;k<n;k++)
	    a[i+k*n] -= scale * a[j+k*n];
	  for (k=0;k<n;k++)
	    inv[i+k*n] -= scale * inv[j+k*n];
	}

    }

  ref_free( a );

  return REF_SUCCESS;
}

REF_STATUS ref_matrix_transpose_gen( REF_INT n, REF_DBL *a, REF_DBL *at )
{
  REF_INT i, j;

  for (j = 0; j<n ; j++ )
    for (i=0;i<n;i++)
      at[j+n*i] = a[i+n*j];

  return REF_SUCCESS;
}
