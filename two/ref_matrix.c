
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_matrix.h"
#include "ref_math.h"

REF_STATUS ref_matrix_diagonalize( REF_DBL *m, 
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

REF_STATUS ref_matrix_log_m( REF_DBL *m_upper_tri,
			     REF_DBL *log_m_upper_tri)
{
  REF_DBL d[12];

  RSS( ref_matrix_diagonalize( m_upper_tri, d ), "diag");

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

  RSS( ref_matrix_diagonalize( m_upper_tri, d ), "diag");

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
