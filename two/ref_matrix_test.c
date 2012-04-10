#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_matrix.h"

int main( void )
{

  { /* diag decom, already diag and decending */
    REF_DBL m[6]={ 3.0, 0.0, 0.0, 
                        2.0, 0.0,
                             1.0};
    REF_DBL d[12];

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RWDS(3.0,d[0],-1,"eig 0");
    RWDS(2.0,d[1],-1,"eig 1");
    RWDS(1.0,d[2],-1,"eig 2");

    RWDS(1.0,d[3],-1,"x0");
    RWDS(0.0,d[4],-1,"y0");
    RWDS(0.0,d[5],-1,"z0");

    RWDS(0.0,d[6],-1,"x1");
    RWDS(1.0,d[7],-1,"y1");
    RWDS(0.0,d[8],-1,"z1");

    RWDS(0.0,d[9],-1,"x2");
    RWDS(0.0,d[10],-1,"y2");
    RWDS(1.0,d[11],-1,"z2");
  }

  { /* diag decom, already diag and decending */
    REF_DBL m[6]={ 3.0, 0.0, 0.0, 
                        2.0, 0.0,
                             1.0};
    REF_DBL d[12];

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RWDS(3.0,d[0],-1,"eig 0");
    RWDS(2.0,d[1],-1,"eig 1");
    RWDS(1.0,d[2],-1,"eig 2");

    RWDS(1.0,d[3],-1,"x0");
    RWDS(0.0,d[4],-1,"y0");
    RWDS(0.0,d[5],-1,"z0");

    RWDS(0.0,d[6],-1,"x1");
    RWDS(1.0,d[7],-1,"y1");
    RWDS(0.0,d[8],-1,"z1");

    RWDS(0.0,d[9],-1,"x2");
    RWDS(0.0,d[10],-1,"y2");
    RWDS(1.0,d[11],-1,"z2");
  }

  { /* diag decom, already diag and repeated */
    REF_DBL m[6]={ 10.0,  0.0,  0.0, 
                         10.0,  0.0,
                               10.0};
    REF_DBL d[12];

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RWDS(10.0,d[0],-1,"eig 0");
    RWDS(10.0,d[1],-1,"eig 1");
    RWDS(10.0,d[2],-1,"eig 2");

    RWDS(1.0,d[3],-1,"x0");
    RWDS(0.0,d[4],-1,"y0");
    RWDS(0.0,d[5],-1,"z0");

    RWDS(0.0,d[6],-1,"x1");
    RWDS(1.0,d[7],-1,"y1");
    RWDS(0.0,d[8],-1,"z1");

    RWDS(0.0,d[9],-1,"x2");
    RWDS(0.0,d[10],-1,"y2");
    RWDS(1.0,d[11],-1,"z2");
  }

  { /* diag decom, already tri diag */
    REF_DBL m[6]={ 13.0, -4.0,  0.0, 
                          7.0,  0.0,
                                1.0};
    REF_DBL d[12];

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RWDS(15.0,d[0],-1,"eig 0");
    RWDS( 5.0,d[1],-1,"eig 1");
    RWDS( 1.0,d[2],-1,"eig 2");

  }

  { /* diag decom, already tri diag */
    REF_DBL m[6]={ 13.0,  0.0, -4.0, 
                          4.0,  0.0,
                                7.0};
    REF_DBL d[12];

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RWDS(15.0,d[0],-1,"eig 0");
    RWDS( 5.0,d[1],-1,"eig 1");
    RWDS( 4.0,d[2],-1,"eig 2");
  }

  { /* diag decom, already diad need sort */
    REF_DBL m[6]={ 1.0,     0.0,     0.0, 
                         1000.0,     0.0,
                               1000000.0};
    REF_DBL d[12];

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RSS( ref_matrix_ascending_eig( d ), "ascend");

    RWDS(1000000.0,d[0],-1,"eig 0");
    RWDS(1000.0,d[1],-1,"eig 1");
    RWDS(1.0,d[2],-1,"eig 2");

    RWDS(0.0,d[3],-1,"x0");
    RWDS(0.0,d[4],-1,"y0");
    RWDS(1.0,d[5],-1,"z0");

    RWDS(0.0,d[6],-1,"x1");
    RWDS(1.0,d[7],-1,"y1");
    RWDS(0.0,d[8],-1,"z1");

    RWDS(1.0,d[9],-1,"x2");
    RWDS(0.0,d[10],-1,"y2");
    RWDS(0.0,d[11],-1,"z2");
  }

  { /* diag decom, exercise self check */
    REF_DBL m[6]={ 1.0,     2.0,     3.0, 
                            4.0,     5.0,
                                     6.0};
    REF_DBL d[12];
    REF_DBL m2[6];
    REF_DBL tol = -1.0;

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RSS( ref_matrix_form_m( d, m2 ), "reform m" );
    RWDS( m[0], m2[0], tol, "m[0]");
    RWDS( m[1], m2[1], tol, "m[1]");
    RWDS( m[2], m2[2], tol, "m[2]");
    RWDS( m[3], m2[3], tol, "m[3]");
    RWDS( m[4], m2[4], tol, "m[4]");
    RWDS( m[5], m2[5], tol, "m[5]");
  }

  { /* diag decom, exercise self check r1 */
    REF_DBL m[6]={ 1345234.0,    245.0,    1700.0, 
                            45.0,     5.0,
                                 24000.0};
    REF_DBL d[12];
    REF_DBL m2[6];
    REF_DBL tol = 1.0e-9;

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RSS( ref_matrix_form_m( d, m2 ), "reform m" );
    RWDS( m[0], m2[0], tol, "m[0]");
    RWDS( m[1], m2[1], tol, "m[1]");
    RWDS( m[2], m2[2], tol, "m[2]");
    RWDS( m[3], m2[3], tol, "m[3]");
    RWDS( m[4], m2[4], tol, "m[4]");
    RWDS( m[5], m2[5], tol, "m[5]");
  }

  { /* diag decom, exercise self check r2 */
    REF_DBL m[6]={ 1345234.0,  -10000.0,    3400.0, 
		   2345.0,     -15.0,
		   24.0};
    REF_DBL d[12];
    REF_DBL m2[6];
    REF_DBL tol = 1.0e-9;

    RSS( ref_matrix_diagonalize( m, d ), "diag");

    RSS( ref_matrix_form_m( d, m2 ), "reform m" );
    RWDS( m[0], m2[0], tol, "m[0]");
    RWDS( m[1], m2[1], tol, "m[1]");
    RWDS( m[2], m2[2], tol, "m[2]");
    RWDS( m[3], m2[3], tol, "m[3]");
    RWDS( m[4], m2[4], tol, "m[4]");
    RWDS( m[5], m2[5], tol, "m[5]");
  }

  { /* det */
    REF_DBL tol = -1.0;
    REF_DBL m[6]={10.0, 0.0, 0.0, 
		        2.0, 0.0,
                             5.0};
    
    RWDS( 100.0, ref_matrix_det_m(m), tol, "det");
  }

  { /* inv diag */
    REF_DBL tol = -1.0;
    REF_DBL m[6]={10.0, 0.0, 0.0, 
		        2.0, 0.0,
                             5.0};
    
    REF_DBL inv[6];

    RSS( ref_matrix_inv_m( m, inv ), "inv");

    RWDS(   0.1, inv[0], tol, "inv[0]");
    RWDS(   0.0, inv[1], tol, "inv[1]");
    RWDS(   0.0, inv[2], tol, "inv[2]");
    RWDS(   0.5, inv[3], tol, "inv[3]");
    RWDS(   0.0, inv[4], tol, "inv[4]");
    RWDS(   0.2, inv[5], tol, "inv[5]");
  }

  { /* solve x = 1*/
    REF_DBL tol = -1.0;
    REF_INT rows = 1, cols = 2;
    REF_DBL ab[2]={ 1.0, 1.0 };

    RSS( ref_matrix_solve_ab( rows, cols, ab ), "solve");

    RWDS( 1.0, ab[0+0*cols+rows*rows], tol, "x[0]");
  }

  { /* solve 0.5x = 1*/
    REF_DBL tol = -1.0;
    REF_INT rows = 1, cols = 2;
    REF_DBL ab[2]={ 0.5, 1.0 };

    RSS( ref_matrix_solve_ab( rows, cols, ab ), "solve");

    RWDS( 2.0, ab[0+0*cols+rows*rows], tol, "x[0]");
  }

  { /* solve Ix = [1,2]^t*/
    REF_DBL tol = -1.0;
    REF_INT rows = 2, cols = 3;
    REF_DBL ab[6]={ 1.0, 0.0, 0.0, 1.0, 1.0, 2.0};

    RSS( ref_matrix_solve_ab( rows, cols, ab ), "solve");

    RWDS( 1.0, ab[0+0*cols+rows*rows], tol, "x[0]");
    RWDS( 2.0, ab[1+0*cols+rows*rows], tol, "x[1]");
  }

  { /* solve flip(I)x = [1,2]^t*/
    REF_DBL tol = -1.0;
    REF_INT rows = 2, cols = 3;
    REF_DBL ab[6]={ 0.0, 1.0, 1.0, 0.0, 1.0, 2.0};

    RSS( ref_matrix_solve_ab( rows, cols, ab ), "solve");

    RWDS( 2.0, ab[0+0*cols+rows*rows], tol, "x[0]");
    RWDS( 1.0, ab[1+0*cols+rows*rows], tol, "x[1]");
  }

  { /* solve singular 1 */
    REF_INT rows = 1, cols = 2;
    REF_DBL ab[6]={ 0.0, 1.0};

    REIS( REF_DIV_ZERO, ref_matrix_solve_ab( rows, cols, ab ), "expect sing");
  }

  { /* solve singular 2 */
    REF_INT rows = 2, cols = 3;
    REF_DBL ab[6]={ 1.0, 0.0, 1.0, 0.0, 1.0, 2.0};

    REIS( REF_DIV_ZERO, ref_matrix_solve_ab( rows, cols, ab ), "expect sing");
  }

#define sqrt3 (1.73205080756888)
#define sqrt6 (2.44948974278318)

  { /* imply m from from iso tet nodes */
    REF_DBL tol = -1.0;
    REF_DBL m[6];
    REF_DBL xyz0[]={ 1.0/3.0*sqrt3,0.0,0.0};
    REF_DBL xyz1[]={-1.0/6.0*sqrt3,0.5,0.0};
    REF_DBL xyz2[]={-1.0/6.0*sqrt3,-0.5,0.0};
    REF_DBL xyz3[]={0.0,0.0, 1.0/3.0*sqrt6};

    RSS( ref_matrix_imply_m( m, xyz0, xyz1, xyz2, xyz3), "imply");

    RWDS( 1.0, m[0], tol, "m[0]");
    RWDS( 0.0, m[1], tol, "m[1]");
    RWDS( 0.0, m[2], tol, "m[2]");
    RWDS( 1.0, m[3], tol, "m[3]");
    RWDS( 0.0, m[4], tol, "m[4]");
    RWDS( 1.0, m[5], tol, "m[5]");
  }

  { /* imply m from from short z iso tet nodes */
    REF_DBL tol = -1.0;
    REF_DBL m[6];
    REF_DBL xyz0[]={ 1.0/3.0*sqrt3,0.0,0.0};
    REF_DBL xyz1[]={-1.0/6.0*sqrt3,0.5,0.0};
    REF_DBL xyz2[]={-1.0/6.0*sqrt3,-0.5,0.0};
    REF_DBL xyz3[]={0.0,0.0, 0.1 * 1.0/3.0*sqrt6};

    RSS( ref_matrix_imply_m( m, xyz0, xyz1, xyz2, xyz3), "imply");

    RWDS(   1.0, m[0], tol, "m[0]");
    RWDS(   0.0, m[1], tol, "m[1]");
    RWDS(   0.0, m[2], tol, "m[2]");
    RWDS(   1.0, m[3], tol, "m[3]");
    RWDS(   0.0, m[4], tol, "m[4]");
    RWDS( 100.0, m[5], tol, "m[5]");
  }

  { /* imply m from from right tet nodes */
    REF_DBL tol = -1.0;
    REF_DBL m[6];
    REF_DBL xyz0[]={0.0,0.0,0.0};
    REF_DBL xyz1[]={1.0,0.0,0.0};
    REF_DBL xyz2[]={0.0,1.0,0.0};
    REF_DBL xyz3[]={0.0,0.0,1.0};

    RSS( ref_matrix_imply_m( m, xyz0, xyz1, xyz2, xyz3), "imply");

    RWDS( 1.0, m[0], tol, "m[0]");
    RWDS( 0.5, m[1], tol, "m[1]");
    RWDS( 0.5, m[2], tol, "m[2]");
    RWDS( 1.0, m[3], tol, "m[3]");
    RWDS( 0.5, m[4], tol, "m[4]");
    RWDS( 1.0, m[5], tol, "m[5]");
  }

  { /* imply m from from short z right tet nodes */
    REF_DBL tol = -1.0;
    REF_DBL m[6];
    REF_DBL xyz0[]={0.0,0.0,0.0};
    REF_DBL xyz1[]={1.0,0.0,0.0};
    REF_DBL xyz2[]={0.0,1.0,0.0};
    REF_DBL xyz3[]={0.0,0.0,0.1};

    RSS( ref_matrix_imply_m( m, xyz0, xyz1, xyz2, xyz3), "imply");

    RWDS(   1.0, m[0], tol, "m[0]");
    RWDS(   0.5, m[1], tol, "m[1]");
    RWDS(   5.0, m[2], tol, "m[2]");
    RWDS(   1.0, m[3], tol, "m[3]");
    RWDS(   5.0, m[4], tol, "m[4]");
    RWDS( 100.0, m[5], tol, "m[5]");
  }

  { /* qr wiki 3 */
    REF_DBL tol = -1.0;
    REF_INT n=3;
    REF_DBL a[9]= {12.0, 6.0, -4.0,  -51.0, 167.0, 24.0,  4.0, -68.0, -41.0};
    REF_DBL q[9], r[9];

    /*
A=[12.0, 6.0, -4.0;  
   -51.0, 167.0, 24.0;  
    4.0, -68.0, -41.0]' , [Q, R] = qr(A);
Q=-Q
R=-R
[6.0/7.0, 3.0/7.0, 
 -2.0/7.0; -69.0/175.0, 158./175.0, 6.0/35.0; 
 -58.0/175.0, 6.0/175.0, -33/35.0]'-Q
     */

    RSS( ref_matrix_qr( n, a, q, r ), "qr");

    RWDS(   6.0/7.0,   q[0+0*3], tol, "q[0,0]");
    RWDS(   3.0/7.0,   q[1+0*3], tol, "q[1,0]");
    RWDS(  -2.0/7.0,   q[2+0*3], tol, "q[2,0]");
    RWDS( -69.0/175.0, q[0+1*3], tol, "q[0,1]");
    RWDS( 158.0/175.0, q[1+1*3], tol, "q[1,1]");
    RWDS(   6.0/35.0,  q[2+1*3], tol, "q[2,1]");
    RWDS( -58.0/175.0, q[0+2*3], tol, "q[0,2]");
    RWDS(   6.0/175.0, q[1+2*3], tol, "q[1,2]");
    RWDS( -33.0/35.0,  q[2+2*3], tol, "q[2,2]");

    RWDS(  14.0, r[0+0*3], tol, "r[0,0]");
    RWDS(   0.0, r[1+0*3], tol, "r[1,0]");
    RWDS(   0.0, r[2+0*3], tol, "r[2,0]");
    RWDS(  21.0, r[0+1*3], tol, "r[0,1]");
    RWDS( 175.0, r[1+1*3], tol, "r[1,1]");
    RWDS(   0.0, r[2+1*3], tol, "r[2,1]");
    RWDS( -14.0, r[0+2*3], tol, "r[0,2]");
    RWDS( -70.0, r[1+2*3], tol, "r[1,2]");
    RWDS(  35.0, r[2+2*3], tol, "r[2,2]");

  }

  { /* qr I2 */
    REF_DBL tol = -1.0;
    REF_INT n=2;
    REF_DBL a[4]= {1.0, 0.0, 0.0, 1.0};
    REF_DBL q[4], r[4];


    RSS( ref_matrix_qr( n, a, q, r ), "qr");

    RWDS(   1.0, q[0], tol, "q[0]");
    RWDS(   0.0, q[1], tol, "q[1]");
    RWDS(   0.0, q[2], tol, "q[2]");
    RWDS(   1.0, q[3], tol, "q[3]");

    RWDS(   1.0, r[0], tol, "r[0]");
    RWDS(   0.0, r[1], tol, "r[1]");
    RWDS(   0.0, r[2], tol, "r[2]");
    RWDS(   1.0, r[3], tol, "r[3]");

  }

  { /* diag gen */
    REF_DBL tol = 1.0e-5;
    REF_INT n=2;
    REF_DBL a[4]= {3.5,-0.5,0.21429,0.35714};
    REF_DBL vectors[4], values[2];

    /*
a1 = [ 3 1 ;
       1 5 ]
a2 = [ 10 1;
        1 2 ]
a3 = inv(a1)*a2
[vector3, value3] = eig(a3)
     */

    RSS( ref_matrix_gen_diag( n, a, values, vectors ), "gen diag");

    RWDS( 3.46553, values[0], tol, "val[0]");
    RWDS( 0.39161, values[1], tol, "val[1]");

    RWDS( -0.987300, vectors[0+0*2], tol, "vec[0,0]");
    RWDS(  0.158814, vectors[1+0*2], tol, "vec[1,0]");
    RWDS( -0.068775, vectors[0+1*2], tol, "vec[0,1]");
    RWDS(  0.997632, vectors[1+1*2], tol, "vec[1,1]");

  }

  return 0;
}
