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

    RWDS(0.0,d[6],-1,"x0");
    RWDS(1.0,d[7],-1,"y0");
    RWDS(0.0,d[8],-1,"z0");

    RWDS(0.0,d[9],-1,"x0");
    RWDS(0.0,d[10],-1,"y0");
    RWDS(1.0,d[11],-1,"z0");
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

    RWDS(0.0,d[6],-1,"x0");
    RWDS(1.0,d[7],-1,"y0");
    RWDS(0.0,d[8],-1,"z0");

    RWDS(0.0,d[9],-1,"x0");
    RWDS(0.0,d[10],-1,"y0");
    RWDS(1.0,d[11],-1,"z0");
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

    RWDS(0.0,d[6],-1,"x0");
    RWDS(1.0,d[7],-1,"y0");
    RWDS(0.0,d[8],-1,"z0");

    RWDS(0.0,d[9],-1,"x0");
    RWDS(0.0,d[10],-1,"y0");
    RWDS(1.0,d[11],-1,"z0");
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

  return 0;
}
