
#include <stdlib.h>
#include <stdio.h>

#include "ref_matrix.h"

REF_STATUS ref_matrix_diagonalize( REF_DBL *m, 
				   REF_DBL *d )
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

  return REF_SUCCESS;
}
