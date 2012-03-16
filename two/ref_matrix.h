
#ifndef REF_MATRIX_H
#define REF_MATRIX_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define ref_matrix_eig( d, e ) ((d)[(e)])
#define ref_matrix_vec( d, xyz, v) ((d)[(xyz)+3+3*(v)])

REF_STATUS ref_matrix_diagonalize( REF_DBL *m_upper_tri, 
				   REF_DBL *diagonal_system );

END_C_DECLORATION

#endif /* REF_MATRIX_H */
