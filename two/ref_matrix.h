
#ifndef REF_MATRIX_H
#define REF_MATRIX_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define ref_matrix_eig( d, e ) ((d)[(e)])
#define ref_matrix_vec( d, xyz, v) ((d)[(xyz)+3+3*(v)])

#define ref_matrix_sqrt_vt_m_v( m, v ) \
  sqrt( (v)[0] * ( (m)[0]*(v)[0] + (m)[1]*(v)[1] + (m)[2]*(v)[2] ) + \
	(v)[1] * ( (m)[1]*(v)[0] + (m)[3]*(v)[1] + (m)[4]*(v)[2] ) + \
	(v)[2] * ( (m)[2]*(v)[0] + (m)[4]*(v)[1] + (m)[5]*(v)[2] ) )

REF_STATUS ref_matrix_diagonalize( REF_DBL *m_upper_tri, 
				   REF_DBL *diagonal_system );

REF_STATUS ref_matrix_ascending_eig( REF_DBL *diagonal_system );

REF_STATUS ref_matrix_form_m( REF_DBL *diagonal_system,
			      REF_DBL *m_upper_tri);

END_C_DECLORATION

#endif /* REF_MATRIX_H */
