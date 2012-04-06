
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

#define ref_matrix_m_determinate(m) 		\
  ( (m)[0]*(m)[3]*(m)[5] +			\
    (m)[1]*(m)[4]*(m)[2]*2.0 -			\
    (m)[0]*(m)[4]*(m)[4] -			\
    (m)[5]*(m)[1]*(m)[1] -			\
    (m)[3]*(m)[2]*(m)[2] ) 

REF_STATUS ref_matrix_diagonalize( REF_DBL *m_upper_tri, 
				   REF_DBL *diagonal_system );

REF_STATUS ref_matrix_ascending_eig( REF_DBL *diagonal_system );

REF_STATUS ref_matrix_form_m( REF_DBL *diagonal_system,
			      REF_DBL *m_upper_tri);

REF_STATUS ref_matrix_log_m( REF_DBL *m_upper_tri,
			     REF_DBL *log_m_upper_tri);
REF_STATUS ref_matrix_exp_m( REF_DBL *m_upper_tri,
			     REF_DBL *exp_m_upper_tri);
REF_STATUS ref_matrix_average_m( REF_DBL *m0_upper_tri,
				 REF_DBL *m1_upper_tri,
				 REF_DBL *avg_m_upper_tri);

REF_STATUS ref_matrix_show_ab( REF_INT rows, REF_INT cols, REF_DBL *ab );
REF_STATUS ref_matrix_solve_ab( REF_INT rows, REF_INT cols, REF_DBL *ab );

REF_STATUS ref_matrix_imply_m( REF_DBL *m_upper_tri, 
			       REF_DBL *xyz0, REF_DBL *xyz1, 
			       REF_DBL *xyz2, REF_DBL *xyz3 );


END_C_DECLORATION

#endif /* REF_MATRIX_H */
