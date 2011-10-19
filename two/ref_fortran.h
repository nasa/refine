
#ifndef REF_FORTRAN_H
#define REF_FORTRAN_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

REF_STATUS ref_init_node__( REF_INT *nnodes, REF_INT *nnodesg,
			    REF_INT *l2g, REF_INT *part, REF_INT *partition, 
			    REF_DBL *x, REF_DBL *y, REF_DBL *z );
REF_STATUS ref_init_node_( REF_INT *nnodes, REF_INT *nnodesg,
			   REF_INT *l2g, REF_INT *part, REF_INT *partition, 
			   REF_DBL *x, REF_DBL *y, REF_DBL *z );

REF_STATUS ref_free_( void );

END_C_DECLORATION

#endif /* REF_FORTRAN_H */
