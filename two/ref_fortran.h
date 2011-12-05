
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

REF_STATUS ref_import_cell__(REF_INT *node_per_cell, REF_INT *ncell,
			     REF_INT *c2n);
REF_STATUS ref_import_cell_(REF_INT *node_per_cell, REF_INT *ncell,
			    REF_INT *c2n);

REF_STATUS ref_import_boundary_flag__(REF_INT *boundary_index, 
				      REF_INT *boundary_flag);
REF_STATUS ref_import_boundary_flag_(REF_INT *boundary_index, 
				     REF_INT *boundary_flag);

REF_STATUS ref_import_boundary__(REF_INT *node_per_face, REF_INT *nface,
				 REF_INT *f2n, REF_INT *boundary_index);
REF_STATUS ref_import_boundary_(REF_INT *node_per_face, REF_INT *nface,
				REF_INT *f2n, REF_INT *boundary_index);

REF_STATUS ref_import_metric__(REF_INT *nnodes, REF_DBL *metric);
REF_STATUS ref_import_metric_(REF_INT *nnodes, REF_DBL *metric);

REF_STATUS ref_viz__( void );
REF_STATUS ref_viz_( void );

REF_STATUS ref_free__( void );
REF_STATUS ref_free_( void );

END_C_DECLORATION

#endif /* REF_FORTRAN_H */
