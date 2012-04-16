
#ifndef REF_FORTRAN_H
#define REF_FORTRAN_H

#include "ref_defs.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#define FC_FUNC_(name,NAME) name ## __
#endif

BEGIN_C_DECLORATION

REF_STATUS FC_FUNC_(ref_init_node,REF_INIT_NODE)( REF_INT *nnodes, REF_INT *nnodesg,
			   REF_INT *l2g, REF_INT *part, REF_INT *partition, 
			   REF_DBL *x, REF_DBL *y, REF_DBL *z );

REF_STATUS FC_FUNC_(ref_import_cell,REF_IMPORT_CELL)(REF_INT *node_per_cell, REF_INT *ncell,
			    REF_INT *c2n);

REF_STATUS FC_FUNC_(ref_import_boundary,REF_IMPORT_BOUNDARY)(REF_INT *node_per_face, REF_INT *nface,
				REF_INT *f2n, REF_INT *boundary_index);

REF_STATUS FC_FUNC_(ref_import_metric,REF_IMPORT_METRIC)(REF_INT *nnodes, REF_DBL *metric);

REF_STATUS FC_FUNC_(ref_viz,REF_VIZ)( void );

REF_STATUS FC_FUNC_(ref_free,REF_FREE)( void );

END_C_DECLORATION

#endif /* REF_FORTRAN_H */
