
#ifndef REF_SORT_H
#define REF_SORT_H

#include "ref_defs.h"

REF_STATUS ref_sort_insertion( REF_INT n, REF_INT *original, REF_INT *sorted );
REF_STATUS ref_sort_unique( REF_INT n, REF_INT *original, 
			    REF_INT *nunique, REF_INT *unique );

END_C_DECLORATION

#endif /* REF_SORT_H */
