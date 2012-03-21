
#ifndef REF_SORT_H
#define REF_SORT_H

#include "ref_defs.h"

REF_STATUS ref_sort_insertion_int( REF_INT n, REF_INT *original, 
				   REF_INT *sorted );

REF_STATUS ref_sort_heap_int( REF_INT n, REF_INT *original, 
			      REF_INT *sorted_index );

REF_STATUS ref_sort_unique_int( REF_INT n, REF_INT *original, 
				REF_INT *nunique, REF_INT *unique );

REF_STATUS ref_sort_search( REF_INT n, REF_INT *ascending_list, 
			    REF_INT target, REF_INT *position );

END_C_DECLORATION

#endif /* REF_SORT_H */
