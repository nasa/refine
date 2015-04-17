
#ifndef REF_HISTOGRAM_H
#define REF_HISTOGRAM_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_HISTOGRAM_STRUCT REF_HISTOGRAM_STRUCT;
typedef REF_HISTOGRAM_STRUCT * REF_HISTOGRAM;
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION
struct REF_HISTOGRAM_STRUCT {
  REF_INT n;
  REF_DBL max, min;
  REF_DBL log_total;
  REF_INT *bins;
};

REF_STATUS ref_histogram_create( REF_HISTOGRAM *ref_histogram );
REF_STATUS ref_histogram_free( REF_HISTOGRAM ref_histogram );

#define ref_histogram_n( ref_histogram ) ((ref_histogram)->n)
#define ref_histogram_max( ref_histogram ) ((ref_histogram)->max)
#define ref_histogram_min( ref_histogram ) ((ref_histogram)->min)
#define ref_histogram_log_total( ref_histogram ) ((ref_histogram)->log_total)
#define ref_histogram_bin( ref_histogram, i ) ((ref_histogram)->bins[(i)])

#define ref_histogram_to_bin(o) \
  ( ((REF_INT)(6.4*log2((o)))) + ref_histogram_n(ref_histogram)/2 )

#define ref_histogram_to_obs(i) \
  ( pow(2.0,0.15625*((REF_DBL)((i)-ref_histogram_n(ref_histogram)/2) ) ) )

REF_STATUS ref_histogram_add( REF_HISTOGRAM ref_histogram, 
			      REF_DBL observation );

REF_STATUS ref_histogram_ratio( REF_GRID ref_grid );
REF_STATUS ref_histogram_quality( REF_GRID ref_grid );

REF_STATUS ref_histogram_gather( REF_HISTOGRAM ref_histogram );
REF_STATUS ref_histogram_print( REF_HISTOGRAM ref_histogram, 
				char *description );
REF_STATUS ref_histogram_export( REF_HISTOGRAM ref_histogram, 
				 char *description );

REF_STATUS ref_histogram_tec_ratio( REF_GRID ref_grid );

END_C_DECLORATION

#endif /* REF_HISTOGRAM_H */
