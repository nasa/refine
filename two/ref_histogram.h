
#ifndef REF_HISTOGRAM_H
#define REF_HISTOGRAM_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_HISTOGRAM_STRUCT REF_HISTOGRAM_STRUCT;
typedef REF_HISTOGRAM_STRUCT * REF_HISTOGRAM;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_HISTOGRAM_STRUCT {
  REF_INT n;
};

REF_STATUS ref_histogram_create( REF_HISTOGRAM *ref_histogram );
REF_STATUS ref_histogram_free( REF_HISTOGRAM ref_histogram );

#define ref_histogram_n( ref_histogram ) ((ref_histogram)->n)

END_C_DECLORATION

#endif /* REF_HISTOGRAM_H */
