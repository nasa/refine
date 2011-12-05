
#ifndef REF_METRIC_H
#define REF_METRIC_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_METRIC_STRUCT REF_METRIC_STRUCT;
typedef REF_METRIC_STRUCT * REF_METRIC;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_METRIC_STRUCT {
  REF_INT max;
  REF_DBL *m;
};

REF_STATUS ref_metric_create( REF_METRIC *ref_metric );
REF_STATUS ref_metric_free( REF_METRIC ref_metric );

#define ref_metric_max( ref_metric ) ((ref_metric)->max)
#define ref_metric_m( ref_metric, entry, node ) ((ref_metric)->m[(entry)+6*(node)])

REF_STATUS ref_metric_set( REF_METRIC ref_metric, REF_INT node, REF_DBL *m );

END_C_DECLORATION

#endif /* REF_METRIC_H */
