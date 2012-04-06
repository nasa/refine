
#ifndef REF_METRIC_H
#define REF_METRIC_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_metric_imply_from( REF_DBL *metric, REF_GRID ref_grid );

REF_STATUS ref_metric_imply_non_tet( REF_DBL *metric, REF_GRID ref_grid );

END_C_DECLORATION

#endif /* REF_METRIC_H */
