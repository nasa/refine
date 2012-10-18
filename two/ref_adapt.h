
#ifndef REF_ADAPT_H
#define REF_ADAPT_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

extern REF_DBL ref_adapt_split_ratio;
extern REF_DBL ref_adapt_split_quality;

extern REF_DBL ref_adapt_collapse_ratio;
extern REF_DBL ref_adapt_collapse_quality;
extern REF_DBL ref_adapt_collapse_ratio_limit;

REF_STATUS ref_adapt_pass( REF_GRID ref_grid );
REF_STATUS ref_adapt_twod_pass( REF_GRID ref_grid );

END_C_DECLORATION

#endif /* REF_ADAPT_H */
