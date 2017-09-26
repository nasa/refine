
#ifndef REF_ADAPT_H
#define REF_ADAPT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_ADAPT_STRUCT REF_ADAPT_STRUCT;
typedef REF_ADAPT_STRUCT * REF_ADAPT;
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_ADAPT_STRUCT {
  REF_DBL split_ratio;
  REF_DBL split_quality_absolute;
  REF_DBL split_quality_relative;

  REF_DBL collapse_ratio;
  REF_DBL collapse_quality_absolute;
  REF_DBL collapse_ratio_limit;

  REF_DBL smooth_min_quality;
};

REF_STATUS ref_adapt_create( REF_ADAPT *ref_adapt );
REF_STATUS ref_adapt_deep_copy( REF_ADAPT *ref_adapt_ptr, REF_ADAPT original );
REF_STATUS ref_adapt_free( REF_ADAPT ref_adapt );

REF_STATUS ref_adapt_pass( REF_GRID ref_grid );
REF_STATUS ref_adapt_threed_pass( REF_GRID ref_grid );
REF_STATUS ref_adapt_twod_pass( REF_GRID ref_grid );

END_C_DECLORATION

#endif /* REF_ADAPT_H */
