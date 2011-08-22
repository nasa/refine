
#ifndef REF_ADJ_H
#define REF_ADJ_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_ADJ_STRUCT REF_ADJ_STRUCT;
typedef REF_ADJ_STRUCT * REF_ADJ;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_ADJ_STRUCT {
  REF_INT blank;
};

REF_STATUS ref_adj_create( REF_ADJ *ref_adj );
REF_STATUS ref_adj_free( REF_ADJ ref_adj );

REF_STATUS ref_adj_inspect( REF_ADJ ref_adj );

END_C_DECLORATION

#endif /* REF_ADJ_H */
