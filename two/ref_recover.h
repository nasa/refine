
#ifndef REF_RECOVER_H
#define REF_RECOVER_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_RECOVER_STRUCT REF_RECOVER_STRUCT;
typedef REF_RECOVER_STRUCT * REF_RECOVER;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_RECOVER_STRUCT {
  REF_INT n;
};

REF_STATUS ref_recover_create( REF_RECOVER *ref_recover );
REF_STATUS ref_recover_free( REF_RECOVER ref_recover );

#define ref_recover_n( ref_recover ) ((ref_recover)->n)

END_C_DECLORATION

#endif /* REF_RECOVER_H */
