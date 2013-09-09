
#ifndef REF_FRONT_H
#define REF_FRONT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_FRONT_STRUCT REF_FRONT_STRUCT;
typedef REF_FRONT_STRUCT * REF_FRONT;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_FRONT_STRUCT {
  REF_INT n, node_per;
};

REF_STATUS ref_front_create( REF_FRONT *ref_front, REF_INT node_per );
REF_STATUS ref_front_free( REF_FRONT ref_front );

#define ref_front_n( ref_front ) ((ref_front)->n)
#define ref_front_node_per( ref_front ) ((ref_front)->node_per)

END_C_DECLORATION

#endif /* REF_FRONT_H */
