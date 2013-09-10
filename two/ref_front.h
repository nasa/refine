
#ifndef REF_FRONT_H
#define REF_FRONT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_FRONT_STRUCT REF_FRONT_STRUCT;
typedef REF_FRONT_STRUCT * REF_FRONT;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_FRONT_STRUCT {
  REF_INT node_per;
  REF_INT n;
  REF_INT max;
  REF_INT blank;
  REF_INT *f2n;
};

REF_STATUS ref_front_create( REF_FRONT *ref_front, REF_INT node_per );
REF_STATUS ref_front_free( REF_FRONT ref_front );

#define ref_front_n( ref_front ) ((ref_front)->n)
#define ref_front_node_per( ref_front ) ((ref_front)->node_per)

#define ref_front_f2n(ref_front,node,front) \
  ((ref_front)->f2n[(node)+ref_front_node_per(ref_front)*(front)])

#define ref_front_max( ref_front ) ((ref_front)->max)
#define ref_front_blank( ref_front ) ((ref_front)->blank)

#define ref_front_valid(ref_front,face) \
  ( (face) >=0 && (face) < ref_front_max(ref_front) &&			\
    REF_EMPTY != ref_front_f2n(ref_front,0,face) )

#define each_ref_front_valid_face( ref_front, face )			\
  for ( (face) = 0 ;							\
	(face) < ref_front_max(ref_front);				\
	(face)++ )							\
    if ( ref_front_valid( ref_front, face ) )

REF_STATUS ref_front_insert( REF_FRONT ref_front, REF_INT *nodes );
REF_STATUS ref_front_find( REF_FRONT ref_front, REF_INT *nodes,
			   REF_INT *found_face, REF_BOOL *reversed);

END_C_DECLORATION

#endif /* REF_FRONT_H */
