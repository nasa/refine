
#ifndef REF_ADJ_H
#define REF_ADJ_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_ADJ_STRUCT REF_ADJ_STRUCT;
typedef REF_ADJ_STRUCT * REF_ADJ;
typedef struct REF_ADJ_ITEM_STRUCT REF_ADJ_ITEM_STRUCT;
typedef REF_ADJ_ITEM_STRUCT * REF_ADJ_ITEM;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_ADJ_STRUCT {
  REF_INT nnode, nitem;
  REF_INT *first;
  REF_ADJ_ITEM item;
  REF_INT blank;
};

struct REF_ADJ_ITEM_STRUCT {
  REF_INT reference;
  REF_INT next;
};

REF_STATUS ref_adj_create( REF_ADJ *ref_adj );
REF_STATUS ref_adj_free( REF_ADJ ref_adj );

REF_STATUS ref_adj_inspect( REF_ADJ ref_adj );

#define ref_adj_nnode( ref_adj ) ((ref_adj)->nnode)
#define ref_adj_nitem( ref_adj ) ((ref_adj)->nitem)

END_C_DECLORATION

#endif /* REF_ADJ_H */
